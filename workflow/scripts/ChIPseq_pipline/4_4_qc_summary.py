#!/usr/bin/env python3
"""ChIP-seq QC summary and gates (module 4).

Metrics implemented here are the ones that can be derived from the alignments
that the pipeline already produces: mapping statistics, duplication rate
(Picard), FRiP (reads in the mark's own peak set) and replicate agreement
(ENCODE-style Jaccard on peak intervals). Cross-correlation NSC/RSC would need
phantompeakqualtools; it is deliberately not faked here.

The rule that calls this script fails when ``qc.fail_on_threshold`` is true and
a sample misses a gate, so a bad library cannot silently enter the integrated
analysis.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common.metadata import read_chip_metadata, samples_by_group


def run(cmd: list[str], capture: bool = True) -> str:
    result = subprocess.run(cmd, check=False, text=True, capture_output=capture)
    if result.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{result.stderr or ''}")
    return result.stdout or ""


def flagstat_counts(bam: str) -> tuple[int, int]:
    """Return (total, mapped) read counts from samtools flagstat."""
    text = run(["samtools", "flagstat", bam])
    total = mapped = 0
    for line in text.splitlines():
        if "in total" in line:
            total = int(line.split()[0])
        elif re.match(r"^\d+\s+\+\s+\d+\s+mapped", line.strip()):
            mapped = int(line.split()[0])
    return total, mapped


def duplication_rate(metrics_file: str) -> float | None:
    """PERCENT_DUPLICATION from a Picard MarkDuplicates metrics file.

    The file is an htsjdk metrics block: ``#`` comment lines first, then a
    ``LIBRARY``-headed column row, then the values. Blank lines separate blocks
    and can appear before the header, so they must not be parsed as data.
    """
    if not os.path.exists(metrics_file):
        return None
    header: list[str] | None = None
    with open(metrics_file) as handle:
        for line in handle:
            if line.startswith("LIBRARY"):
                header = line.rstrip("\n").split("\t")
                continue
            if header is None or not line.strip() or line.startswith("#"):
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) == len(header):
                return float(values[header.index("PERCENT_DUPLICATION")])
    return None


def frip(bam: str, peaks: list[str]) -> float | None:
    """Fraction of reads in the union of the sample's own peaks."""
    if not peaks:
        return None
    _, mapped = flagstat_counts(bam)
    if mapped == 0:
        return None
    in_peaks = 0
    for peak in peaks:
        if not os.path.exists(peak) or os.path.getsize(peak) == 0:
            continue
        in_peaks += int(run(["samtools", "view", "-c", "-L", peak, bam]).strip() or 0)
    return in_peaks / mapped


def jaccard(a: str, b: str) -> float | None:
    text = run(["bedtools", "jaccard", "-a", a, "-b", b])
    for line in text.splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) >= 3:
            return float(fields[2])
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="ChIP-seq QC summary")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--peakdir", required=True)
    parser.add_argument("--bamdir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--frip-threshold", type=float, default=0.05)
    parser.add_argument("--min-mapped", type=int, default=10_000_000)
    parser.add_argument("--fail-on-threshold", action="store_true")
    args = parser.parse_args()

    samples = read_chip_metadata(args.metadata)
    by_group = samples_by_group(samples)

    with open(args.manifest) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        manifest = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in handle if line.strip()]
    bams = {(row["sample"], row["role"]): row["bam"] for row in manifest}

    rows = ["sample\tmark\tcondition\treplicate\tstatus\ttotal_reads\tmapped_reads\t"
            "duplication\tfrip\treplicate_jaccard\tgate"]
    failures: list[str] = []

    for group, members in by_group.items():
        rep_peaks = [
            os.path.join(args.peakdir, m.sample, f"{m.sample}_peaks.bed") for m in members
        ]
        pair_jaccard = None
        if len(rep_peaks) >= 2 and all(os.path.exists(p) for p in rep_peaks[:2]):
            pair_jaccard = jaccard(rep_peaks[0], rep_peaks[1])

        for member in members:
            bam = bams[(member.sample, "IP")]
            total, mapped = flagstat_counts(bam)
            dup = duplication_rate(
                os.path.join(args.bamdir, member.sample, f"{member.sample}.metricsFile")
            )
            sample_frip = frip(bam, [rep_peaks[members.index(member)]])
            gates = []
            if mapped < args.min_mapped:
                gates.append("mapped_reads")
            if sample_frip is not None and sample_frip < args.frip_threshold:
                gates.append("frip")
            gate = "PASS" if not gates else "FAIL:" + ",".join(gates)
            if gates:
                failures.append(f"{member.sample}({gate})")
            rows.append(
                f"{member.sample}\t{member.mark}\t{member.condition}\t{member.replicate}\t"
                f"{'IP'}\t{total}\t{mapped}\t"
                f"{'' if dup is None else f'{dup:.4f}'}\t"
                f"{'' if sample_frip is None else f'{sample_frip:.4f}'}\t"
                f"{'' if pair_jaccard is None else f'{pair_jaccard:.4f}'}\t{gate}"
            )

    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, "qc_summary.tsv")
    with open(out_path, "w") as handle:
        handle.write("\n".join(rows) + "\n")

    print(f"[qc] wrote {out_path}")
    if failures:
        print("[qc] samples failing gates: " + ", ".join(failures))
        if args.fail_on_threshold:
            raise SystemExit("[qc] QC gates failed and qc.fail_on_threshold is true")


if __name__ == "__main__":
    main()
