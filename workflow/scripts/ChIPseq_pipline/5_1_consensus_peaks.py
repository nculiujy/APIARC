#!/usr/bin/env python3
"""Build a per-mark consensus peak set with replicate support (module 5).

Differential binding cannot be computed on a peak set that was produced by
pooling replicates: doing so hides replicate variability and makes any test
anti-conservative. Here peaks are merged across all libraries of a mark and the
number of libraries supporting each interval is recorded, so downstream counts
can be restricted to reproducible regions.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common.metadata import comparison_samples, read_chip_metadata


def run(cmd: list[str], stdout_path: str | None = None) -> str:
    handle = open(stdout_path, "w") if stdout_path else subprocess.PIPE
    try:
        result = subprocess.run(cmd, stdout=handle, stderr=subprocess.PIPE, text=True)
    finally:
        if stdout_path:
            handle.close()
    if result.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{result.stderr}")
    return result.stdout or ""


def main() -> None:
    parser = argparse.ArgumentParser(description="consensus peaks with replicate support")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--peakdir", required=True)
    parser.add_argument("--mark", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--treatment", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    samples = comparison_samples(
        read_chip_metadata(args.metadata), args.reference, args.treatment, args.mark
    )
    if not samples:
        raise ValueError(
            f"no samples for {args.mark}: {args.treatment}_vs_{args.reference}"
        )

    peak_files = {
        s.sample: os.path.join(args.peakdir, s.sample, f"{s.sample}_peaks.bed")
        for s in samples
    }
    missing = [p for p in peak_files.values() if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"missing peak files: {missing}")

    os.makedirs(args.outdir, exist_ok=True)

    combined = os.path.join(args.outdir, "peaks_combined.bed")
    with open(combined, "w") as out:
        for sample, path in peak_files.items():
            with open(path) as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{sample}\n")
    sorted_bed = combined + ".sorted"
    run(["bedtools", "sort", "-i", combined], stdout_path=sorted_bed)
    merged = os.path.join(args.outdir, "peaks_merged.bed")
    run(["bedtools", "merge", "-i", sorted_bed], stdout_path=merged)

    support = {s.sample: [] for s in samples}
    for sample, path in peak_files.items():
        counts = run(["bedtools", "intersect", "-a", merged, "-b", path, "-c", "-nonamecheck"])
        support[sample] = [int(line.split("\t")[3]) for line in counts.strip().splitlines()]

    consensus_path = os.path.join(args.outdir, "consensus_peaks.bed")
    support_path = os.path.join(args.outdir, "replicate_support.tsv")

    with open(consensus_path, "w") as bed, open(support_path, "w") as table:
        table.write("region_id\t" + "\t".join(s.sample for s in samples) + "\tn_support\twidth\n")
        for index, line in enumerate(open(merged)):
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            region_id = f"{chrom}:{start}-{end}"
            flags = [support[s.sample][index] for s in samples]
            n_support = sum(1 for value in flags if value > 0)
            bed.write(f"{chrom}\t{start}\t{end}\t{region_id}\t{n_support}\n")
            table.write(
                f"{region_id}\t" + "\t".join(str(1 if v > 0 else 0) for v in flags)
                + f"\t{n_support}\t{int(end) - int(start)}\n"
            )

    os.remove(combined)
    os.remove(sorted_bed)
    print(f"[consensus] {args.mark} {args.treatment}_vs_{args.reference}: "
          f"{sum(1 for _ in open(consensus_path))} regions")


if __name__ == "__main__":
    main()
