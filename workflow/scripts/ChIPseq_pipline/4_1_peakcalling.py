#!/usr/bin/env python3
"""MACS2 peak calling for APIARC (module 4).

Changes with respect to the previous implementation, each of them addressing a
reviewer request:

* ``--keep-dup all`` is no longer hard-coded. Duplicates are removed earlier by
  Picard MarkDuplicates; re-enabling them in MACS2 is what reviewers read as
  "PCR artefacts retained".
* Peaks are called per biological replicate, not on pooled BAMs, so replicate
  support can be quantified downstream.
* Broad marks (H3K27ac, H3K4me1, ...) use ``--broad``; the mode is configurable
  per mark.
* MAPQ filtering and ENCODE blacklist removal happen once, and the filtered BAM
  is recorded in ``bam_manifest.tsv`` so every downstream module (coverage,
  counting) consumes the same alignments.
* A MACS2 failure aborts the run. The previous version silently created empty
  peak files, which made a failed sample look like a successful one.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from statistics import median

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common import annotation
from common.metadata import ChipSample, input_ids, read_chip_metadata, uses_group_names

GENOME_SIZE = {"mm": "mm", "homo": "hs", "hs": "hs", "tair": "1.2e8"}


def log(message: str) -> None:
    print(f"[peakcalling] {message}", flush=True)


def run(cmd: list[str], *, capture: bool = False) -> str:
    log(" ".join(cmd))
    result = subprocess.run(
        cmd, check=False, text=True, capture_output=capture
    )
    if result.returncode != 0:
        stderr = (result.stderr or "").strip() if capture else ""
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(cmd)}\n{stderr}")
    return (result.stdout or "") if capture else ""


def resolve_blacklist(config: dict, info: dict, resource_root: str) -> str | None:
    setting = config.get("chip", {}).get("blacklist", "auto")
    if setting in (None, "", "none", False):
        return None
    if setting != "auto":
        if not os.path.exists(str(setting)):
            raise FileNotFoundError(f"chip.blacklist={setting} does not exist")
        return str(setting)
    candidates = [
        os.path.join(info["chip_dir"], "blacklist.bed"),
        os.path.join(resource_root, "anno", info["assembly"], "blacklist.bed"),
    ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    log(
        f"WARNING: no blacklist found at {candidates[0]}; peak calling will run without "
        "blacklist filtering (record this in the Methods section)"
    )
    return None


def filter_bam(
    bam: str, out_bam: str, mapq: int, blacklist: str | None, threads: int
) -> tuple[str, bool]:
    """MAPQ filtering with samtools and blacklist removal with bedtools."""
    if os.path.exists(out_bam) and os.path.exists(out_bam + ".bai"):
        return out_bam, blacklist is not None
    if mapq <= 0 and blacklist is None:
        return bam, False

    mapq_bam = out_bam + ".mapq.bam"
    run([
        "samtools", "view", "-@", str(threads), "-b", "-q", str(max(mapq, 0)),
        "-o", mapq_bam, bam,
    ])

    if not blacklist:
        os.replace(mapq_bam, out_bam)
    else:
        with open(out_bam, "wb") as handle:
            result = subprocess.run(
                ["bedtools", "intersect", "-v", "-abam", mapq_bam, "-b", blacklist,
                 "-nonamecheck"],
                stdout=handle, stderr=subprocess.PIPE,
            )
        if result.returncode != 0:
            raise RuntimeError(
                "bedtools intersect failed: " + result.stderr.decode()
            )
        os.remove(mapq_bam)

    run(["samtools", "index", "-@", str(threads), out_bam])
    return out_bam, blacklist is not None


def peak_mode_for(config: dict, mark: str) -> str:
    """Narrow or broad for one library.

    The five-column schema has no mark column, so ``mark`` is a group name such
    as ``H3K27ac_veh``; a ``chip.peak_mode`` key matches when the group name is
    that key or starts with ``<key>_``.
    """
    modes = config.get("chip", {}).get("peak_mode", {})
    key = mark
    if key not in modes:
        # Two keys can prefix the same group name (``H3K4`` and ``H3K4me3``);
        # the more specific one is the one the user meant.
        prefixes = [name for name in modes if mark.startswith(f"{name}_")]
        key = max(prefixes, key=len) if prefixes else mark
    mode = str(modes.get(key, modes.get("default", "narrow"))).lower()
    if mode not in {"narrow", "broad"}:
        raise ValueError(f"chip.peak_mode for {mark} must be narrow or broad, got {mode}")
    return mode


def effective_genome(config: dict, info: dict) -> str:
    """MACS2 ``-g`` value: a shortcut for real genomes, an explicit size otherwise.

    ``chip.genome_size`` keeps synthetic references (a few hundred kilobases)
    from being called against the mm10/hg19 effective genome size, which would
    make every P-value meaningless.
    """
    override = config.get("chip", {}).get("genome_size")
    if override not in (None, "", "auto"):
        return str(override)
    species = str(info["species"]).lower()
    return GENOME_SIZE.get(species, species)


def extra_macs2_args(config: dict) -> list[str]:
    """Project-level MACS2 switches, appended verbatim to every callpeak call.

    Synthetic or otherwise tiny references cannot satisfy the fragment model
    (MACS2 pairs plus/minus strand peaks to infer the fragment length), so such
    projects have to state ``--nomodel --extsize`` themselves instead of having
    the pipeline guess.
    """
    extra = config.get("chip", {}).get("macs_extra_args") or []
    if isinstance(extra, str):
        extra = extra.split()
    return [str(item) for item in extra]


def macs2_command(
    sample: ChipSample,
    treatment_bam: str,
    control_bam: str,
    outdir: str,
    config: dict,
    info: dict,
    qvalue: float,
    cutoff_type: str,
) -> tuple[list[str], str]:
    mode = peak_mode_for(config, sample.mark)
    extension = "broadPeak" if mode == "broad" else "narrowPeak"
    cmd = [
        "macs2", "callpeak",
        "-t", treatment_bam,
        "-c", control_bam,
        "-g", effective_genome(config, info),
        "-n", sample.sample,
        "--outdir", outdir,
        "-q" if cutoff_type == "qvalue" else "-p", str(qvalue),
    ]
    if mode == "broad":
        cutoff = config.get("chip", {}).get("broad_cutoff", 0.1)
        cmd += ["--broad", "--broad-cutoff", str(cutoff)]
    if config.get("chip", {}).get("keep_duplicates", False):
        cmd += ["--keep-dup", "all"]
    cmd += extra_macs2_args(config)
    return cmd, extension


def normalise_peak_table(raw_peak: str, out_bed: str) -> tuple[int, int, float]:
    """Write chrom/start/end/peak_id/signal/-log10q and return summary stats."""
    widths, qvals = [], []
    with open(raw_peak) as source, open(out_bed, "w") as target:
        for index, line in enumerate(source, start=1):
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            signal = float(fields[6]) if len(fields) > 6 and fields[6] else 0.0
            qval = float(fields[8]) if len(fields) > 8 and fields[8] else 0.0
            widths.append(end - start)
            qvals.append(qval)
            target.write(
                f"{chrom}\t{start}\t{end}\t{chrom}:{start}-{end}\t{signal:.3f}\t{qval:.3f}\n"
            )
    if not widths:
        raise RuntimeError(f"MACS2 produced no peaks: {raw_peak}")
    return len(widths), sum(widths), median(widths)


def main() -> None:
    parser = argparse.ArgumentParser(description="MACS2 peak calling per replicate")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--bamdir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--assembly", default=None, help="assembly name, e.g. GRCh37")
    parser.add_argument("--species", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--qvalue", type=float, default=0.05)
    parser.add_argument("--cutoff-type", default="qvalue", choices=["qvalue", "pvalue"])
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--resource-root", default="workflow/resources")
    args = parser.parse_args()

    import yaml

    with open(args.config) as handle:
        config = yaml.safe_load(handle)

    samples = read_chip_metadata(args.metadata)
    # ``--assembly`` identifies the reference; ``--species`` is kept so an older
    # command line still resolves instead of failing on an unknown option.
    info = annotation.single(args.assembly or args.species or "GRCm38")
    blacklist = resolve_blacklist(config, info, args.resource_root)
    mapq = int(config.get("chip", {}).get("mapq", 30))
    keep_dups = bool(config.get("chip", {}).get("keep_duplicates", False))

    bam_dir = args.bamdir
    peak_root = args.outdir
    os.makedirs(peak_root, exist_ok=True)

    manifest_rows: list[str] = [
        "sample\trole\tmark\tcondition\treplicate\tbam\tfiltered\tname"
    ]
    summary_rows: list[str] = [
        "sample\tmark\tcondition\treplicate\tmode\tn_peaks\ttotal_bp\tmedian_width\tmean_log10q\tkeep_duplicates\tblacklist"
    ]

    filtered_inputs: dict[str, str] = {}
    # The five-column schema names each group; older schemas are labelled by
    # library id so existing plots keep the names they always had.
    named = uses_group_names(samples)
    input_names: dict[str, str] = {}
    for item in samples:
        input_names.setdefault(item.input_id, item.input_name or item.input_id)
    for input_id in input_ids(samples):
        raw = os.path.join(bam_dir, input_id, "accepted_hits.sorted.unique.bam")
        if not os.path.exists(raw):
            raise FileNotFoundError(f"missing Input BAM: {raw}")
        out_bam = os.path.join(bam_dir, input_id, "accepted_hits.sorted.unique.filt.bam")
        filtered, used_blacklist = filter_bam(raw, out_bam, mapq, blacklist, args.threads)
        filtered_inputs[input_id] = filtered
        name = input_names.get(input_id, input_id) if named else input_id
        manifest_rows.append(
            f"{input_id}\tInput\t-\t-\t0\t{filtered}\t{str(used_blacklist).lower()}\t{name}"
        )

    for sample in samples:
        raw = os.path.join(bam_dir, sample.sample, "accepted_hits.sorted.unique.bam")
        if not os.path.exists(raw):
            raise FileNotFoundError(f"missing IP BAM: {raw}")
        out_bam = os.path.join(bam_dir, sample.sample, "accepted_hits.sorted.unique.filt.bam")
        treatment_bam, used_blacklist = filter_bam(
            raw, out_bam, mapq, blacklist, args.threads
        )
        manifest_rows.append(
            f"{sample.sample}\tIP\t{sample.mark}\t{sample.condition}\t{sample.replicate}"
            f"\t{treatment_bam}\t{str(used_blacklist).lower()}"
            f"\t{sample.sample_name if named else sample.sample}"
        )

        sample_dir = os.path.join(peak_root, sample.sample)
        os.makedirs(sample_dir, exist_ok=True)
        control_bam = filtered_inputs[sample.input_id]

        cmd, extension = macs2_command(
            sample, treatment_bam, control_bam, sample_dir, config,
            info, args.qvalue, args.cutoff_type,
        )
        run(cmd)

        raw_peak = os.path.join(sample_dir, f"{sample.sample}_peaks.{extension}")
        if not os.path.exists(raw_peak):
            raise RuntimeError(f"MACS2 did not create {raw_peak}")
        out_bed = os.path.join(sample_dir, f"{sample.sample}_peaks.bed")
        n_peaks, total_bp, median_width = normalise_peak_table(raw_peak, out_bed)

        qvals = []
        with open(out_bed) as handle:
            for line in handle:
                qvals.append(float(line.split("\t")[5]))
        mode = peak_mode_for(config, sample.mark)
        summary_rows.append(
            f"{sample.sample}\t{sample.mark}\t{sample.condition}\t{sample.replicate}\t"
            f"{mode}\t{n_peaks}\t{total_bp}\t{median_width}\t"
            f"{sum(qvals) / len(qvals):.3f}\t{str(keep_dups).lower()}\t"
            f"{str(used_blacklist).lower()}"
        )

    with open(os.path.join(peak_root, "peak_calling_summary.tsv"), "w") as handle:
        handle.write("\n".join(summary_rows) + "\n")
    with open(os.path.join(peak_root, "bam_manifest.tsv"), "w") as handle:
        handle.write("\n".join(manifest_rows) + "\n")

    log(f"finished: {len(samples)} libraries, {len(filtered_inputs)} inputs")


if __name__ == "__main__":
    main()
