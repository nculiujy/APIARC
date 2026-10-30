#!/usr/bin/env python3
"""deepTools coverage and profile/heatmap generation (module 4).

Compared with the previous version this script:

* consumes the filtered BAMs recorded in ``bam_manifest.tsv`` so coverage and
  peak calling are computed on identical alignments;
* groups samples by ``mark`` + ``condition`` instead of string-stripping
  ``_vs_Input`` from file names, which broke as soon as a real contrast existed;
* emits one matrix per mark/condition, which is what the statistics layer reads.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common.metadata import read_chip_metadata, samples_by_group


def run(cmd: list[str]) -> None:
    print("[deeptools] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def read_manifest(path: str) -> list[dict]:
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in handle if line.strip()]




def main() -> None:
    parser = argparse.ArgumentParser(description="bamCoverage + computeMatrix per mark/condition")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--merged-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--norm", default="BPM", choices=["BPM", "RPKM", "CPM", "none"])
    args = parser.parse_args()

    samples = read_chip_metadata(args.metadata)
    by_group = samples_by_group(samples)
    manifest = {(row["sample"], row["role"]): row for row in read_manifest(args.manifest)}

    for group, members in by_group.items():
        out_dir = os.path.join(args.outdir, group)
        os.makedirs(out_dir, exist_ok=True)

        ip_rows = [manifest[(m.sample, "IP")] for m in members]
        input_ids = sorted({m.input_id for m in members})
        input_rows = [manifest[(i, "Input")] for i in input_ids]
        ip_bams = [row["bam"] for row in ip_rows]
        control_bams = [row["bam"] for row in input_rows]
        # ``name`` carries the group name in the role-aware schema; manifests
        # written before that column existed keep the library id.
        labels = [row.get("name") or m.sample for row, m in zip(ip_rows, members)]
        input_names = [row.get("name") or f"Input_{i}" for row, i in zip(input_rows, input_ids)]
        # Replicates of one Input group share a label; the library id keeps the
        # legend from showing the same name twice.
        repeated = Counter(input_names)
        labels += [
            f"{name} ({library})" if repeated[name] > 1 else name
            for name, library in zip(input_names, input_ids)
        ]

        bigwigs = []
        libraries = [m.sample for m in members] + [f"Input_{i}" for i in input_ids]
        # Bigwigs are keyed by library id: replicates of one group share a label
        # and would otherwise overwrite each other.
        for sample_id, bam in zip(libraries, ip_bams + control_bams):
            bw = os.path.join(out_dir, f"{sample_id}.bw")
            cmd = ["bamCoverage", "-b", bam, "-o", bw, "-p", str(args.threads)]
            if args.norm != "none":
                cmd += ["--normalizeUsing", args.norm]
            run(cmd)
            bigwigs.append(bw)

        region_bed = os.path.join(args.merged_dir, group, f"{group}_peaks_merged.bed")
        if not os.path.exists(region_bed):
            raise FileNotFoundError(
                f"missing merged region set for {group}: {region_bed} "
                "(chipseq_4_merged_peaks must run first)"
            )

        matrix = os.path.join(out_dir, f"matrix_{group}.gz")
        run([
            "computeMatrix", "reference-point",
            "-R", region_bed, "-S", *bigwigs,
            "--referencePoint", "center",
            "-b", "2000", "-a", "2000", "--binSize", "10",
            "-p", str(args.threads),
            "--samplesLabel", *labels,
            "-o", matrix,
        ])

        run([
            "plotProfile", "-m", matrix, "-out", os.path.join(out_dir, "profile.pdf"),
            "--plotTitle", group, "--samplesLabel", *labels,
        ])
        run([
            "plotHeatmap", "-m", matrix, "-out", os.path.join(out_dir, "heatmap.pdf"),
            "--colorMap", "RdYlBu_r", "--zMin", "0",
            "--plotTitle", group, "--samplesLabel", *labels,
        ])


if __name__ == "__main__":
    main()
