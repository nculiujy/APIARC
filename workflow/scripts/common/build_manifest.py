#!/usr/bin/env python3
"""Build the project manifest that records *what was compared*.

The manifest exists because the previous version of this project had no record
of what was compared. A reader could see peak files and a DEG table but not
which samples belonged to which condition, which Input belonged to which IP, or
whether a second condition existed at all.

The manifest is derived from the dataset metadata the pipeline runs on, so it
cannot drift from it: rerun this after editing ``workflow/data/<omics>/<dataset>/metadata.csv``.

Read counts are optional. Without a FASTQ directory the column is written as
"unknown" rather than guessed.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common import workspace
from common.metadata import read_chip_metadata

COLUMNS = [
    "dataset", "species", "modality", "mark", "condition", "replicate",
    "sample", "input_id", "paired_rna_sample", "reads", "design_role",
    "supports_differential", "notes",
]


def read_rna_metadata(path: str) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def count_reads(fastq: str) -> str:
    """Number of reads in a gzip-compressed FASTQ, via line count / 4."""
    if not os.path.exists(fastq):
        return "unknown"
    decompressor = "pigz" if subprocess.run(
        ["bash", "-lc", "command -v pigz"], capture_output=True
    ).returncode == 0 else "gzip"
    result = subprocess.run(
        ["bash", "-lc", f"{decompressor} -dc {fastq!s} | wc -l"],
        capture_output=True, text=True, check=True,
    )
    return str(int(result.stdout.strip()) // 4)


def find_fastq(fastq_dir: str, sample: str) -> str:
    matches = sorted(glob.glob(os.path.join(fastq_dir, f"{sample}*.fastq.gz")))
    return matches[0] if len(matches) == 1 else ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--chip-metadata", default=None)
    parser.add_argument("--rna-metadata", default=None)
    parser.add_argument("--dataset", default="")
    parser.add_argument("--species", default=None)
    parser.add_argument("--chip-fastq-dir", default="")
    parser.add_argument("--rna-fastq-dir", default="")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = {}
    if os.path.exists(args.config):
        import yaml

        with open(args.config) as handle:
            config = yaml.safe_load(handle) or {}

    rna_name = workspace.default_name(config, "RNAseq")
    chip_name = workspace.default_name(config, "ChIPseq")
    joint = workspace.integrated_id(rna_name, chip_name)
    chip_metadata = args.chip_metadata or workspace.metadata_from_config(config, "ChIPseq")
    rna_metadata = args.rna_metadata or workspace.metadata_from_config(config, "RNAseq")
    species = args.species or config.get("species", "mm")
    out = args.out or str(Path("result") / "provenance" / joint / "dataset_manifest.csv")
    dataset = args.dataset or f"{rna_name}__{chip_name}"

    chip = read_chip_metadata(chip_metadata)
    rna = read_rna_metadata(rna_metadata) if os.path.exists(rna_metadata) else []

    conditions_by_mark: dict[str, set] = {}
    for sample in chip:
        conditions_by_mark.setdefault(sample.mark, set()).add(sample.condition)
    differential_marks = {
        mark for mark, conditions in conditions_by_mark.items() if len(conditions) > 1
    }

    rna_conditions = sorted({row.get("group", "") for row in rna if row.get("group")})

    rows = []
    for sample in sorted(chip, key=lambda s: (s.mark, s.condition, s.replicate)):
        supports = sample.mark in differential_marks
        if supports:
            note = "two or more conditions present for this mark"
        else:
            note = (
                "single condition: descriptive mode only, no between-condition "
                "difference is defined for this mark"
            )
        reads = (
            count_reads(find_fastq(args.chip_fastq_dir, sample.sample))
            if args.chip_fastq_dir else "unknown"
        )
        rows.append({
            "dataset": dataset,
            "species": species,
            "modality": "ChIP-seq",
            "mark": sample.mark,
            "condition": sample.condition,
            "replicate": sample.replicate,
            "sample": sample.sample,
            "input_id": sample.input_id,
            "paired_rna_sample": "",
            "reads": reads,
            "design_role": "IP",
            "supports_differential": "yes" if supports else "no",
            "notes": note,
        })

    for row in sorted(rna, key=lambda r: (r.get("group", ""), r.get("sample", ""))):
        reads = (
            count_reads(find_fastq(args.rna_fastq_dir, row.get("sample", "")))
            if args.rna_fastq_dir else "unknown"
        )
        rows.append({
            "dataset": dataset,
            "species": species,
            "modality": "RNA-seq",
            "mark": "",
            "condition": row.get("group", ""),
            "replicate": "",
            "sample": row.get("sample", ""),
            "input_id": "",
            "paired_rna_sample": "",
            "reads": reads,
            "design_role": "expression",
            "supports_differential": "yes" if len(rna_conditions) > 1 else "no",
            "notes": f"RNA conditions: {','.join(rna_conditions) or 'none'}",
        })

    outdir = os.path.dirname(out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[manifest] {len(rows)} rows -> {out}\n"
        f"[manifest] marks with >=2 conditions: "
        f"{', '.join(sorted(differential_marks)) or 'none (descriptive mode)'}"
    )


if __name__ == "__main__":
    main()
