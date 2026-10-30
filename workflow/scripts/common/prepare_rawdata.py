#!/usr/bin/env python3
"""Normalise the raw reads of one dataset into the stage-1 hand-off directory.

The dataset directory is the source of truth: it holds ``metadata.csv`` plus the
FASTQ files themselves. Downstream stages expect a flat directory whose file
names match the ``sample`` column, so this step validates presence, checks that
paired-end mates come in twos, and materialises one relative symlink per read
file. Nothing is copied and nothing outside ``--outdir`` is written.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common import workspace


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--omics", required=True, choices=list(workspace.OMICS))
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--data-root", default=str(workspace.DATA_ROOT))
    parser.add_argument("--allow-missing", action="store_true")
    return parser.parse_args(argv)


def canonical_name(sample: str, read: int, suffix: str, paired: bool) -> str:
    if paired:
        return f"{sample}_{read}{suffix}"
    return f"{sample}{suffix}"


def clear_previous_links(outdir: Path) -> None:
    for entry in outdir.iterdir():
        if not entry.is_symlink():
            continue
        if any(entry.name.endswith(suffix) for suffix in workspace.FASTQ_SUFFIXES):
            entry.unlink()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    data_root = Path(args.data_root)
    dataset_root = workspace.dataset_dir(args.omics, args.dataset, data_root)
    metadata = dataset_root / workspace.META_NAME

    if not metadata.is_file():
        print(f"[ERROR] missing metadata: {metadata}", file=sys.stderr)
        return 2

    rows = workspace.read_metadata_rows(metadata)
    if not rows:
        print(f"[ERROR] no samples declared in {metadata}", file=sys.stderr)
        return 2

    table = workspace.index_fastq(dataset_root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    clear_previous_links(outdir)

    missing: list[str] = []
    unpaired: list[str] = []
    records: list[dict] = []
    linked: set[str] = set()

    for row in rows:
        sample = (row.get("sample") or "").strip()
        if not sample:
            continue
        key = (row.get("fastq_prefix") or "").strip() or sample
        slot = workspace.lookup_fastq(table, key)
        if not slot or "r1" not in slot:
            missing.append(f"{sample} (looked for '{key}')")
            continue
        paired = "r2" in slot
        if paired and slot["r1"].stat().st_size != slot["r2"].stat().st_size:
            print(f"[WARN] {sample}: read mates differ in size", file=sys.stderr)
        links = {}
        for read, field in ((1, "r1"), (2, "r2")):
            source = slot.get(field)
            if source is None:
                continue
            suffix = next(s for s in workspace.FASTQ_SUFFIXES if source.name.endswith(s))
            target = outdir / canonical_name(sample, read, suffix, paired)
            if target.exists() or target.is_symlink():
                target.unlink()
            target.symlink_to(os.path.relpath(source, outdir))
            links[read] = target
        if paired and 2 not in links:
            unpaired.append(sample)
        records.append(
            {
                "sample": sample,
                "source_prefix": key,
                "paired": "yes" if paired else "no",
                "read1": links[1].name,
                "read2": links.get(2).name if paired else "",
                "bytes_read1": links[1].stat().st_size,
                "bytes_read2": links[2].stat().st_size if paired else 0,
            }
        )
        linked.add(sample)

    # A role-aware table lists its Inputs as rows, which the loop above already
    # links. Older files name them only through ``input_sample``/``input_id``,
    # yet peak calling opens them from this stage; link those here.
    for input_id in workspace.declared_inputs(args.omics, rows):
        if input_id in linked:
            continue
        slot = workspace.lookup_fastq(table, input_id)
        if not slot or "r1" not in slot:
            missing.append(f"{input_id} (declared as input_id)")
            continue
        paired = "r2" in slot
        links = {}
        for read, field in ((1, "r1"), (2, "r2")):
            source = slot.get(field)
            if source is None:
                continue
            suffix = next(s for s in workspace.FASTQ_SUFFIXES if source.name.endswith(s))
            target = outdir / canonical_name(input_id, read, suffix, paired)
            if target.exists() or target.is_symlink():
                target.unlink()
            target.symlink_to(os.path.relpath(source, outdir))
            links[read] = target
        records.append(
            {
                "sample": input_id,
                "source_prefix": input_id,
                "paired": "yes" if paired else "no",
                "read1": links[1].name,
                "read2": links.get(2).name if paired else "",
                "bytes_read1": links[1].stat().st_size,
                "bytes_read2": links[2].stat().st_size if paired else 0,
            }
        )

    manifest = outdir / workspace.MANIFEST_NAME
    with open(manifest, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()) if records else ["sample"])
        writer.writeheader()
        writer.writerows(records)

    paired_count = sum(1 for record in records if record["paired"] == "yes")
    print(f"[Info] omics={args.omics} dataset={args.dataset}")
    print(f"[Info] samples={len(records)} paired={paired_count} single={len(records) - paired_count}")
    print(f"[Info] manifest={manifest}")

    if unpaired:
        print(f"[ERROR] incomplete read pairs: {', '.join(unpaired)}", file=sys.stderr)
        return 1
    if missing:
        message = f"FASTQ not found for {len(missing)} sample(s): " + "; ".join(missing)
        if args.allow_missing:
            print(f"[WARN] {message}", file=sys.stderr)
        else:
            print(f"[ERROR] {message}", file=sys.stderr)
            return 1
    if not records:
        print("[ERROR] no FASTQ file could be assigned to a sample", file=sys.stderr)
        return 1

    (outdir / workspace.FINISH_NAME).write_text(
        f"omics\t{args.omics}\ndataset\t{args.dataset}\nsamples\t{len(records)}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
