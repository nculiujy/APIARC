#!/usr/bin/env python3
"""Read counts over consensus regions for every IP and Input library (module 5).

Inputs are counted as well, not because they enter the differential model, but
because they are needed for the negative control (Input vs Input must not
produce differential regions).
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common.metadata import comparison_samples, input_ids, read_chip_metadata


def read_manifest(path: str) -> dict[tuple[str, str], str]:
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in handle if line.strip()]
    return {(row["sample"], row["role"]): row["bam"] for row in rows}


def multicov(bams: list[str], bed: str) -> str:
    cmd = ["bedtools", "multicov", "-bams", *bams, "-bed", bed]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"bedtools multicov failed: {result.stderr}")
    return result.stdout


def count_regions(bams: list[str], consensus: str, threads: int) -> str:
    """Count every region, slicing the file across the threads the rule grants.

    multicov seeks the index once per region, so a consensus of ~150k regions
    x 6 libraries takes the better part of half an hour on one core. The slices
    are independent, and their outputs are stitched back in file order.
    """
    with open(consensus) as handle:
        regions = [line for line in handle if line.strip()]
    slices = max(1, min(int(threads), len(regions) // 1000))
    if slices < 2:
        return multicov(bams, consensus)

    size = len(regions) // slices + 1
    with tempfile.TemporaryDirectory() as tmp:
        paths = []
        for index in range(slices):
            chunk = regions[index * size:(index + 1) * size]
            if not chunk:
                continue
            path = os.path.join(tmp, f"slice{index}.bed")
            with open(path, "w") as handle:
                handle.writelines(chunk)
            paths.append(path)
        with ThreadPoolExecutor(max_workers=len(paths)) as pool:
            return "".join(pool.map(lambda path: multicov(bams, path), paths))


def main() -> None:
    parser = argparse.ArgumentParser(description="counts over consensus regions")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--mark", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--treatment", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    samples = comparison_samples(
        read_chip_metadata(args.metadata), args.reference, args.treatment, args.mark
    )
    if not samples:
        raise ValueError(
            f"no samples for {args.mark}: {args.treatment}_vs_{args.reference}"
        )
    manifest = read_manifest(args.manifest)

    ordered = [s.sample for s in samples]
    inputs = input_ids(samples)
    bams = [manifest[(s, "IP")] for s in ordered] + [manifest[(i, "Input")] for i in inputs]
    labels = ordered + [f"Input_{i}" for i in inputs]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    stdout = count_regions(bams, args.consensus, args.threads)

    # consensus_peaks.bed carries chrom/start/end/region_id/n_support, so the
    # counts start after that many columns, not after the usual three.
    with open(args.consensus) as handle:
        bed_width = len(handle.readline().split("\t"))

    with open(args.out, "w") as out:
        out.write("region_id\t" + "\t".join(labels) + "\n")
        for line in stdout.splitlines():
            fields = line.split("\t")
            region_id = f"{fields[0]}:{fields[1]}-{fields[2]}"
            counts = fields[bed_width:]
            if len(counts) != len(labels):
                raise ValueError(
                    f"{region_id}: {len(counts)} count columns for {len(labels)} libraries"
                )
            out.write(region_id + "\t" + "\t".join(counts) + "\n")

    print(f"[counts] {args.mark} {args.treatment}_vs_{args.reference}: "
          f"{len(stdout.splitlines())} regions x {len(labels)} libraries")


if __name__ == "__main__":
    main()
