#!/usr/bin/env python3
"""Merge replicate peaks into one BED per mark/condition group.

Descriptive mode needs a bound-region set per group, and the integrated layer
needs to know how many replicates actually support each region. The previous
version merged peaks without keeping that count, so a region seen in one of two
replicates was indistinguishable from a region seen in all of them.

The output is BED5: ``chrom start end region_id n_supporting_replicates``.
Merging is done with ``bedtools merge`` over the concatenated replicate BEDs
after tagging each interval with its sample id, which keeps the support count
exact rather than inferred.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common.metadata import read_chip_metadata


def run(cmd, stdout=None):
    print("[merge] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=stdout)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--peakdir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-support", type=int, default=1,
                        help="minimum number of replicates supporting a merged region")
    args = parser.parse_args()

    if shutil.which("bedtools") is None:
        raise SystemExit("bedtools not found on PATH")

    samples = read_chip_metadata(args.metadata)
    groups: dict[str, list] = defaultdict(list)
    for sample in samples:
        groups[sample.group].append(sample)

    summary = ["group\tn_replicates\tmin_support\tn_regions\ttotal_bp"]
    for group in sorted(groups):
        members = sorted(groups[group], key=lambda s: s.replicate)
        group_dir = os.path.join(args.outdir, group)
        os.makedirs(group_dir, exist_ok=True)
        tagged = os.path.join(group_dir, f"{group}_replicates_tagged.bed")
        missing = []
        with open(tagged, "w") as out:
            for sample in members:
                bed = os.path.join(args.peakdir, sample.sample, f"{sample.sample}_peaks.bed")
                if not os.path.exists(bed):
                    missing.append(bed)
                    continue
                with open(bed) as handle:
                    for line in handle:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 3:
                            continue
                        out.write("\t".join(fields[:3] + [sample.sample]) + "\n")
        if missing:
            raise SystemExit("[merge] missing replicate peak beds: " + ", ".join(missing))

        sorted_bed = tagged + ".sorted"
        with open(sorted_bed, "w") as out:
            run(["sort", "-k1,1", "-k2,2n", tagged], stdout=out)

        merged = os.path.join(group_dir, f"{group}_peaks_merged.bed")
        merged_raw = merged + ".raw"
        with open(merged_raw, "w") as out:
            run([
                "bedtools", "merge", "-i", sorted_bed,
                "-c", "4", "-o", "distinct,count",
            ], stdout=out)

        n_regions = 0
        total_bp = 0
        with open(merged_raw) as handle, open(merged, "w") as out:
            for index, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise SystemExit(f"[merge] unexpected bedtools output: {line!r}")
                support = int(fields[4])
                if support < args.min_support:
                    continue
                start, end = int(fields[1]), int(fields[2])
                n_regions += 1
                total_bp += end - start
                out.write(
                    f"{fields[0]}\t{start}\t{end}\t{group}_region_{index}\t{support}\n"
                )
        os.remove(sorted_bed)
        os.remove(merged_raw)
        os.remove(tagged)
        summary.append(
            f"{group}\t{len(members)}\t{args.min_support}\t{n_regions}\t{total_bp}"
        )

    with open(os.path.join(args.outdir, "merged_peaks_summary.tsv"), "w") as out:
        out.write("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
