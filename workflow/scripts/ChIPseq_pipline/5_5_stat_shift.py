#!/usr/bin/env python3
"""Peak-strength shift between conditions (module 5, criterion 3 of 3).

For every consensus region we take the maximum -log10(q) reported by MACS2 in
each replicate (0 when the region is not called) and test the shift between
conditions with Welch's t-test followed by BH correction. This is deliberately a
sensitivity analysis: it answers "does the confidence of peak calling change
between conditions", which is a different question from "does read depth
change".
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common.metadata import comparison_samples, read_chip_metadata


def peak_scores(consensus: str, peak_file: str) -> dict[str, float]:
    """Maximum -log10(q) per consensus region for one library."""
    result = subprocess.run(
        ["bedtools", "intersect", "-a", consensus, "-b", peak_file, "-wa", "-wb", "-nonamecheck"],
        text=True, capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"bedtools intersect failed: {result.stderr}")

    scores: dict[str, float] = {}
    qvalue_column = 10
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        region_id = fields[3]
        qvalue = float(fields[qvalue_column]) if len(fields) > qvalue_column and fields[qvalue_column] else 0.0
        if qvalue > scores.get(region_id, 0.0):
            scores[region_id] = qvalue
    return scores


def bh(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [1.0] * n
    running = 1.0
    for rank in range(n - 1, -1, -1):
        index = order[rank]
        running = min(running, pvalues[index] * n / (rank + 1))
        adjusted[index] = running
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser(description="peak-strength shift test")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--peakdir", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--mark", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--treatment", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    import numpy as np
    from scipy.stats import ttest_ind

    samples = comparison_samples(
        read_chip_metadata(args.metadata), args.reference, args.treatment, args.mark
    )
    ref = [s for s in samples if s.comparison_key == args.reference]
    treat = [s for s in samples if s.comparison_key == args.treatment]
    if len(ref) < 2 or len(treat) < 2:
        raise ValueError("at least two replicates per condition are required")

    region_ids = []
    with open(args.consensus) as handle:
        for line in handle:
            region_ids.append(line.rstrip("\n").split("\t")[3])

    matrices = {}
    for sample in samples:
        peak_file = os.path.join(args.peakdir, sample.sample, f"{sample.sample}_peaks.bed")
        if not os.path.exists(peak_file):
            raise FileNotFoundError(peak_file)
        matrices[sample.sample] = peak_scores(args.consensus, peak_file)

    rows = []
    for region_id in region_ids:
        ref_values = np.array([matrices[s.sample].get(region_id, 0.0) for s in ref])
        treat_values = np.array([matrices[s.sample].get(region_id, 0.0) for s in treat])
        if ref_values.sum() == 0 and treat_values.sum() == 0:
            continue
        log2fc = float(np.log2((treat_values.mean() + 0.1) / (ref_values.mean() + 0.1)))
        pvalue = float(ttest_ind(treat_values, ref_values, equal_var=False).pvalue)
        rows.append({"region_id": region_id, "log2FC": log2fc, "pvalue": pvalue})

    fdr = bh([row["pvalue"] for row in rows])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as out:
        out.write("region_id\tlog2FC\tpvalue\tFDR\tsignificant\tmethod\n")
        for row, qvalue in zip(rows, fdr):
            significant = qvalue < args.fdr
            out.write(
                f"{row['region_id']}\t{row['log2FC']:.4f}\t{row['pvalue']:.6g}\t"
                f"{qvalue:.6g}\t{str(significant).lower()}\tstat\n"
            )
    print(
        f"[stat-shift] {args.mark} {args.treatment}_vs_{args.reference}: "
        f"{sum(1 for q in fdr if q < args.fdr)} significant regions"
    )


if __name__ == "__main__":
    main()
