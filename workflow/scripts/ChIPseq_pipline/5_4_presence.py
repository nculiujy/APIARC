#!/usr/bin/env python3
"""Presence/absence differential binding (module 5, criterion 2 of 3).

A region that is reproducibly called in one condition and never in the other is
evidence of differential binding regardless of fold change. Testing replicate
support with Fisher's exact test keeps this criterion honest instead of reading
a hard threshold off an enrichment ratio.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common.metadata import comparison_samples, read_chip_metadata

PSEUDOCOUNT = 0.5


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
    parser = argparse.ArgumentParser(description="presence/absence differential binding")
    parser.add_argument("--support", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--mark", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--treatment", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    from scipy.stats import fisher_exact

    samples = comparison_samples(
        read_chip_metadata(args.metadata), args.reference, args.treatment, args.mark
    )
    ref_samples = [s.sample for s in samples if s.comparison_key == args.reference]
    treat_samples = [s.sample for s in samples if s.comparison_key == args.treatment]
    if not ref_samples or not treat_samples:
        raise ValueError("missing samples for one of the conditions")

    rows = []
    with open(args.support) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            region_id = fields[index["region_id"]]
            a = sum(int(fields[index[s]]) for s in ref_samples)
            b = sum(int(fields[index[s]]) for s in treat_samples)
            table = [[b, len(treat_samples) - b], [a, len(ref_samples) - a]]
            pvalue = fisher_exact(table, alternative="two-sided")[1]
            rate_treat = (b + PSEUDOCOUNT) / (len(treat_samples) + 1)
            rate_ref = (a + PSEUDOCOUNT) / (len(ref_samples) + 1)
            rows.append(
                {
                    "region_id": region_id,
                    "support_ref": a,
                    "support_treat": b,
                    "log2FC": math.log2(rate_treat / rate_ref),
                    "pvalue": pvalue,
                }
            )

    fdr = bh([row["pvalue"] for row in rows])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as out:
        out.write("region_id\tsupport_ref\tsupport_treat\tlog2FC\tpvalue\tFDR\tsignificant\tmethod\n")
        for row, qvalue in zip(rows, fdr):
            significant = qvalue < args.fdr
            out.write(
                f"{row['region_id']}\t{row['support_ref']}\t{row['support_treat']}\t"
                f"{row['log2FC']:.4f}\t{row['pvalue']:.6g}\t{qvalue:.6g}\t"
                f"{str(significant).lower()}\tpresence\n"
            )
    print(
        f"[presence] {args.mark} {args.treatment}_vs_{args.reference}: "
        f"{sum(1 for q in fdr if q < args.fdr)} significant regions"
    )


if __name__ == "__main__":
    main()
