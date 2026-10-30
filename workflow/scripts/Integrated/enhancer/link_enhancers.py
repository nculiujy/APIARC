#!/usr/bin/env python3
"""Distal-region to gene links with a distance-matched random control.

The previous enhancer module measured the distance from each distal peak to a
file of *promoter peaks*, so every reported "enhancer-gene distance" was really
the distance between two binding sites. Distances here are measured to the
transcription start site of the nearest gene in the gene model, cached once per
chromosome.

Two statistical guards are applied, because they are the difference between a
result and an artefact:

* Regions are aggregated per gene before any test. A gene with twelve distal
  regions is one observation, not twelve, and treating peaks as independent
  observations is the standard way to manufacture significance.
* The null is distance matched. The probability that a gene is differentially
  expressed depends strongly on how far the region sits from a TSS, so a
  genome-wide uniform sample would attribute that trend to the enhancer. Random
  positions are therefore drawn genome-wide, binned by their TSS distance, and
  compared with the observed regions *within each bin*.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from common.annotate import assign_gene, load_genes


def read_tsv(path: str) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def bh(pvalues: list[float]) -> list[float]:
    indexed = sorted(
        ((value, position) for position, value in enumerate(pvalues) if value is not None),
        key=lambda item: item[0],
    )
    n = len(indexed)
    adjusted = [None] * len(pvalues)
    previous = 1.0
    for rank in range(n, 0, -1):
        value, position = indexed[rank - 1]
        previous = min(previous, value * n / rank)
        adjusted[position] = previous
    return adjusted


def fisher_two_sided(a: int, b: int, c: int, d: int) -> tuple[float, float, float, float]:
    """Fisher exact test on a 2x2 table, with odds ratio and Wald CI."""
    n = a + b + c + d
    if n == 0:
        return (1.0, float("nan"), float("nan"), float("nan"))
    def log_choose(n_, k_):
        return math.lgamma(n_ + 1) - math.lgamma(k_ + 1) - math.lgamma(n_ - k_ + 1)

    row1, col1 = a + b, a + c
    log_denominator = log_choose(n, col1)
    def probability(x):
        return math.exp(
            log_choose(row1, x) + log_choose(n - row1, col1 - x) - log_denominator
        )

    observed = probability(a)
    low, high = max(0, col1 - (n - row1)), min(row1, col1)
    p = sum(
        value for value in (probability(x) for x in range(low, high + 1))
        if value <= observed * (1 + 1e-7)
    )
    p = min(1.0, p)

    zero = 0.5
    odds = ((a + zero) * (d + zero)) / ((b + zero) * (c + zero))
    if min(a, b, c, d) > 0:
        log_or = math.log(odds)
        standard_error = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        ci_low = math.exp(log_or - 1.96 * standard_error)
        ci_high = math.exp(log_or + 1.96 * standard_error)
    else:
        ci_low = ci_high = float("nan")
    return (p, odds, ci_low, ci_high)


def distance_bin(distance: int, bins: list[int]) -> str:
    for low, high in zip(bins, bins[1:]):
        if low <= distance < high:
            return f"{low}-{high}"
    return f">={bins[-1]}"


def main() -> None:
    parser = argparse.ArgumentParser(description="enhancer-gene links with distance-matched controls")
    parser.add_argument("--regions", required=True, help="regions.tsv from 0_stats")
    parser.add_argument("--genes", required=True, help="GENCODE BED")
    parser.add_argument("--rna", action="append", default=[], metavar="LABEL=PATH",
                        help="DESeq2 table for one comparison, as comparison=path")
    parser.add_argument("--bins", default="0,10000,50000,100000,500000,1000000")
    parser.add_argument("--promoter-window", type=int, default=2000)
    parser.add_argument("--max-distance", type=int, default=1000000)
    parser.add_argument("--random-per-region", type=int, default=5)
    parser.add_argument("--min-base-mean", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    random.seed(args.seed)
    bins = [int(value) for value in args.bins.split(",")]
    index = load_genes(args.genes)

    rna: dict[str, dict] = {}
    rna_by_comparison: dict[str, dict[str, dict]] = {}
    for item in args.rna:
        label, _, path = item.partition("=")
        table: dict[str, dict] = {}
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                symbol = (row.get("SYMBOL") or "").strip()
                if symbol and symbol not in table:
                    table[symbol] = row
        rna_by_comparison[label] = table
        for symbol, row in table.items():
            rna.setdefault(symbol, row)
    if not rna:
        raise SystemExit("[enhancer] no --rna table was given")

    def de_table(comparison: str) -> dict[str, dict]:
        # Each comparison is paired with its own contrast; pooling them would
        # import the 2 h expression calls into the 12 h binding result.
        return rna_by_comparison.get(comparison) or next(iter(rna_by_comparison.values()))

    def is_de(gene: str, comparison: str = "") -> bool:
        row = de_table(comparison).get(gene) if comparison else rna.get(gene)
        if row is None:
            return False
        padj = to_float(row.get("padj"))
        return padj is not None and padj < 0.05

    def is_de_any(gene: str) -> bool:
        """DE in any paired contrast - the bin test pools every comparison."""
        return any(is_de(gene, label) for label in rna_by_comparison)

    def is_expressed(gene: str) -> bool:
        row = rna.get(gene)
        if row is None:
            return False
        base_mean = to_float(row.get("baseMean"), -1.0)
        return base_mean is not None and base_mean >= args.min_base_mean

    regions = read_tsv(args.regions)
    distal = [row for row in regions if row.get("region_type") == "distal" and row.get("gene")]
    if not distal:
        raise SystemExit("[enhancer] no distal regions with an assigned gene")

    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in distal:
        grouped[(row["gene"], row["mark"], row["comparison"])].append(row)

    links = []
    for (gene, mark, comparison), items in grouped.items():
        categories = {item["category"] for item in items}
        category = (
            "diff_up" if "diff_up" in categories else
            "diff_down" if "diff_down" in categories else
            "discordant" if "discordant" in categories else
            "bound" if "bound" in categories else "static"
        )
        links.append(
            {
                "gene": gene,
                "mark": mark,
                "comparison": comparison,
                "n_distal_regions": len(items),
                "min_distance": min(int(item["distance"]) for item in items),
                "median_distance": sorted(int(item["distance"]) for item in items)[
                    len(items) // 2
                ],
                "category": category,
                "expressed": "true" if is_expressed(gene) else "false",
                "differential": "true" if is_de(gene, comparison) else "false",
            }
        )
    links.sort(key=lambda row: (row["mark"], row["comparison"], row["gene"]))

    os.makedirs(args.outdir, exist_ok=True)
    link_columns = [
        "gene", "mark", "comparison", "n_distal_regions", "min_distance",
        "median_distance", "category", "expressed", "differential",
    ]
    with open(os.path.join(args.outdir, "enhancer_gene_links.tsv"), "w") as out:
        out.write("\t".join(link_columns) + "\n")
        for row in links:
            out.write("\t".join(str(row[column]) for column in link_columns) + "\n")

    chrom_span = {}
    for chrom, genes in index.genes.items():
        chrom_span[chrom] = max(gene.end for gene in genes)
    spans = [(chrom, length) for chrom, length in chrom_span.items() if length > 0]
    if not spans:
        raise SystemExit("[enhancer] gene model has no usable chromosome spans")
    total_span = sum(length for _, length in spans)

    def random_position() -> tuple[str, int]:
        position = random.randrange(total_span)
        for chrom, length in spans:
            if position < length:
                return chrom, position
            position -= length
        return spans[-1][0], spans[-1][1] - 1

    null_gene: dict[str, list[str]] = defaultdict(list)
    for _ in range(len(distal) * args.random_per_region):
        chrom, position = random_position()
        symbol, distance, region_type = assign_gene(
            index, chrom, position, position + 1, args.promoter_window, args.max_distance
        )
        if symbol and region_type == "distal":
            null_gene[distance_bin(distance, bins)].append(symbol)

    observed_gene: dict[str, set[str]] = defaultdict(set)
    for row in distal:
        observed_gene[distance_bin(int(row["distance"]), bins)].add(row["gene"])

    bin_rows = []
    pvalues = []
    for bin_label in sorted(observed_gene, key=lambda value: bins[0] if value.startswith(">=") else int(value.split("-")[0])):
        observed = {gene for gene in observed_gene[bin_label] if is_expressed(gene)}
        null = {gene for gene in null_gene.get(bin_label, []) if is_expressed(gene)}
        if len(observed) < 10 or len(null) < 10:
            continue
        a = sum(1 for gene in observed if is_de_any(gene))
        b = len(observed) - a
        c = sum(1 for gene in null if is_de_any(gene))
        d = len(null) - c
        p, odds, ci_low, ci_high = fisher_two_sided(a, b, c, d)
        pvalues.append(p)
        bin_rows.append(
            {
                "distance_bin": bin_label,
                "n_observed": len(observed),
                "n_null": len(null),
                "de_observed": a,
                "de_null": c,
                "odds_ratio": round(odds, 4) if odds == odds else "NA",
                "ci_low": round(ci_low, 4) if ci_low == ci_low else "NA",
                "ci_high": round(ci_high, 4) if ci_high == ci_high else "NA",
                "p": p,
            }
        )

    adjusted = bh(pvalues)
    for row, fdr in zip(bin_rows, adjusted):
        row["FDR"] = fdr
    with open(os.path.join(args.outdir, "distance_bins.tsv"), "w") as out:
        columns = [
            "distance_bin", "n_observed", "n_null", "de_observed", "de_null",
            "odds_ratio", "ci_low", "ci_high", "p", "FDR",
        ]
        out.write("\t".join(columns) + "\n")
        for row in bin_rows:
            out.write("\t".join(str(row[column]) for column in columns) + "\n")

    promoter = [row for row in regions if row.get("region_type") == "promoter" and row.get("gene")]
    promoter_genes = {row["gene"] for row in promoter if is_expressed(row["gene"])}
    distal_genes = {row["gene"] for row in distal if is_expressed(row["gene"])}
    summary_rows = []
    if distal_genes and promoter_genes:
        a = sum(1 for gene in distal_genes if is_de_any(gene))
        b = len(distal_genes) - a
        c = sum(1 for gene in promoter_genes if is_de_any(gene))
        d = len(promoter_genes) - c
        p, odds, ci_low, ci_high = fisher_two_sided(a, b, c, d)
        summary_rows.append(
            {
                "analysis": "distal_vs_promoter_de_fraction",
                "n_distal_genes": len(distal_genes),
                "n_promoter_genes": len(promoter_genes),
                "de_distal": a,
                "de_promoter": c,
                "odds_ratio": round(odds, 4) if odds == odds else "NA",
                "ci_low": round(ci_low, 4) if ci_low == ci_low else "NA",
                "ci_high": round(ci_high, 4) if ci_high == ci_high else "NA",
                "p": p,
            }
        )
    with open(os.path.join(args.outdir, "enhancer_summary.tsv"), "w") as out:
        columns = [
            "analysis", "n_distal_genes", "n_promoter_genes", "de_distal",
            "de_promoter", "odds_ratio", "ci_low", "ci_high", "p",
        ]
        out.write("\t".join(columns) + "\n")
        for row in summary_rows:
            out.write("\t".join(str(row[column]) for column in columns) + "\n")

    print(
        f"[enhancer] {len(links)} gene links, {len(bin_rows)} distance bins, "
        f"{len(distal)} distal regions -> {args.outdir}"
    )


if __name__ == "__main__":
    main()
