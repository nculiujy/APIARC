#!/usr/bin/env python3
"""Aggregate region-level calls to gene level (integrated analysis, step 2).

A gene usually carries several regions. Aggregating with an explicit rule
(strongest evidence wins) avoids the double counting that inflates gene-level
statistics when every peak is treated as an independent observation.

Each comparison is joined to the RNA contrast of the same time point
(``--rna comparison=table``); a single table would compare 12 h binding with
whatever expression contrast happened to be written first.
"""

from __future__ import annotations

import argparse
import csv
import os

CATEGORY_PRIORITY = {"diff_up": 3, "diff_down": 3, "discordant": 2, "static": 1, "bound": 1}


def to_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_regions(path: str) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_rna(path: str) -> dict[str, dict]:
    table = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            symbol = (row.get("SYMBOL") or "").strip()
            if not symbol:
                continue
            table[symbol] = row
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description="gene-level binding table")
    parser.add_argument("--regions", required=True)
    parser.add_argument("--rna", action="append", default=[], metavar="LABEL=PATH",
                        help="DESeq2 table for one comparison, as comparison=path")
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-base-mean", type=float, default=10.0)
    args = parser.parse_args()

    regions = read_regions(args.regions)

    tables: dict[str, dict[str, dict]] = {}
    for item in args.rna:
        label, _, path = item.partition("=")
        tables[label] = read_rna(path)
    if not tables:
        raise SystemExit("[gene-table] no --rna table was given")
    default_table = next(iter(tables.values()))

    grouped: dict[tuple, list[dict]] = {}
    for row in regions:
        gene = (row.get("gene") or "").strip()
        if not gene or row.get("region_type") == "unassigned":
            continue
        key = (gene, row["mark"], row["comparison"], row["region_type"])
        grouped.setdefault(key, []).append(row)

    rows = []
    for (gene, mark, comparison, region_type), items in grouped.items():
        categories = {item["category"] for item in items}
        if "diff_up" in categories:
            category = "diff_up"
        elif "diff_down" in categories:
            category = "diff_down"
        elif "discordant" in categories:
            category = "discordant"
        else:
            category = "bound" if "bound" in categories else "static"

        log2fcs = [to_float(item["log2FC"]) for item in items if item["log2FC"] not in ("", None)]
        fdrs = [to_float(item["FDR"], 1.0) for item in items if item["FDR"] not in ("", None)]
        strongest = max(log2fcs, key=abs) if log2fcs else ""

        rna_row = tables.get(comparison, default_table).get(gene, {})
        base_mean = to_float(rna_row.get("baseMean", ""), -1.0)

        rows.append(
            {
                "gene": gene,
                "mark": mark,
                "comparison": comparison,
                "region_type": region_type,
                "n_regions": len(items),
                "log2FC": strongest,
                "min_FDR": min(fdrs) if fdrs else "",
                "category": category,
                "rna_baseMean": base_mean,
                "rna_log2FC": rna_row.get("log2FoldChange", ""),
                "rna_padj": rna_row.get("padj", ""),
                "expressed": "true" if base_mean >= args.min_base_mean else "false",
            }
        )

    columns = [
        "gene", "mark", "comparison", "region_type", "n_regions", "log2FC",
        "min_FDR", "category", "rna_baseMean", "rna_log2FC", "rna_padj", "expressed",
    ]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")

    print(f"[gene-table] {len(rows)} gene/mark/comparison rows -> {args.out}")


if __name__ == "__main__":
    main()
