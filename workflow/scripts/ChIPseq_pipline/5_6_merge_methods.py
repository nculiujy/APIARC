#!/usr/bin/env python3
"""Combine the three differential-binding criteria into one call set (module 5).

Reviewers asked for more than one definition of differential enrichment. Rather
than picking one and hiding the rest, all criteria are computed independently and
merged here, with the contributing methods recorded per region so the manuscript
can report agreement (and disagreement) explicitly.
"""

from __future__ import annotations

import argparse
import os

METHODS = ("count", "presence", "stat")


def read_table(path: str) -> dict[str, dict]:
    if not os.path.exists(path):
        return {}
    table = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            table[row["region_id"]] = row
    return table


def as_bool(value: str | bool | None) -> bool:
    return str(value).lower() in {"true", "1", "yes"}


def as_float(value, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def main() -> None:
    parser = argparse.ArgumentParser(description="merge differential-binding criteria")
    parser.add_argument("--count", required=True)
    parser.add_argument("--presence", required=True)
    parser.add_argument("--stat", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--mark", required=True)
    parser.add_argument("--comparison", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--lfc", type=float, default=1.0)
    parser.add_argument("--methods", default=",".join(METHODS),
                        help="criteria that contribute a call, comma separated")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    active = [m.strip() for m in args.methods.split(",") if m.strip() in METHODS]
    if not active:
        raise SystemExit(f"--methods must contain at least one of {METHODS}")

    tables = {
        "count": read_table(args.count),
        "presence": read_table(args.presence),
        "stat": read_table(args.stat),
    }

    regions = []
    with open(args.consensus) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            regions.append((fields[3], fields[0], int(fields[1]), int(fields[2]), int(fields[4])))

    rows = []
    for region_id, chrom, start, end, support in regions:
        significant_methods = []
        for method in active:
            row = tables[method].get(region_id)
            if row and as_bool(row.get("significant")):
                significant_methods.append(method)

        primary = tables["count"].get(region_id)
        log2fc = as_float(primary["log2FC"]) if primary else 0.0
        fdr = as_float(primary["FDR"], 1.0) if primary else 1.0
        if not significant_methods:
            category = "static"
        else:
            direction_values = [
                as_float(tables[m][region_id]["log2FC"]) for m in significant_methods
                if tables[m].get(region_id)
            ]
            positive = sum(1 for value in direction_values if value > 0)
            negative = sum(1 for value in direction_values if value < 0)
            if positive and negative:
                category = "discordant"
            else:
                category = "diff_up" if positive >= negative else "diff_down"

        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "region_id": region_id,
                "mark": args.mark,
                "comparison": args.comparison,
                "n_support": support,
                "log2FC": log2fc,
                "FDR": fdr,
                "methods": ";".join(significant_methods),
                "n_methods": len(significant_methods),
                "category": category,
            }
        )

    os.makedirs(args.outdir, exist_ok=True)
    out_path = os.path.join(args.outdir, "differential_regions.tsv")
    with open(out_path, "w") as out:
        columns = [
            "chrom", "start", "end", "region_id", "mark", "comparison",
            "n_support", "log2FC", "FDR", "methods", "n_methods", "category",
        ]
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in columns) + "\n")

    significant_sets = {
        method: {rid for rid, row in tables[method].items() if as_bool(row.get("significant"))}
        for method in METHODS
    }
    with open(os.path.join(args.outdir, "method_overlap.tsv"), "w") as out:
        out.write("comparison\tmethod_a\tmethod_b\tn_a\tn_b\tintersection\tjaccard\n")
        for i, method_a in enumerate(METHODS):
            for method_b in METHODS[i + 1:]:
                a, b = significant_sets[method_a], significant_sets[method_b]
                out.write(
                    f"{args.comparison}\t{method_a}\t{method_b}\t{len(a)}\t{len(b)}\t"
                    f"{len(a & b)}\t{jaccard(a, b):.4f}\n"
                )

    counts = {}
    for row in rows:
        counts[row["category"]] = counts.get(row["category"], 0) + 1
    with open(os.path.join(args.outdir, "diff_summary.tsv"), "w") as out:
        out.write("comparison\tcategory\tn_regions\n")
        for category in sorted(counts):
            out.write(f"{args.comparison}\t{category}\t{counts[category]}\n")
        for method in METHODS:
            out.write(f"{args.comparison}\tsig_{method}\t{len(significant_sets[method])}\n")

    print(f"[merge] {args.comparison}: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
