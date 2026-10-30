#!/usr/bin/env python3
"""Assign every candidate region to a gene (integrated analysis, step 1).

In differential mode the input is ``differential_regions.tsv`` from module 5.
In descriptive mode (a single ChIP condition) the input is the merged peak set of
each mark/condition group. The output is a single region table that keeps the
region type (promoter vs distal) explicit, so downstream tests and figures never
mix the two regulatory layers.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from common.annotate import assign_gene, load_genes

COLUMNS = [
    "region_id", "chrom", "start", "end", "mark", "comparison", "region_type",
    "gene", "distance", "log2FC", "FDR", "category", "n_methods",
]


def read_differential(path: str) -> list[dict]:
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values))
            rows.append(
                {
                    "region_id": row["region_id"],
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "mark": row["mark"],
                    "comparison": row["comparison"],
                    "log2FC": row.get("log2FC", ""),
                    "FDR": row.get("FDR", ""),
                    "category": row.get("category", "static"),
                    "n_methods": row.get("n_methods", "0"),
                }
            )
    return rows


def read_bound(group_bed: str, group: str) -> list[dict]:
    mark, _, condition = group.rpartition("_")
    rows = []
    with open(group_bed) as handle:
        for index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "region_id": f"{fields[0]}:{fields[1]}-{fields[2]}",
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "mark": mark or group,
                    "comparison": f"bound_{condition}",
                    "log2FC": "",
                    "FDR": "",
                    "category": "bound",
                    "n_methods": "0",
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="annotate regions with genes")
    parser.add_argument("--mode", required=True, choices=["differential", "descriptive"])
    parser.add_argument("--genes", required=True, help="GENCODE BED")
    parser.add_argument("--promoter-window", type=int, default=2000)
    parser.add_argument("--max-distance", type=int, default=1000000)
    # Written as `--differential a.tsv --differential b.tsv` by the rule, one
    # flag per comparison. With nargs="*" and the default store action every
    # repeat replaced the previous list, so only the last mark was annotated
    # and the other two marks vanished from region_sets.tsv.
    parser.add_argument("--differential", action="append", nargs="+", default=[],
                        help="one or more differential_regions.tsv tables")
    parser.add_argument("--bound-bed", action="append", default=[],
                        help="group=TABseparated path to merged peaks")
    parser.add_argument("--out", required=True)
    parser.add_argument("--beddir", required=True)
    args = parser.parse_args()

    genes = load_genes(args.genes)

    records: list[dict] = []
    if args.mode == "differential":
        for group in args.differential:
            for path in group:
                records.extend(read_differential(path))
    else:
        for item in args.bound_bed:
            group, _, path = item.partition("=")
            records.extend(read_bound(path, group))

    if not records:
        raise SystemExit("[annotate] no regions to annotate")

    os.makedirs(args.beddir, exist_ok=True)
    exported: dict[tuple, list[str]] = {}
    with open(args.out, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in records:
            symbol, distance, region_type = assign_gene(
                genes, row["chrom"], row["start"], row["end"],
                args.promoter_window, args.max_distance,
            )
            row["region_type"] = region_type
            row["gene"] = symbol
            row["distance"] = distance
            out.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")

            if symbol and region_type != "unassigned":
                key = (row["mark"], row["comparison"], row["category"], region_type)
                exported.setdefault(key, []).append(
                    f"{row['chrom']}\t{row['start']}\t{row['end']}\t{symbol}\n"
                )

    manifest_rows = [
        "set_id\tbed\tn_regions\tn_genes\tmark\tcomparison\tcategory\tregion_type"
    ]
    for (mark, comparison, category, region_type), lines in exported.items():
        filename = f"{mark}__{comparison}__{category}__{region_type}.bed"
        path = os.path.join(args.beddir, filename)
        with open(path, "w") as handle:
            handle.writelines(lines)
        genes_in_set = {line.rstrip("\n").split("\t")[3] for line in lines}
        manifest_rows.append(
            "\t".join([
                os.path.splitext(filename)[0], path, str(len(lines)),
                str(len(genes_in_set)), mark, comparison, category, region_type,
            ])
        )

    with open(os.path.join(os.path.dirname(args.out), "region_sets.tsv"), "w") as handle:
        handle.write("\n".join(manifest_rows) + "\n")

    print(
        f"[annotate] {len(records)} regions -> {args.out}; "
        f"{len(exported)} BED sets"
    )


if __name__ == "__main__":
    main()
