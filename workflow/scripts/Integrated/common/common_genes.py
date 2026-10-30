#!/usr/bin/env python3
"""Common genes of the joint analysis: differential binding plus differential expression.

The gene table written by ``0_2_gene_table`` collapses every region of a gene to
the single strongest log2 fold change, which is the right summary for testing but
the wrong one for a scatter plot: a gene with eight peaks is then represented by
its most extreme one. This step keeps three aggregates (strongest, equal-weight
mean, signal-weighted mean) and joins the per-sample values that the figures
need.

Two definitions of "common" are written side by side. ``concordant`` requires
binding and expression to move in the same direction; the full set stays in the
table behind ``concordant=false`` so a plot can be redrawn under either rule
without rerunning the stage.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict


def to_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def is_number(value) -> bool:
    return to_float(value) is not None


def read_tsv(path: str) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_metadata(path: str) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def read_rna(path: str) -> tuple[list[str], dict[str, dict]]:
    """Return the sample columns and a symbol-keyed table.

    SYMBOL repeats for some Ensembl ids; the first row wins, matching the
    expression universe used by the enrichment stage.
    """
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        fields = [name for name in (reader.fieldnames or []) if name not in (
            "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue",
            "padj", "SYMBOL",
        )]
        table: dict[str, dict] = {}
        for row in reader:
            symbol = (row.get("SYMBOL") or "").strip()
            if not symbol or symbol in table:
                continue
            table[symbol] = row
    return fields, table


def rankdata(values: list[float]) -> list[float]:
    """Average ranks, the tie handling Spearman's rho requires."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position
        while end + 1 < len(order) and values[order[end + 1]] == values[order[position]]:
            end += 1
        mean_rank = (position + end) / 2.0 + 1.0
        for index in range(position, end + 1):
            ranks[order[index]] = mean_rank
        position = end + 1
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 3:
        return float("nan")
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sxy = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    sxx = sum((a - mean_x) ** 2 for a in x)
    syy = sum((b - mean_y) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def spearman(pairs: list[tuple[float, float]]) -> tuple[float, int]:
    if len(pairs) < 3:
        return float("nan"), len(pairs)
    x = [pair[0] for pair in pairs]
    y = [pair[1] for pair in pairs]
    return pearson(rankdata(x), rankdata(y)), len(pairs)


def timepoint_of(group: str, mark: str) -> str:
    """H3K27ac_TGFb12h -> TGFb12h. The mark never contains an underscore."""
    prefix = f"{mark}_"
    return group[len(prefix):] if group.startswith(prefix) else group


def main() -> None:
    parser = argparse.ArgumentParser(description="common genes of the joint analysis")
    parser.add_argument("--gene-table", required=True)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--rna", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--chip-metadata", required=True)
    parser.add_argument("--rna-metadata", required=True)
    parser.add_argument("--normalized-counts", action="append", default=[],
                        metavar="COMPARISON=PATH",
                        help="region x sample normalized counts of one comparison")
    parser.add_argument("--enhancer-links", default="")
    parser.add_argument("--rna-alpha", type=float, default=0.05)
    parser.add_argument("--min-base-mean", type=float, default=10.0)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # --- sample design -----------------------------------------------------
    rna_samples = read_metadata(args.rna_metadata)
    chip_samples = read_metadata(args.chip_metadata)
    # The RNA metadata names a time point and a treatment group but carries no
    # per-sample condition column, so the replicate number comes from the order
    # in which the samples of one time point appear. Without it every column of
    # a time point collapses onto the same name.
    rna_label_of: dict[str, str] = {}
    rna_label_order: list[str] = []
    seen_name: dict[str, int] = defaultdict(int)
    for row in rna_samples:
        sample = (row.get("sample") or "").strip()
        name = (row.get("sample_name") or "").strip()
        if not sample or not name:
            continue
        condition = (row.get("condition") or "").strip()
        if not condition:
            seen_name[name] += 1
            condition = f"rep{seen_name[name]}"
        rna_label_of[sample] = f"{name}_{condition}"
        if rna_label_of[sample] not in rna_label_order:
            rna_label_order.append(rna_label_of[sample])
    rna_columns: list[str] = []

    chip_by_mark: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    chip_group_of: dict[str, str] = {}
    for row in chip_samples:
        if (row.get("role") or "").strip() != "IP":
            continue
        sample = (row.get("sample") or "").strip()
        group = (row.get("sample_name") or "").strip()
        mark = group.split("_")[0]
        chip_by_mark[mark].append((sample, group, (row.get("condition") or "").strip()))
        chip_group_of[sample] = group
    for mark in chip_by_mark:
        chip_by_mark[mark].sort(key=lambda item: (timepoint_of(item[1], mark), item[2]))

    # --- expression --------------------------------------------------------
    expression: dict[str, dict[str, str]] = {}
    rna_tables: dict[str, dict[str, dict]] = {}
    for item in args.rna:
        label, _, path = item.partition("=")
        fields, table = read_rna(path)
        rna_tables[label] = table
        # One DESeq2 table holds only the samples of its own contrast, so the
        # column order is the union over the tables, not the metadata order.
        for sample in fields:
            rna_label = rna_label_of.get(sample, sample)
            if rna_label not in rna_columns:
                rna_columns.append(rna_label)
        for symbol, row in table.items():
            target = expression.setdefault(symbol, {})
            for sample in fields:
                target[rna_label_of.get(sample, sample)] = (row.get(sample) or "").strip()
    if not rna_tables:
        raise SystemExit("[common] no --rna table was given")
    # Each DESeq2 table lists its own contrast, so the tables arrive in an
    # arbitrary order (12 h first here). The metadata order is the design order
    # of the time course; the figures draw the x axis in column order.
    rna_rank = {label: index for index, label in enumerate(rna_label_order)}
    rna_columns.sort(key=lambda label: rna_rank.get(label, len(rna_rank)))

    # --- per-sample ChIP signal, restricted to regions of interest ----------
    signal_of_region: dict[tuple, dict[str, float]] = {}
    for item in args.normalized_counts:
        comparison, _, path = item.partition("=")
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (comparison, row["region_id"])
                signal_of_region[key] = {
                    sample: to_float(row.get(sample), 0.0) or 0.0
                    for sample in row if sample != "region_id"
                }

    # --- diff regions, grouped by gene -------------------------------------
    diff_regions: dict[tuple, list[dict]] = defaultdict(list)
    with open(args.regions) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("category") not in ("diff_up", "diff_down"):
                continue
            gene = (row.get("gene") or "").strip()
            if not gene:
                continue
            if not is_number(row.get("log2FC")) or not is_number(row.get("FDR")):
                continue
            key = (gene, row["mark"], row["comparison"], row["region_type"])
            diff_regions[key].append(row)

    enhancer_distance: dict[tuple, float] = {}
    if args.enhancer_links and os.path.exists(args.enhancer_links):
        for row in read_tsv(args.enhancer_links):
            value = to_float(row.get("median_distance"))
            if value is None:
                continue
            enhancer_distance[(row["gene"], row["mark"], row["comparison"])] = value

    # --- gene table: the RNA side decided by 0_stats ---------------------
    gene_table_rows = read_tsv(args.gene_table)
    rna_by_key = {
        (row["gene"], row["mark"], row["comparison"], row["region_type"]): row
        for row in gene_table_rows
    }

    # --- gene-level rows ---------------------------------------------------
    columns = [
        "gene", "mark", "comparison", "region_type", "chip_direction",
        "rna_direction", "concordant", "n_regions", "chip_lfc_strongest",
        "chip_lfc_mean", "chip_lfc_weighted", "chip_fdr", "chip_methods",
        "rna_baseMean", "rna_log2FC", "rna_padj", "enhancer_distance",
    ]
    rows = []
    signal_rows = []
    for key, items in diff_regions.items():
        gene, mark, comparison, region_type = key
        if comparison not in rna_tables:
            continue
        # The RNA values come from the gene table, not from a lookup by symbol in
        # the DESeq2 table: two Ensembl ids can share one symbol, and 0_stats
        # already decided which of them the gene row describes. Relooking the
        # symbol up would pick the other id for those genes and the row count
        # would stop matching chip_gene_table.tsv.
        gene_row = rna_by_key.get(key, {})
        base_mean = to_float(gene_row.get("rna_baseMean"), -1.0) or -1.0
        padj = to_float(gene_row.get("rna_padj"))
        rna_lfc = to_float(gene_row.get("rna_log2FC"))
        if padj is None or rna_lfc is None:
            continue
        if base_mean < args.min_base_mean or padj >= args.rna_alpha:
            continue

        chip_values = [to_float(item["log2FC"]) for item in items]
        # Up wins when a gene carries both an up and a down region, matching the
        # priority order of the gene table.
        chip_direction = "up" if any(item["category"] == "diff_up" for item in items) else "down"
        rna_direction = "up" if rna_lfc > 0 else "down"

        weights = []
        for item in items:
            sample_values = signal_of_region.get((comparison, item["region_id"]), {})
            observed = [value for value in sample_values.values() if value > 0]
            weights.append(sum(observed) / len(observed) if observed else 0.0)
        if sum(weights) > 0:
            weighted = sum(lfc * weight for lfc, weight in zip(chip_values, weights)) / sum(weights)
        else:
            weighted = sum(chip_values) / len(chip_values)

        rows.append(
            {
                "gene": gene,
                "mark": mark,
                "comparison": comparison,
                "region_type": region_type,
                "chip_direction": chip_direction,
                "rna_direction": rna_direction,
                "concordant": "true" if chip_direction == rna_direction else "false",
                "n_regions": len(items),
                "chip_lfc_strongest": max(chip_values, key=abs),
                "chip_lfc_mean": sum(chip_values) / len(chip_values),
                "chip_lfc_weighted": weighted,
                "chip_fdr": min(to_float(item["FDR"]) for item in items),
                "chip_methods": max(int(to_float(item.get("n_methods"), 0) or 0) for item in items),
                "rna_baseMean": base_mean,
                "rna_log2FC": rna_lfc,
                "rna_padj": padj,
                "enhancer_distance": "" if region_type == "promoter"
                else enhancer_distance.get((gene, mark, comparison), ""),
            }
        )

        for item in items:
            sample_values = signal_of_region.get((comparison, item["region_id"]), {})
            for sample, value in sample_values.items():
                group = chip_group_of.get(sample, "")
                signal_rows.append(
                    {
                        "gene": gene,
                        "mark": mark,
                        "comparison": comparison,
                        "region_type": region_type,
                        "region_id": item["region_id"],
                        "chrom": item["chrom"],
                        "start": item["start"],
                        "end": item["end"],
                        "sample": sample,
                        "group": group,
                        "timepoint": timepoint_of(group, mark) if group else "",
                        "value": value,
                    }
                )

    rows.sort(key=lambda row: (row["mark"], row["comparison"], row["region_type"],
                               -abs(row["chip_lfc_mean"]), row["gene"]))

    def write(path: str, header: list[str], payload: list[dict]) -> None:
        with open(path, "w") as out:
            out.write("\t".join(header) + "\n")
            for row in payload:
                out.write("\t".join(str(row.get(name, "")) for name in header) + "\n")

    write(os.path.join(args.outdir, "common_genes.tsv"), columns, rows)
    write(
        os.path.join(args.outdir, "signal_by_region.tsv"),
        ["gene", "mark", "comparison", "region_type", "region_id", "chrom", "start",
         "end", "sample", "group", "timepoint", "value"],
        signal_rows,
    )

    # --- expression matrix -------------------------------------------------
    genes = sorted({row["gene"] for row in rows})
    write(
        os.path.join(args.outdir, "expression_matrix.tsv"),
        ["gene", *rna_columns],
        [{"gene": gene, **expression.get(gene, {})} for gene in genes],
    )

    # --- signal matrix, one row per common-gene row ------------------------
    signal_matrix = []
    for row in rows:
        key = (row["gene"], row["mark"], row["comparison"], row["region_type"])
        items = diff_regions[key]
        per_sample: dict[str, list[float]] = defaultdict(list)
        for item in items:
            for sample, value in signal_of_region.get((row["comparison"], item["region_id"]), {}).items():
                per_sample[sample].append(value)
        record = {
            "gene": row["gene"], "mark": row["mark"], "comparison": row["comparison"],
            "region_type": row["region_type"], "n_regions": row["n_regions"],
            # The figure splits by direction; carrying it here avoids a second
            # join in R, which has no merge helper in this environment.
            "chip_direction": row["chip_direction"], "concordant": row["concordant"],
        }
        for sample, group, condition in chip_by_mark.get(row["mark"], []):
            timepoint = timepoint_of(group, row["mark"])
            values = per_sample.get(sample, [])
            record[f"{timepoint}_{condition}"] = sum(values) / len(values) if values else ""
        signal_matrix.append(record)

    signal_columns = []
    if chip_by_mark:
        first_mark = rows[0]["mark"] if rows else ""
        for _sample, group, condition in chip_by_mark.get(first_mark, []):
            signal_columns.append(f"{timepoint_of(group, first_mark)}_{condition}")
    write(
        os.path.join(args.outdir, "signal_matrix.tsv"),
        ["gene", "mark", "comparison", "region_type", "n_regions",
         "chip_direction", "concordant", *signal_columns],
        signal_matrix,
    )

    # --- sample map --------------------------------------------------------
    sample_map = []
    for row in rna_samples:
        sample = (row.get("sample") or "").strip()
        name, _, condition = rna_label_of.get(sample, "").rpartition("_")
        fallback = (row.get("sample_name") or "").strip()
        sample_map.append({
            "omics": "RNAseq", "sample": sample,
            "group": name or fallback, "condition": condition,
            "timepoint": name or fallback,
        })
    sample_map += [
        {"omics": "ChIPseq", "sample": sample, "group": group, "condition": condition,
         "timepoint": timepoint_of(group, group.split("_")[0])}
        for mark in chip_by_mark for sample, group, condition in chip_by_mark[mark]
    ]
    write(os.path.join(args.outdir, "sample_map.tsv"),
          ["omics", "sample", "group", "condition", "timepoint"], sample_map)

    # --- correlation reference --------------------------------------------
    # The common set is filtered on both axes, so a correlation computed inside
    # it is not comparable to the genome-wide one. Both are written here; the
    # figure prints the genome-wide value in its subtitle.
    stats_rows = []
    grouped_all: dict[tuple, list[dict]] = defaultdict(list)
    for row in gene_table_rows:
        if row.get("expressed") != "true":
            continue
        grouped_all[(row["mark"], row["comparison"], row["region_type"])].append(row)

    for key, items in sorted(grouped_all.items()):
        mark, comparison, region_type = key
        pairs = [
            (to_float(item["log2FC"]), to_float(item["rna_log2FC"]))
            for item in items
            if is_number(item.get("log2FC")) and is_number(item.get("rna_log2FC"))
        ]
        rho_all, n_all = spearman(pairs)
        subset = [row for row in rows
                  if (row["mark"], row["comparison"], row["region_type"]) == key]
        for label, selected in (
            ("common_concordant", [row for row in subset if row["concordant"] == "true"]),
            ("common_any", subset),
        ):
            rho, n = spearman([(row["chip_lfc_mean"], row["rna_log2FC"]) for row in selected])
            stats_rows.append(
                {
                    "mark": mark, "comparison": comparison, "region_type": region_type,
                    "set": label, "n": n, "spearman": rho,
                    "spearman_all_genes": rho_all, "n_all_genes": n_all,
                    "note": "common-set rho is computed after two filters; "
                            "cite spearman_all_genes as the unbiased reference",
                }
            )
    write(
        os.path.join(args.outdir, "common_stats.tsv"),
        ["mark", "comparison", "region_type", "set", "n", "spearman",
         "spearman_all_genes", "n_all_genes", "note"],
        stats_rows,
    )

    n_concordant = sum(1 for row in rows if row["concordant"] == "true")
    print(
        f"[common] {len(rows)} common gene rows ({n_concordant} concordant), "
        f"{len(genes)} genes, {len(signal_rows)} region/sample values -> {args.outdir}"
    )


if __name__ == "__main__":
    main()
