#!/usr/bin/env python3
"""TF-gene edges from per-region motif hits.

The previous edges were an annotation, not a regulation claim: for every motif
that was enriched in a peak set, the TF was connected to *every* gene of that
set, which produced 6.59M rows. A set-level AME result cannot say which peak
carries the site, so this step scans the individual sequences.

Design:

* only the common genes of the stage (DEG x differential binding) become edge
  targets, because that is the gene list the report is about;
* only motifs that survived both backgrounds (``motif_background_compare.tsv``)
  and are not present in almost every sequence (``tp_pct``) are scanned - a
  matrix that hits 95% of all peaks cannot discriminate between them;
* each sequence is a ``--window`` bp window centred on the region, which is what
  the plan promised (fimo on 200 bp, not on the whole 2 kb promoter);
* an edge means: this gene has at least one region with a fimo hit for this
  motif. ``n_hit_regions`` and ``best_p`` are what the figure weights.
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import subprocess
from collections import defaultdict

EDGE_COLUMNS = [
    "tf", "gene", "mark", "comparison", "region_type", "category",
    "motif_id", "n_hit_regions", "n_regions_gene", "best_p", "fdr", "evidence_level",
]


def read_tsv(path: str) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value, default: float = 1.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def load_id_to_name(motif_dir: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    path = os.path.join(os.path.dirname(motif_dir.rstrip("/")), "JASPAR_CORE_ID_to_NAME.txt")
    if not os.path.exists(path):
        return mapping
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                mapping[fields[0]] = fields[1].split("::")[0].strip()
    return mapping


def read_meme_block(motif_dir: str, motif_id: str) -> str:
    """``MOTIF`` block of one JASPAR matrix, header stripped."""
    path = os.path.join(motif_dir, f"{motif_id}.meme")
    if not os.path.exists(path):
        return ""
    with open(path) as handle:
        lines = handle.readlines()
    start = next((i for i, line in enumerate(lines) if line.startswith("MOTIF")), None)
    if start is None:
        return ""
    return "".join(lines[start:]).rstrip("\n") + "\n"


def write_meme(blocks: list[str], path: str) -> None:
    with open(path, "w") as out:
        out.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n")
        out.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        out.writelines(blocks)


def robust_motifs(compare_path: str, alpha: float, max_tp_pct: float,
                  max_motifs: int) -> dict[str, list[dict]]:
    """Per set: motifs that beat both backgrounds and still discriminate."""
    selected: dict[str, list[dict]] = defaultdict(list)
    for row in read_tsv(compare_path):
        if row["robust"] != "true":
            continue
        tp_values = [as_float(value, 1.0) for value in
                     (row["tp_pct_static"], row["tp_pct_shuffle"]) if value not in ("", "NA")]
        if tp_values and min(tp_values) >= max_tp_pct:
            continue
        margins = [as_float(value, 1.0) for value in
                   (row["q_static"], row["q_shuffle"]) if value not in ("", "NA")]
        if not margins or max(margins) >= alpha:
            continue
        selected[row["set_id"]].append(
            {
                "motif_id": row["motif_id"],
                "motif_alt_id": row["motif_alt_id"],
                "q": max(margins),
                "category": row["category"],
            }
        )
    for set_id, rows in selected.items():
        rows.sort(key=lambda row: row["q"])
        selected[set_id] = rows[:max_motifs]
    return selected


def centered_window(chrom: str, start: int, end: int, window: int) -> tuple[int, int]:
    middle = (start + end) // 2
    left = max(0, middle - window // 2)
    return left, left + window


def subset_regions(rows: list[dict], genes: set[str], limit: int, seed: int) -> list[dict]:
    rows = [row for row in rows if row["gene"] in genes]
    if limit > 0 and len(rows) > limit:
        keep = sorted(random.Random(seed).sample(range(len(rows)), limit))
        rows = [rows[index] for index in keep]
    return rows


def read_bed(path: str) -> list[dict]:
    rows = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or not fields[3].strip():
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            rows.append({"chrom": fields[0], "start": start, "end": end,
                         "gene": fields[3].strip()})
    return rows


def write_window_bed(regions: list[dict], window: int, path: str) -> dict[str, dict]:
    """BED of centred windows, keyed by the sequence name fimo will report.

    Sequence names are plain ids on purpose. MEME tools read ``chr1:100-200``
    as a coordinate range and then report every hit against that parsed
    interval rather than the window actually scanned, so a name carrying
    coordinates silently moved each hit off its window.

    ``window_start`` is kept with the region: fimo reports offsets inside the
    window it was given, and those have to be shifted back onto the genome.
    """
    index: dict[str, dict] = {}
    seen: set[str] = set()
    with open(path, "w") as out:
        for region in regions:
            left, right = centered_window(region["chrom"], region["start"], region["end"], window)
            label = f"{region['chrom']}:{region['start']}-{region['end']}"
            if label in seen:
                continue
            seen.add(label)
            name = f"w{len(index):06d}"
            index[name] = dict(region, window_start=left, label=label)
            out.write(f"{region['chrom']}\t{left}\t{right}\t{name}\n")
    return index


def region_of(index: dict[str, dict], sequence: str):
    """Region behind a fimo sequence name.

    ``bedtools getfasta -name`` prefixes the BED name and appends the real
    interval when the name is not itself a coordinate range, so fimo reports
    ``chr1:100-200::chr1:120-320`` while the index is keyed on the bare name.
    Without the fallback every hit is dropped and the stage ends with an empty
    network.
    """
    return index.get(sequence) or index.get(sequence.split("::", 1)[0])


def run(cmd: list[str], log_path: str) -> None:
    print("[tf-gene] " + " ".join(cmd), flush=True)
    with open(log_path, "a") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)


def main() -> None:
    parser = argparse.ArgumentParser(description="TF-gene edges from fimo hits")
    parser.add_argument("--motif-hits", required=True)
    parser.add_argument("--background-compare", required=True)
    parser.add_argument("--common", required=True)
    parser.add_argument("--region-sets", required=True)
    parser.add_argument("--motif-dir", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--max-tp-pct", type=float, default=0.9)
    parser.add_argument("--max-motifs", type=int, default=20)
    parser.add_argument("--max-regions", type=int, default=4000)
    parser.add_argument("--window", type=int, default=200)
    parser.add_argument("--thresh", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    for tool in ("fimo", "bedtools"):
        if shutil.which(tool) is None:
            raise SystemExit(f"{tool} not found on PATH")

    os.makedirs(args.outdir, exist_ok=True)
    log_path = os.path.join(args.outdir, "tf_gene_network.log")

    # --- common genes per (key, category) ---------------------------------
    common: dict[tuple, set[str]] = defaultdict(set)
    for row in read_tsv(args.common):
        if row["concordant"] != "true":
            continue
        category = f"diff_{row['chip_direction']}"
        common[(row["mark"], row["comparison"], row["region_type"], category)].add(row["gene"])
    if not common:
        raise SystemExit("[tf-gene] common_genes.tsv has no concordant row")

    sets = {(row["mark"], row["comparison"], row["region_type"], row["category"]): row
            for row in read_tsv(args.region_sets)}
    id_to_name = load_id_to_name(args.motif_dir)
    motifs = robust_motifs(args.background_compare, args.alpha, args.max_tp_pct, args.max_motifs)
    if not motifs:
        raise SystemExit("[tf-gene] no motif survived both backgrounds; nothing to scan")

    scanned = os.path.join(args.outdir, "fimo")
    os.makedirs(scanned, exist_ok=True)

    hit_rows: list[dict] = []
    edge_rows: list[dict] = []
    summary_rows: list[dict] = []
    for key, genes in sorted(common.items()):
        mark, comparison, region_type, category = key
        set_row = sets.get(key)
        if set_row is None:
            continue
        selected = motifs.get(set_row["set_id"], [])
        if not selected:
            continue
        regions = subset_regions(read_bed(set_row["bed"]), genes, args.max_regions, args.seed)
        if not regions:
            continue

        set_dir = os.path.join(scanned, set_row["set_id"])
        os.makedirs(set_dir, exist_ok=True)
        bed_path = os.path.join(set_dir, "windows.bed")
        index = write_window_bed(regions, args.window, bed_path)
        fasta_path = os.path.join(set_dir, "windows.fa")
        # -nameOnly, not -name: -name appends "::chrom:start-end" to every
        # header and MEME reads that suffix as sequence coordinates, replacing
        # the window id with the chromosome in its output.
        run(["bedtools", "getfasta", "-fi", args.genome, "-bed", bed_path,
             "-fo", fasta_path, "-nameOnly"], log_path)

        blocks = [block for block in
                  (read_meme_block(args.motif_dir, motif["motif_id"]) for motif in selected)
                  if block]
        if not blocks:
            continue
        meme_path = os.path.join(set_dir, "motifs.meme")
        write_meme(blocks, meme_path)

        fimo_table = os.path.join(set_dir, "fimo.tsv")
        fimo_out = os.path.join(set_dir, "fimo")
        os.makedirs(fimo_out, exist_ok=True)
        with open(fimo_table, "w") as out:
            print(f"[tf-gene] fimo {set_row['set_id']} {len(index)} sequences "
                  f"x {len(blocks)} motifs", flush=True)
            subprocess.run(
                ["fimo", "--text", "--thresh", str(args.thresh), "--verbosity", "1",
                 "--oc", fimo_out, meme_path, fasta_path],
                check=True, stdout=out, stderr=subprocess.DEVNULL,
            )

        q_of = {motif["motif_id"]: motif["q"] for motif in selected}
        alt_of = {motif["motif_id"]: motif["motif_alt_id"] for motif in selected}
        per_edge: dict[tuple, dict] = {}
        with open(fimo_table) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                motif_id = row.get("motif_id") or ""
                sequence = row.get("sequence_name") or ""
                region = region_of(index, sequence)
                if not motif_id or region is None:
                    continue
                tf = (id_to_name.get(motif_id) or alt_of.get(motif_id) or motif_id).strip()
                p_value = as_float(row.get("p-value"), 1.0)
                try:
                    hit_start = region["window_start"] + int(row.get("start", 0)) - 1
                    hit_stop = region["window_start"] + int(row.get("stop", 0))
                except (TypeError, ValueError):
                    hit_start, hit_stop = "", ""
                hit_rows.append(
                    {
                        "motif_id": motif_id, "motif_alt_id": alt_of.get(motif_id, ""),
                        "tf": tf, "gene": region["gene"], "set_id": set_row["set_id"],
                        "mark": mark, "comparison": comparison,
                        "region_type": region_type, "category": category,
                        "sequence": region["label"], "start": hit_start,
                        "stop": hit_stop, "strand": row.get("strand", ""),
                        "score": row.get("score", ""), "p": p_value,
                    }
                )
                edge_key = (tf, region["gene"])
                entry = per_edge.get(edge_key)
                if entry is None:
                    per_edge[edge_key] = {
                        "regions": {sequence}, "best_p": p_value,
                        "motif_id": motif_id, "q": q_of.get(motif_id, 1.0),
                    }
                else:
                    entry["regions"].add(sequence)
                    if p_value < entry["best_p"]:
                        entry["best_p"] = p_value
                        entry["motif_id"] = motif_id
                        entry["q"] = q_of.get(motif_id, 1.0)
                    else:
                        entry["q"] = min(entry["q"], q_of.get(motif_id, 1.0))

        gene_regions: dict[str, int] = defaultdict(int)
        for region in regions:
            gene_regions[region["gene"]] += 1
        for (tf, gene), entry in per_edge.items():
            edge_rows.append(
                {
                    "tf": tf, "gene": gene, "mark": mark, "comparison": comparison,
                    "region_type": region_type, "category": category,
                    "motif_id": entry["motif_id"],
                    "n_hit_regions": len(entry["regions"]),
                    "n_regions_gene": gene_regions[gene],
                    "best_p": entry["best_p"], "fdr": entry["q"],
                    "evidence_level": "region_hit",
                }
            )
        summary_rows.append(
            {
                "mark": mark, "comparison": comparison, "region_type": region_type,
                "category": category, "set_id": set_row["set_id"],
                "n_motifs": len(blocks), "n_sequences": len(index),
                "n_hits": sum(1 for row in hit_rows if row["set_id"] == set_row["set_id"]),
                "n_edges": len(per_edge),
                "n_tfs": len({tf for tf, _ in per_edge}),
                "n_genes": len({gene for _, gene in per_edge}),
            }
        )

    if not edge_rows:
        raise SystemExit("[tf-gene] no edge reached a common gene; check --window/--thresh")

    hit_columns = ["motif_id", "motif_alt_id", "tf", "gene", "set_id", "mark", "comparison",
                   "region_type", "category", "sequence", "start", "stop", "strand",
                   "score", "p"]
    with open(os.path.join(args.outdir, "fimo_hits.tsv"), "w") as out:
        out.write("\t".join(hit_columns) + "\n")
        for row in hit_rows:
            out.write("\t".join(str(row[column]) for column in hit_columns) + "\n")

    edge_rows.sort(key=lambda row: (row["mark"], row["comparison"], row["region_type"],
                                    row["category"], -row["n_hit_regions"], row["best_p"]))
    with open(os.path.join(args.outdir, "tf_gene_edges.tsv"), "w") as out:
        out.write("\t".join(EDGE_COLUMNS) + "\n")
        for row in edge_rows:
            out.write("\t".join(str(row[column]) for column in EDGE_COLUMNS) + "\n")

    summary_columns = ["mark", "comparison", "region_type", "category", "set_id",
                       "n_motifs", "n_sequences", "n_hits", "n_edges", "n_tfs", "n_genes"]
    with open(os.path.join(args.outdir, "tf_gene_summary.tsv"), "w") as out:
        out.write("\t".join(summary_columns) + "\n")
        for row in summary_rows:
            out.write("\t".join(str(row[column]) for column in summary_columns) + "\n")

    with open(os.path.join(args.outdir, "tf_gene_network.version"), "w") as out:
        out.write(f"fimo_thresh\t{args.thresh}\n")
        out.write(f"window\t{args.window}\n")
        out.write(f"alpha\t{args.alpha}\n")
        out.write(f"max_tp_pct\t{args.max_tp_pct}\n")
        out.write(f"max_motifs\t{args.max_motifs}\n")
        out.write(f"max_regions\t{args.max_regions}\n")
        out.write("evidence_level\tregion_hit\n")
        out.write("background_rule\tq < alpha under static and shuffle\n")

    # The log only receives failed commands otherwise, which leaves the file at
    # zero bytes and trips the acceptance check for placeholder outputs.
    with open(log_path, "a") as log:
        log.write(f"sets\t{len(summary_rows)}\n")
        log.write(f"hits\t{len(hit_rows)}\tsequences\t{sum(r['n_sequences'] for r in summary_rows)}\n")
        log.write(f"edges\t{len(edge_rows)}\n")

    print(
        f"[tf-gene] {len(hit_rows)} fimo hits, {len(edge_rows)} edges over "
        f"{len({row['tf'] for row in edge_rows})} TFs -> {args.outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
