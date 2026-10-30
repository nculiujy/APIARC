#!/usr/bin/env python3
"""Motif enrichment (AME) whose background answers the question being asked.

Two defects motivated this rewrite:

* The previous driver read ``ame.tsv``, kept only the first motif's adjusted
  p value and wrote it into every row of the output table, so the reported q
  value belonged to the top-ranked motif only.
* The background was ``--shuffle--`` alone, i.e. composition-unmatched random
  sequences. That tests "is this peak set enriched for the motif relative to
  random DNA", which is almost always yes. The question of interest is whether
  the motif is specific to *differential* binding, so the control is the
  static-bound peak set of the same mark, with a dinucleotide-preserving
  shuffle as the fallback when no static set exists.

Every set is now scored against both controls when a static set exists. A motif
that beats the static background only because the whole peak set is GC-rich also
beats shuffled DNA, so requiring both keeps the GC-driven matrices out of the
network; ``motif_background_compare.tsv`` keeps the two q values side by side.

Per-region evidence (fimo) and the TF-gene edges are built in
``tf_gene_network.py``: a set-level enrichment does not say which gene's peak
carries the site, and that is the only thing an edge should mean.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing
import os
import random
import shutil
import subprocess
from collections import defaultdict

COLUMNS = [
    "mark", "comparison", "region_type", "category", "set_id",
    "motif_id", "motif_alt_id", "consensus", "n_foreground", "n_background",
    "background_mode", "tp_pct", "p", "q",
]


def run(cmd: list[str], log=None) -> None:
    print("[motif] " + " ".join(cmd), flush=True)
    stream = open(log, "a") if log else None
    try:
        subprocess.run(cmd, check=True, stdout=stream, stderr=subprocess.STDOUT if stream else None)
    finally:
        if stream:
            stream.close()


def ame_table_ok(path: str) -> bool:
    """True when an AME table looks finished: one header, whole last row."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    with open(path, "rb") as handle:
        handle.seek(-1, os.SEEK_END)
        if handle.read(1) != b"\n":
            return False
    with open(path) as handle:
        lines = handle.readlines()
    if len(lines) < 2:
        return False
    return len(lines[-1].rstrip("\n").split("\t")) == len(lines[0].rstrip("\n").split("\t"))


def run_ame_task(task: dict) -> tuple[str, str, list[dict]]:
    """One AME run in a worker process.

    AME is single-threaded and scores every matrix against every sequence, so a
    peak set takes minutes and 34 runs in a row take hours. The sets do not
    depend on each other, which is the only parallelism available: the parent
    writes the FASTA files first because the two backgrounds of one set share
    them, and each run keeps its own table and log so nothing interleaves.

    ``--reuse-ame`` keeps a finished table from an earlier invocation; the
    aggregation below is cheap, the scan is not.
    """
    os.makedirs(task["set_dir"], exist_ok=True)
    table = os.path.join(task["set_dir"], f"ame_{task['mode']}.tsv")
    if not (task.get("reuse") and ame_table_ok(table)):
        with open(os.path.join(task["set_dir"], f"ame_{task['mode']}.log"), "a") as err:
            with open(table, "w") as out:
                subprocess.run(
                    ame_command(task["foreground"], task["database"], task["mode"],
                                task["background"], task["kmer"]),
                    check=True, stdout=out, stderr=err,
                )
    rows = []
    for row in read_ame(table):
        # AME reports TP/FP against primary and control, not FN: the share of
        # primary sequences carrying the motif is the %TP column. Without it the
        # fraction fell back to tp/(tp+0) = 1.0 and every matrix looked useless.
        pct = row["tp_pct"]
        if pct in ("", "NA"):
            tp = as_float(row["tp"], 0.0)
            fn = as_float(row["fn"], 0.0)
            total = tp + fn
            tp_pct = round(tp / total, 4) if total > 0 else ""
        else:
            tp_pct = round(as_float(pct, 0.0) / 100.0, 4)
        rows.append(
            {
                "mark": task["mark"],
                "comparison": task["comparison"],
                "region_type": task["region_type"],
                "category": task["category"],
                "set_id": task["set_id"],
                "motif_id": row["motif_id"],
                "motif_alt_id": row["motif_alt_id"],
                "consensus": row["consensus"],
                "n_foreground": task["n_foreground"],
                "n_background": task["n_background"] if task["mode"] == "static" else "",
                "background_mode": task["mode"],
                "tp_pct": tp_pct,
                "p": row["p"],
                "q": row["q"],
            }
        )
    return task["set_id"], task["mode"], rows


def read_sets(path: str) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_id_to_name(motif_dir: str) -> dict[str, str]:
    """JASPAR matrix id -> TF symbol, taken from the file beside the matrices."""
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


def species_tf_symbols(path: str) -> set[str]:
    symbols = set()
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            symbol = (row.get("Symbol") or "").strip().upper()
            if symbol:
                symbols.add(symbol)
    return symbols


def build_motif_database(motif_dir: str, out_path: str, allowed=None) -> int:
    """Merge the individual JASPAR MEME files into one searchable database.

    Concatenating the files would repeat the ``MEME version``/``ALPHABET`` header
    and produce an invalid file, so the header is written once and only the
    ``MOTIF`` blocks are appended.

    ``allowed`` keeps only matrices for TFs that are both in the species TF
    list and expressed here. AME scores every matrix against every sequence, so
    the full core collection (2346 matrices) put a single peak set at forty
    minutes; the discarded matrices could never reach the TF-gene edges either,
    because those keep expressed TFs only.
    """
    files = sorted(
        os.path.join(motif_dir, name)
        for name in os.listdir(motif_dir)
        if name.endswith(".meme")
    )
    if not files:
        raise SystemExit(f"no .meme files under {motif_dir}")
    n_available = len(files)
    if allowed is not None:
        id_to_name = load_id_to_name(motif_dir)
        files = [
            path for path in files
            if id_to_name.get(os.path.basename(path)[: -len(".meme")], "").upper() in allowed
        ]
        if not files:
            raise SystemExit(
                f"no matrix in {motif_dir} belongs to an expressed TF; check --tf-list"
            )
    blocks: list[str] = []
    for path in files:
        with open(path) as handle:
            lines = handle.readlines()
        start = next((i for i, line in enumerate(lines) if line.startswith("MOTIF")), None)
        if start is None:
            continue
        blocks.append("".join(lines[start:]).rstrip("\n") + "\n")
    with open(out_path, "w") as out:
        out.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n")
        out.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        out.writelines(blocks)
    with open(out_path + ".version", "w") as out:
        out.write(f"source_dir\t{motif_dir}\nn_matrices\t{len(blocks)}\n")
        out.write(f"n_candidates\t{n_available}\ntf_filter\t{'species+expressed' if allowed is not None else 'none'}\n")
        for path in files:
            out.write(f"matrix\t{os.path.basename(path)}\n")
    return len(blocks)


def ensure_fasta_index(genome: str) -> None:
    index = genome + ".fai"
    if not os.path.exists(index):
        print(f"[motif] building FASTA index for {genome}", flush=True)
        subprocess.run(["samtools", "faidx", genome], check=True)


def subset_bed(path: str, limit: int, seed: int, out_path: str) -> str:
    """Random subset of BED lines, with a fixed seed.

    AME scores every primary and control sequence against every matrix, so a
    130k-region control set turns one run into half an hour. Capping the
    control keeps the distance/mark-matched design while making the stage
    finish; the sampled lines are reproducible from the seed.
    """
    if limit <= 0:
        return path
    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if len(lines) <= limit:
        return path
    keep = sorted(random.Random(seed).sample(range(len(lines)), limit))
    with open(out_path, "w") as out:
        out.write("\n".join(lines[index] for index in keep) + "\n")
    return out_path


def bed_to_fasta(genome: str, bed: str, fasta: str, limit: int = 0, seed: int = 1) -> int:
    if limit:
        bed = subset_bed(bed, limit, seed, fasta + ".subset.bed")
    run(["bedtools", "getfasta", "-fi", genome, "-bed", bed, "-fo", fasta])
    with open(fasta) as handle:
        return sum(1 for line in handle if line.startswith(">"))


def read_ame(path: str) -> list[dict]:
    """Parse an AME result table by header name.

    Column order differs between MEME releases, and the previous implementation
    indexed by position, which is how the wrong column ended up as a q value.
    """
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            identifiers = {key.lower(): key for key in row}
            def pick(*names, default=""):
                for name in names:
                    key = identifiers.get(name)
                    if key is not None and row[key] not in (None, ""):
                        return row[key]
                return default
            rows.append(
                {
                    "motif_id": pick("motif_id", "motif"),
                    "motif_alt_id": pick("motif_alt_id", "motif_alt_id "),
                    "consensus": pick("consensus"),
                    "p": pick("p-value", "p_value", "pvalue"),
                    "q": pick("adj_p-value", "adj_p_value", "q-value", "qvalue"),
                    "tp": pick("tp", default="0"),
                    "fn": pick("fn", default="0"),
                    "tp_pct": pick("%tp", "percent_tp", "tp_pct"),
                    "fp": pick("fp", default="0"),
                }
            )
    return rows


def as_float(value, default: float = 1.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def ame_command(foreground: str, database: str, mode: str, background: str, kmer: int) -> list[str]:
    """``ame --text`` writes the result table to stdout and skips the HTML
    report. The report needs the MEME template directory, which this conda
    package does not ship, and AME aborted on it right after writing the table.
    ``ame`` has no --threads flag either - passing one made its parser read the
    next option as a motif file."""
    cmd = ["ame", "--text"]
    if mode == "static":
        cmd += ["--control", background]
    else:
        cmd += ["--control", "--shuffle--", "--kmer", str(kmer)]
    return cmd + [foreground, database]


def main() -> None:
    parser = argparse.ArgumentParser(description="motif enrichment with a matched background")
    parser.add_argument("--region-sets", required=True)
    parser.add_argument("--rna", action="append", default=[],
                        help="DESeq2 tables, one per contrast; the TF filter uses their union")
    parser.add_argument("--motif-dir", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--background", default="both", choices=["both", "static", "shuffle"])
    parser.add_argument("--kmer", type=int, default=2)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-regions", type=int, default=50)
    parser.add_argument("--min-base-mean", type=float, default=10.0)
    parser.add_argument("--max-background", type=int, default=20000)
    parser.add_argument("--tf-list", default="",
                        help="species TF table; keeps only matrices of expressed TFs")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1,
                        help="worker processes; one AME run per worker")
    parser.add_argument("--reuse-ame", action="store_true",
                        help="keep finished ame_*.tsv files and only re-aggregate")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    if shutil.which("ame") is None:
        raise SystemExit("ame not found on PATH; MEME suite is a hard requirement")
    if shutil.which("bedtools") is None:
        raise SystemExit("bedtools not found on PATH")

    expressed = set()
    for path in args.rna:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                symbol = (row.get("SYMBOL") or "").strip()
                try:
                    base_mean = float(row.get("baseMean") or 0)
                except ValueError:
                    continue
                if symbol and base_mean >= args.min_base_mean:
                    expressed.add(symbol.upper())
    if not expressed:
        raise SystemExit("[motif] no --rna table produced an expressed gene list")

    allowed = None
    if args.tf_list:
        allowed = species_tf_symbols(args.tf_list) & expressed
        print(
            f"[motif] {len(allowed)} expressed TFs in {os.path.basename(args.tf_list)}",
            flush=True,
        )

    os.makedirs(args.outdir, exist_ok=True)
    log_path = os.path.join(args.outdir, "motif.log")
    database = os.path.join(args.outdir, "motif_database.meme")
    n_matrices = build_motif_database(args.motif_dir, database, allowed)
    print(f"[motif] {n_matrices} matrices kept for scanning", flush=True)
    ensure_fasta_index(args.genome)

    sets = read_sets(args.region_sets)
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in sets:
        grouped[(row["mark"], row["comparison"], row["region_type"])].append(row)

    modes = ["static", "shuffle"] if args.background == "both" else [args.background]

    # Inputs are prepared here: the foreground FASTA is shared by the two
    # backgrounds of one set, so writing it from workers would race.
    tasks: list[dict] = []
    for (mark, comparison, region_type), members in sorted(grouped.items()):
        static = next((m for m in members if m["category"] == "static"), None)
        background_fasta = None
        n_background = 0
        if static is not None and int(static["n_regions"]) >= args.min_regions:
            background_fasta = os.path.join(args.outdir, f"{static['set_id']}__background.fa")
            n_background = bed_to_fasta(
                args.genome, static["bed"], background_fasta,
                args.max_background, args.seed,
            )
        if background_fasta is None and "static" in modes:
            if args.background == "static":
                print(
                    f"[motif] {mark} {comparison} {region_type} has no usable static set; "
                    "skipping while --background is static",
                    flush=True,
                )
                continue
            modes_for_set = ["shuffle"]
        else:
            modes_for_set = modes

        for member in members:
            if member["category"] == "static":
                continue
            if int(member["n_regions"]) < args.min_regions:
                continue
            set_dir = os.path.join(args.outdir, member["set_id"])
            os.makedirs(set_dir, exist_ok=True)
            foreground_fasta = os.path.join(set_dir, "foreground.fa")
            n_foreground = bed_to_fasta(args.genome, member["bed"], foreground_fasta)
            if n_foreground < args.min_regions:
                continue
            for mode in modes_for_set:
                tasks.append(
                    {
                        "set_dir": set_dir, "foreground": foreground_fasta,
                        "database": database, "mode": mode,
                        "background": background_fasta, "kmer": args.kmer,
                        "mark": mark, "comparison": comparison,
                        "region_type": region_type, "category": member["category"],
                        "set_id": member["set_id"], "n_foreground": n_foreground,
                        "n_background": n_background, "reuse": args.reuse_ame,
                    }
                )

    if not tasks:
        raise SystemExit("[motif] no region set reached --min-regions")
    with open(log_path, "a") as log:
        log.write(f"tasks\t{len(tasks)}\nprocesses\t{args.threads}\n")
    print(f"[motif] {len(tasks)} AME run(s) over {args.threads} process(es)", flush=True)

    hits: list[dict] = []
    with multiprocessing.Pool(processes=max(1, args.threads)) as pool:
        for _set_id, _mode, rows in pool.imap_unordered(run_ame_task, tasks):
            hits.extend(rows)

    if not hits:
        raise SystemExit("[motif] no motif result was produced; check the region sets")

    with open(os.path.join(args.outdir, "motif_hits.tsv"), "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in hits:
            out.write("\t".join(str(row[column]) for column in COLUMNS) + "\n")

    # The two backgrounds side by side. A motif worth reporting beats both: it
    # is not merely enriched relative to static binding, and not merely
    # enriched relative to random DNA of the same composition.
    by_set: dict[tuple, dict[str, dict]] = defaultdict(dict)
    for row in hits:
        by_set[(row["set_id"], row["motif_id"])][row["background_mode"]] = row
    compare_path = os.path.join(args.outdir, "motif_background_compare.tsv")
    with open(compare_path, "w") as out:
        out.write(
            "set_id\tmark\tcomparison\tregion_type\tcategory\tmotif_id\tmotif_alt_id\t"
            "q_static\tq_shuffle\ttp_pct_static\ttp_pct_shuffle\trobust\n"
        )
        for (set_id, motif_id), modes_present in sorted(by_set.items()):
            static_row = modes_present.get("static")
            shuffle_row = modes_present.get("shuffle")
            reference = static_row or shuffle_row
            q_static = as_float(static_row["q"], 1.0) if static_row else ""
            q_shuffle = as_float(shuffle_row["q"], 1.0) if shuffle_row else ""
            margins = [value for value in (q_static, q_shuffle) if value != ""]
            robust = bool(margins) and max(margins) < args.alpha
            out.write(
                "\t".join(
                    str(value) for value in [
                        set_id, reference["mark"], reference["comparison"],
                        reference["region_type"], reference["category"], motif_id,
                        reference["motif_alt_id"], q_static, q_shuffle,
                        static_row["tp_pct"] if static_row else "",
                        shuffle_row["tp_pct"] if shuffle_row else "",
                        str(robust).lower(),
                    ]
                )
                + "\n"
            )

    summary: dict[tuple, dict] = {}
    for row in hits:
        key = (row["mark"], row["comparison"], row["region_type"], row["category"],
               row["background_mode"])
        entry = summary.setdefault(key, {"tested": 0, "significant": 0})
        entry["tested"] += 1
        if as_float(row["q"], 1.0) < args.alpha:
            entry["significant"] += 1
    with open(os.path.join(args.outdir, "motif_summary.tsv"), "w") as out:
        out.write("mark\tcomparison\tregion_type\tcategory\tbackground_mode\t"
                  "n_motifs_tested\tn_motifs_q05\n")
        for key, entry in sorted(summary.items()):
            out.write("\t".join([*key, str(entry["tested"]), str(entry["significant"])]) + "\n")

    print(
        f"[motif] {n_matrices} matrices, {len(hits)} motif/set rows over "
        f"{len(modes)} background(s) -> {args.outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
