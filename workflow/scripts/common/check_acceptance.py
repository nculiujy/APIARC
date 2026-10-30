#!/usr/bin/env python3
"""Post-run acceptance checks.

The refactor has one goal that is not about features: a claim must not be able
to reach a manuscript without a test, a sample size and an adjusted p value
behind it. These checks encode the acceptance criteria as code so that "the
pipeline works" is a non-zero exit status rather than a judgement call.

Checks (numbering follows the delivery plan):

  1. placeholders    no zero-byte or placeholder outputs anywhere under result/
  2. negative_control  the label-permutation / Input-vs-Input controls must pass
  3. traceability    every row of stats_summary.tsv carries n1 and, when it has
                     a p value, an FDR
  4. mode            descriptive runs must not leave a single diff_* artefact;
                     differential runs must produce at least one
  5. design          the design validation report must be free of errors
  6. layout          result/ contains no directory the current graph cannot
                     produce, i.e. no mixture of two pipeline versions
  7. common_genes    the common-gene table, its columns and its four figures
  8. common_enrichment  both enrichment tables, and no negative odds ratio
  9. common_network  the term-gene network tables and their figure
 10. motif_network   per-region TF-gene edges that beat both backgrounds

Run after a real (not dry) run::

    python workflow/scripts/common/check_acceptance.py

Exit status is 0 only when every check passes. A check that cannot be
evaluated (missing table, unreadable file) is a failure: an absent control is
not a passing control.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common import annotation, workspace
from common.metadata import read_chip_metadata, resolve_mode

STATS_COLUMNS = [
    "analysis", "comparison", "test", "n1", "n2", "statistic",
    "effect_size", "ci_low", "ci_high", "p", "FDR", "note",
]

PLACEHOLDER_MARKERS = ("Peak file not found", "placeholder")

EMPTY_BY_DESIGN = (".snakemake_timestamp",)

SELF_REPORTED = ("acceptance_report.tsv",)


def load_config(path: str) -> dict:
    import yaml

    with open(path) as handle:
        return yaml.safe_load(handle) or {}


def expected_tree(config: dict) -> dict:
    """Layout the current graph produces, annotation segment included.

    A node is either a dict of child directories or a set of allowed names;
    ``None`` stops the descent (provenance keeps free-form files).
    """
    rna = workspace.default_name(config, "RNAseq")
    chip = workspace.default_name(config, "ChIPseq")
    joint = workspace.integrated_id(rna, chip)
    anno_id = annotation.annotation_id(*annotation.declared(config))
    return {
        "RNAseq": {rna: {anno_id: set(workspace.STAGE_DIRS["RNAseq"].values())}},
        "ChIPseq": {chip: {anno_id: set(workspace.STAGE_DIRS["ChIPseq"].values())}},
        "Integrated": {joint: {anno_id: set(workspace.INTEGRATED_STAGES.values())}},
        "provenance": {joint: {anno_id: None}},
    }


class Report:
    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str]] = []

    def add(self, check: str, status: str, detail: str) -> None:
        self.rows.append((check, status, detail))
        print(f"[{status}] {check}: {detail}", flush=True)

    @property
    def failed(self) -> list[tuple[str, str, str]]:
        return [row for row in self.rows if row[1] == "FAIL"]

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as handle:
            handle.write("check\tstatus\tdetail\n")
            for check, status, detail in self.rows:
                handle.write(f"{check}\t{status}\t{detail}\n")


def stale_prefixes(root: Path, tree: dict) -> set[tuple[str, ...]]:
    """Directories in ``result/`` that the current layout does not define."""
    found: set[tuple[str, ...]] = set()
    if not root.exists():
        return found

    def walk(node, base: Path, prefix: tuple[str, ...]) -> None:
        for child in sorted(base.iterdir()):
            if not child.is_dir():
                continue
            if isinstance(node, set):
                if child.name not in node:
                    found.add(prefix + (child.name,))
                continue
            if child.name not in node:
                found.add(prefix + (child.name,))
                continue
            child_node = node[child.name]
            if isinstance(child_node, (dict, set)):
                walk(child_node, child, prefix + (child.name,))

    walk(tree, root, ())
    return found


def check_layout(root: Path, tree: dict, report: Report) -> None:
    """Reject results left behind by a previous version of the pipeline."""
    if not root.exists():
        report.add("layout", "FAIL", f"{root} does not exist; nothing was produced")
        return
    stale = ["/".join(prefix) for prefix in sorted(stale_prefixes(root, tree))]
    if stale:
        report.add(
            "layout",
            "FAIL",
            "output from a previous pipeline version (rerun into a clean result dir): "
            + ", ".join(stale[:8]),
        )
    else:
        report.add("layout", "PASS", "result tree matches the current module layout")


def check_placeholders(root: Path, tree: dict, report: Report) -> None:
    empty = []
    marked = []
    stale = stale_prefixes(root, tree)
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        if any(rel.parts[: len(prefix)] == prefix for prefix in stale):
            continue
        if path.name in EMPTY_BY_DESIGN or path.name in SELF_REPORTED:
            continue
        if path.name.endswith("finished.txt"):
            continue
        if path.stat().st_size == 0:
            empty.append(str(rel))
            continue
        if path.suffix in {".bw", ".bigwig", ".bam", ".gz", ".rds"}:
            continue
        try:
            head = path.read_text(errors="ignore")[:4096]
        except OSError:
            continue
        hit = [marker for marker in PLACEHOLDER_MARKERS if marker in head]
        if hit:
            marked.append(f"{rel} ({hit[0]})")

    if not root.exists():
        report.add("placeholders", "FAIL", f"{root} does not exist; nothing was produced")
        return
    if empty:
        report.add("placeholders", "FAIL", f"zero-byte outputs: {', '.join(empty[:8])}")
    elif marked:
        report.add("placeholders", "FAIL", f"placeholder content: {', '.join(marked[:8])}")
    else:
        report.add("placeholders", "PASS", "no zero-byte or placeholder outputs")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))
            rows.append(dict(zip(header, fields)))
    return rows


def check_negative_control(root: Path, chip_dir: Path, mode: str, report: Report) -> None:
    tables = sorted(chip_dir.glob("5_diffbind/*/*/negative_control.tsv"))
    if mode == "descriptive":
        if tables:
            report.add(
                "negative_control",
                "FAIL",
                f"descriptive run wrote {len(tables)} negative_control.tsv file(s); "
                "a run without a contrast has nothing to control",
            )
        else:
            report.add("negative_control", "PASS", "descriptive mode: no contrast, no control needed")
        return

    if not tables:
        report.add("negative_control", "FAIL", "differential mode but no negative_control.tsv found")
        return

    problems = []
    for table in tables:
        label = str(table.relative_to(root))
        for row in _read_tsv(table):
            control = row.get("control", "?")
            status = row.get("status", "?")
            if status not in {"PASS", "SKIPPED"}:
                problems.append(f"{label}:{control}={status}")
                continue
            if status == "SKIPPED":
                problems.append(f"{label}:{control}=SKIPPED ({row.get('note', '')})")
                continue
            value = row.get("max_false_positive", "NA")
            threshold = row.get("threshold", "NA")
            try:
                if float(value) > float(threshold):
                    problems.append(f"{label}:{control}={value} > {threshold}")
            except ValueError:
                problems.append(f"{label}:{control}=unparsable max_false_positive '{value}'")

    hard = [item for item in problems if "=SKIPPED" not in item]
    soft = [item for item in problems if "=SKIPPED" in item]
    if hard:
        report.add("negative_control", "FAIL", "; ".join(hard[:6]))
    elif soft:
        report.add("negative_control", "PASS", f"no failures; skipped: {'; '.join(soft[:4])}")
    else:
        report.add("negative_control", "PASS", f"{len(tables)} control table(s), all PASS")


def check_traceability(root: Path, joint_dir: Path, report: Report) -> None:
    path = joint_dir / "0_stats" / "stats_summary.tsv"
    if not path.exists():
        report.add("traceability", "FAIL", f"{path} missing; no claim can be traced")
        return

    with path.open(newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header != STATS_COLUMNS:
        report.add("traceability", "FAIL", f"unexpected columns: {header}")
        return

    rows = _read_tsv(path)
    if not rows:
        report.add("traceability", "FAIL", "stats_summary.tsv has a header but no rows")
        return

    missing_n = []
    missing_fdr = []
    tested = 0
    for row in rows:
        analysis = row.get("analysis", "?")
        if analysis == "negative_control":
            continue
        if not row.get("n1") or row.get("n1") == "NA":
            missing_n.append(f"{analysis}/{row.get('test', '?')}")
        p = row.get("p", "NA")
        if p and p != "NA":
            tested += 1
            fdr = row.get("FDR", "NA")
            if not fdr or fdr == "NA":
                missing_fdr.append(f"{analysis}/{row.get('test', '?')}")

    problems = []
    if missing_n:
        problems.append(f"rows without n1: {', '.join(sorted(set(missing_n))[:6])}")
    if missing_fdr:
        problems.append(f"rows with p but no FDR: {', '.join(sorted(set(missing_fdr))[:6])}")
    if problems:
        report.add("traceability", "FAIL", "; ".join(problems))
    else:
        report.add(
            "traceability",
            "PASS",
            f"{len(rows)} rows, {tested} with a p value and a family-wise FDR",
        )


def check_common(config: dict, joint_dir: Path, mode: str, report: Report) -> None:
    """The common-gene stage, when the graph is supposed to have produced it.

    Older result sets predate the stage, so absence only fails when the stage is
    switched on and the run is differential. When it does exist, the four
    figures and both direction sets must be present: a table without the plots
    is the state the six manuscript requirements cannot use.
    """
    enabled = bool((config.get("Integrated_modules") or {}).get("4_common", True))
    stage = joint_dir / "4_common"
    if mode != "differential" or not enabled:
        report.add("common_genes", "PASS", "not applicable: stage disabled or descriptive mode")
        return
    path = stage / "tables" / "common_genes.tsv"
    if not path.is_file():
        report.add("common_genes", "FAIL", f"{path} missing; the common-gene stage did not run")
        return
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)
    required = [
        "gene", "mark", "comparison", "region_type", "chip_direction",
        "concordant", "chip_lfc_mean", "rna_log2FC", "rna_padj",
    ]
    missing = [name for name in required if name not in header]
    if missing:
        report.add("common_genes", "FAIL", f"missing columns: {', '.join(missing)}")
        return
    if not rows:
        report.add("common_genes", "FAIL", "common_genes.tsv has a header but no rows")
        return
    figures = (
        "scatter_change.pdf",
        "line_expression.pdf",
        "heatmap_expression.pdf",
        "box_signal.pdf",
    )
    absent = [name for name in figures if not (stage / "plots" / name).is_file()]
    if absent:
        report.add("common_genes", "FAIL", f"missing figures: {', '.join(absent)}")
        return
    concordant = sum(1 for row in rows if row.get("concordant") == "true")
    report.add(
        "common_genes",
        "PASS",
        f"{len(rows)} gene rows ({concordant} concordant), {len(figures)} figures",
    )


def _stage_enabled(config: dict, name: str) -> bool:
    return bool((config.get("Integrated_modules") or {}).get(name, True))


def check_common_enrichment(config: dict, joint_dir: Path, mode: str, report: Report) -> None:
    """Both enrichment tables, and no odds ratio the four-fold table cannot give.

    A negative odds ratio only appears when the cells come from different
    denominators, so the value doubles as a regression test for the universe
    fix that the tables were rebuilt for.
    """
    if mode != "differential" or not _stage_enabled(config, "4_common"):
        report.add("common_enrichment", "PASS", "not applicable: stage disabled or descriptive mode")
        return
    tables = joint_dir / "4_common" / "tables"
    common = tables / "go_kegg_common.tsv"
    reference = tables / "go_kegg_all_diff.tsv"
    absent = [path.name for path in (common, reference) if not path.is_file()]
    if absent:
        report.add("common_enrichment", "FAIL", f"missing tables: {', '.join(absent)}")
        return
    rows = _read_tsv(common)
    if not rows:
        report.add("common_enrichment", "FAIL", "go_kegg_common.tsv has no rows")
        return
    negative = 0
    for row in rows:
        try:
            value = float(row.get("odds_ratio") or "")
        except ValueError:
            continue
        if value < 0:
            negative += 1
    if negative:
        report.add("common_enrichment", "FAIL", f"{negative} row(s) carry a negative odds ratio")
        return
    report.add("common_enrichment", "PASS", f"{len(rows)} term rows, no negative odds ratio")


def check_common_network(config: dict, joint_dir: Path, mode: str, report: Report) -> None:
    """The bipartite term-gene tables and their figure."""
    if mode != "differential" or not _stage_enabled(config, "4_common"):
        report.add("common_network", "PASS", "not applicable: stage disabled or descriptive mode")
        return
    stage = joint_dir / "4_common"
    edges = stage / "tables" / "network_edges.tsv"
    nodes = stage / "tables" / "network_nodes.tsv"
    figure = stage / "plots" / "network.pdf"
    absent = [path.name for path in (edges, nodes, figure) if not path.is_file()]
    if absent:
        report.add("common_network", "FAIL", f"missing: {', '.join(absent)}")
        return
    edge_rows = _read_tsv(edges)
    node_rows = _read_tsv(nodes)
    if not edge_rows or not node_rows:
        report.add("common_network", "FAIL", "network tables are empty")
        return
    kinds = {row.get("type") for row in node_rows}
    if not {"gene", "term"} <= kinds:
        report.add("common_network", "FAIL", "nodes lack gene or term rows")
        return
    report.add("common_network", "PASS", f"{len(edge_rows)} edges over {len(node_rows)} nodes")


def check_motif(config: dict, joint_dir: Path, mode: str, report: Report) -> None:
    """Per-region motif evidence, not set-level annotation.

    The previous network joined every TF enriched in a set to every gene of
    that set: millions of edges that say nothing about a particular peak. The
    row count is a guard as well, because that shape cannot climb back under a
    per-region scan.
    """
    if mode != "differential" or not _stage_enabled(config, "2_motif"):
        report.add("motif_network", "PASS", "not applicable: stage disabled or descriptive mode")
        return
    stage = joint_dir / "2_motif"
    compare = stage / "motif_background_compare.tsv"
    edges = stage / "tf_gene_edges.tsv"
    figure = stage / "plots" / "tf_gene_network.pdf"
    absent = [path.name for path in (compare, edges, figure) if not path.is_file()]
    if absent:
        report.add("motif_network", "FAIL", f"missing: {', '.join(absent)}")
        return
    robust = sum(1 for row in _read_tsv(compare) if row.get("robust") == "true")
    if robust == 0:
        report.add(
            "motif_network", "FAIL",
            "no matrix beats both the static and the shuffled background",
        )
        return
    rows = _read_tsv(edges)
    if not rows:
        report.add("motif_network", "FAIL", "tf_gene_edges.tsv has no rows")
        return
    required = {"tf", "gene", "n_hit_regions", "n_regions_gene", "best_p", "evidence_level"}
    missing = sorted(required - set(rows[0]))
    if missing:
        report.add("motif_network", "FAIL", f"missing columns: {', '.join(missing)}")
        return
    if len(rows) > 200_000:
        report.add(
            "motif_network", "FAIL",
            f"{len(rows)} edges; set-level annotation is back (per-region hits only)",
        )
        return
    report.add(
        "motif_network", "PASS",
        f"{len(rows)} per-region edges, {len({row['tf'] for row in rows})} TFs, "
        f"{robust} matrix/set pairs robust to both backgrounds",
    )


def check_mode(root: Path, mode: str, report: Report) -> None:
    differential = sorted(root.glob("**/differential_regions.tsv"))
    diff_named = [p for p in root.rglob("diff*") if p.is_file()] if root.exists() else []

    if mode == "descriptive":
        offenders = sorted({str(p.relative_to(root)) for p in diff_named} | {
            str(p.relative_to(root)) for p in differential
        })
        if offenders:
            report.add(
                "mode",
                "FAIL",
                f"descriptive run produced differential artefacts: {', '.join(offenders[:6])}",
            )
        else:
            report.add("mode", "PASS", "descriptive: no diff_* artefact on disk")
        return

    if not differential:
        report.add("mode", "FAIL", "differential mode but no differential_regions.tsv produced")
        return
    report.add("mode", "PASS", f"{len(differential)} differential_regions.tsv table(s)")


def check_design(root: Path, joint_id: str, anno_id: str, report: Report) -> None:
    path = root / "provenance" / joint_id / anno_id / "validation_report.tsv"
    if not path.exists():
        report.add("design", "FAIL", f"{path} missing; design was never validated")
        return
    rows = _read_tsv(path)
    errors = [row for row in rows if row.get("level") == "error"]
    warnings = [row for row in rows if row.get("level") == "warning"]
    if errors:
        first = errors[0]
        report.add(
            "design",
            "FAIL",
            f"{len(errors)} design error(s), first: {first.get('check', '?')}: "
            f"{first.get('message', '')}",
        )
    else:
        report.add("design", "PASS", f"no design errors ({len(warnings)} warning(s))")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", default="result", help="result root")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--metadata", default=None, help="ChIP metadata; default comes from --config")
    parser.add_argument("--mode", default="auto", choices=["auto", "differential", "descriptive"])
    parser.add_argument(
        "--report", default=None,
        help="where to write the machine-readable summary",
    )
    args = parser.parse_args()

    root = Path(args.result)
    config = load_config(args.config)
    tree = expected_tree(config)
    rna = workspace.default_name(config, "RNAseq")
    chip = workspace.default_name(config, "ChIPseq")
    joint = workspace.integrated_id(rna, chip)
    anno_id = annotation.annotation_id(*annotation.declared(config))
    chip_dir = workspace.result_dir("ChIPseq", chip, root, anno_id)
    joint_dir = workspace.integrated_dir(rna, chip, root, anno_id)
    metadata = Path(args.metadata) if args.metadata else Path(
        workspace.metadata_from_config(config, "ChIPseq")
    )
    report_path = (
        Path(args.report)
        if args.report
        else root / "provenance" / joint / anno_id / "acceptance_report.tsv"
    )

    report = Report()

    mode = args.mode
    if mode == "auto":
        try:
            samples = read_chip_metadata(metadata)
            mode = resolve_mode({}, samples)
        except (OSError, ValueError) as error:
            report.add("mode", "FAIL", f"cannot resolve mode: {error}")
            report.write(report_path)
            raise SystemExit(1)

    print(f"[info] resolved mode: {mode}", flush=True)
    print(f"[info] dataset: {rna} x {chip} -> {joint} [{anno_id}]", flush=True)

    check_design(root, joint, anno_id, report)
    check_layout(root, tree, report)
    check_mode(root, mode, report)
    check_negative_control(root, chip_dir, mode, report)
    check_traceability(root, joint_dir, report)
    check_common(config, joint_dir, mode, report)
    check_common_enrichment(config, joint_dir, mode, report)
    check_common_network(config, joint_dir, mode, report)
    check_motif(config, joint_dir, mode, report)
    check_placeholders(root, tree, report)

    report.write(report_path)
    failed = report.failed
    if failed:
        names = ", ".join(check for check, _status, _detail in failed)
        raise SystemExit(f"[acceptance] {len(failed)} check(s) failed: {names}")
    print("[acceptance] all checks passed", flush=True)


if __name__ == "__main__":
    main()
