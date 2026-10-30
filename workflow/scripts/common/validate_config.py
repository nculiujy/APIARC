#!/usr/bin/env python3
"""Fail fast when the project metadata cannot support the requested analysis.

Everything checked here is something that previously failed silently or failed
much later with an unhelpful message:

* a mark whose replicates are not all paired to the same Input;
* a mark with a single replicate, which makes a two-group comparison
  unidentifiable;
* a requested differential analysis on a single-condition design;
* RNA-seq metadata whose sample names do not appear in the count matrix;
* a declared reference module whose output files are absent.

The script exits non-zero with a list of problems. It writes a machine readable
report so the same checks can be asserted in CI.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common import annotation, workspace
from common.metadata import (
    read_chip_inputs,
    read_chip_metadata,
    read_comparisons,
    uses_group_names,
)


def read_rna(path: str) -> list[dict]:
    if not os.path.exists(path):
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def check_chip(
    samples, comparisons_path: str = "", input_rows=None
) -> list[tuple[str, str, str]]:
    problems = []
    if not samples:
        problems.append(("chip_metadata", "error", "no ChIP-seq samples declared"))
        return problems
    if uses_group_names(samples):
        return check_chip_groups(samples, comparisons_path, input_rows)
    return check_chip_marks(samples)


def check_chip_groups(
    samples, comparisons_path: str, input_rows=None
) -> list[tuple[str, str, str]]:
    """Design checks for the role-aware schema.

    Each ``sample_name`` is one group and ``condition`` only numbers its
    replicates; comparisons are declared in ``comparisons.csv`` instead of
    being derived from mark/condition pairs. ``input_rows`` is ``None`` for
    files that predate the role column, which cannot declare an Input row.
    """
    problems = []
    by_group: dict[str, list] = {}
    for sample in samples:
        by_group.setdefault(sample.sample_name, []).append(sample)

    declared = {}
    if input_rows is not None:
        declared = {item.sample: item.sample_name for item in input_rows}
        cited = {sample.input_id for sample in samples}
        for sample in samples:
            if sample.input_id not in declared:
                problems.append((
                    "chip_metadata", "error",
                    f"{sample.sample}: input_sample '{sample.input_id}' "
                    "has no row with role Input",
                ))
        for input_id in sorted(set(declared) - cited):
            problems.append((
                "chip_design", "info",
                f"Input {input_id} ({declared[input_id]}) is not cited by any IP row",
            ))

    for group, members in sorted(by_group.items()):
        replicates = sorted(s.replicate for s in members)
        if len(set(replicates)) != len(replicates):
            problems.append((
                "chip_metadata", "error",
                f"{group}: duplicate replicate indices {replicates}",
            ))
        if min(replicates) < 1:
            problems.append((
                "chip_metadata", "error",
                f"{group}: replicate indices must start at 1, found {replicates}",
            ))
        if len(members) < 2:
            problems.append((
                "chip_design", "error",
                f"{group}: {len(members)} replicate(s); a comparison needs at least 2",
            ))
        used_inputs = sorted({s.input_id for s in members})
        if not used_inputs or "" in used_inputs:
            problems.append((
                "chip_metadata", "error",
                f"{group}: every IP needs an input_sample",
            ))
        elif len(used_inputs) > 1:
            problems.append((
                "chip_design", "info",
                f"{group}: {len(used_inputs)} Input libraries inside one group; "
                "each replicate is paired to its own Input",
            ))

    pairs = read_comparisons(comparisons_path) if comparisons_path else []
    if not pairs:
        problems.append((
            "chip_comparisons", "warning",
            "no comparison declared; only descriptive outputs are possible",
        ))
    for reference, treatment in pairs:
        for label in (reference, treatment):
            if label not in by_group:
                problems.append((
                    "chip_comparisons", "error",
                    f"{treatment}_vs_{reference}: unknown sample_name '{label}'",
                ))
            elif len(by_group[label]) < 2:
                problems.append((
                    "chip_comparisons", "error",
                    f"{treatment}_vs_{reference}: '{label}' has "
                    f"{len(by_group[label])} replicate(s), need at least 2",
                ))
        shared = {s.sample for s in by_group.get(reference, [])} & {
            s.sample for s in by_group.get(treatment, [])
        }
        if shared:
            problems.append((
                "chip_comparisons", "error",
                f"{treatment}_vs_{reference}: share libraries {sorted(shared)}",
            ))
    return problems


def check_chip_marks(samples) -> list[tuple[str, str, str]]:
    """Design checks for the legacy mark/condition schema."""
    problems = []
    by_mark: dict[str, list] = {}
    for sample in samples:
        by_mark.setdefault(sample.mark, []).append(sample)

    for mark, members in sorted(by_mark.items()):
        # Replicate numbering restarts in every condition, so 1,2 per condition
        # is a complete between-condition design and not a duplicate.
        for condition in sorted({s.condition for s in members}):
            replicates = sorted(s.replicate for s in members if s.condition == condition)
            if len(set(replicates)) != len(replicates):
                problems.append((
                    "chip_metadata", "error",
                    f"{mark}/{condition}: duplicate replicate indices {replicates}",
                ))
            if min(replicates) < 1:
                problems.append((
                    "chip_metadata", "error",
                    f"{mark}/{condition}: replicate indices must start at 1, found {replicates}",
                ))

        conditions = sorted({s.condition for s in members})
        inputs = sorted({s.input_id for s in members})
        if len(conditions) > 1:
            per_condition_inputs = {
                condition: {s.input_id for s in members if s.condition == condition}
                for condition in conditions
            }
            shared = set.intersection(*per_condition_inputs.values())
            if shared and len(inputs) < len(conditions):
                problems.append((
                    "chip_design", "error",
                    f"{mark}: conditions {conditions} share Input(s) {sorted(shared)}; "
                    "each condition needs its own Input",
                ))
            for condition, condition_inputs in per_condition_inputs.items():
                n_reps = len([s for s in members if s.condition == condition])
                if n_reps < 2:
                    problems.append((
                        "chip_design", "error",
                        f"{mark}/{condition}: {n_reps} replicate(s); a between-condition "
                        "test needs at least 2",
                    ))
        else:
            problems.append((
                "chip_design", "warning",
                f"{mark}: single condition ({conditions[0]}); only descriptive "
                "outputs are possible for this mark",
            ))
        if len(inputs) > len(conditions) and len(conditions) == 1:
            problems.append((
                "chip_design", "info",
                f"{mark}: {len(inputs)} Input libraries for {len(members)} IPs",
            ))
    return problems


def check_rna(rna: list[dict], count_matrix: str) -> list[tuple[str, str, str]]:
    problems = []
    if not rna:
        problems.append(("rna_metadata", "error", "no RNA-seq samples declared"))
        return problems

    groups = {}
    for row in rna:
        groups.setdefault((row.get("group") or "").strip(), []).append(
            row.get("sample_name", "")
        )
    unknown = sorted(set(groups) - set(workspace.RNA_GROUPS))
    if unknown:
        problems.append((
            "rna_design", "error",
            f"unknown RNA group(s) {unknown}; expected control or treatA-D",
        ))
    if "control" not in groups:
        problems.append((
            "rna_design", "error",
            "RNA-seq needs a 'control' group as the DESeq2 reference",
        ))
    if not set(groups) - {"control", ""}:
        problems.append((
            "rna_design", "error",
            "RNA-seq needs at least one treat group besides control",
        ))
    for group, members in sorted(groups.items()):
        if len(members) < 2:
            problems.append((
                "rna_design", "error",
                f"RNA group '{group}' has {len(members)} replicate(s); need at least 2",
            ))

    if os.path.exists(count_matrix):
        with open(count_matrix, newline="") as handle:
            header = next(csv.reader(handle), [])
        columns = {name.strip().strip('"') for name in header[1:]}
        # Replicates of one group share a sample_name, so only the run id can
        # appear in the count matrix; DESeq2 keys its colData the same way.
        # Comparing sample_name here rejected every valid metadata table whose
        # sample_name is a condition label.
        key = "sample" if any("sample" in row for row in rna) else "sample_name"
        declared = {(row.get(key) or "").strip() for row in rna}
        missing = sorted(declared - columns)
        if missing:
            problems.append((
                "rna_metadata", "error",
                f"RNA samples missing from {count_matrix}: {', '.join(missing)}",
            ))
    else:
        problems.append((
            "rna_matrix", "warning",
            f"count matrix not found yet: {count_matrix}",
        ))
    return problems


def check_mode(
    samples, requested_mode: str, comparisons_path: str = ""
) -> list[tuple[str, str, str]]:
    if uses_group_names(samples):
        differential = bool(read_comparisons(comparisons_path)) if comparisons_path else False
        names = sorted({s.sample_name for s in samples})
    else:
        conditions_by_mark: dict[str, set] = {}
        for sample in samples:
            conditions_by_mark.setdefault(sample.mark, set()).add(sample.condition)
        differential = any(len(c) > 1 for c in conditions_by_mark.values())
        names = sorted(conditions_by_mark)
    if requested_mode == "differential" and not differential:
        return [(
            "integrated_mode", "error",
            "integrated.mode is 'differential' but no comparison is declared; "
            "declare pairs in comparisons.csv or add a second condition "
            f"(groups found: {names})",
        )]
    if not differential:
        return [(
            "integrated_mode", "warning",
            "no comparison declared; the pipeline will run in descriptive mode "
            "and will not produce diff_* outputs",
        )]
    return []


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--chip-metadata", default=None)
    parser.add_argument("--rna-metadata", default=None)
    parser.add_argument("--count-matrix", default=None)
    parser.add_argument("--mode", default="auto")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    config = {}
    if os.path.exists(args.config):
        import yaml

        with open(args.config) as handle:
            config = yaml.safe_load(handle) or {}

    rna_name = workspace.default_name(config, "RNAseq")
    anno_id = annotation.resolve(config)["id"]
    chip_metadata = args.chip_metadata or workspace.metadata_from_config(config, "ChIPseq")
    rna_metadata = args.rna_metadata or workspace.metadata_from_config(config, "RNAseq")
    count_matrix = args.count_matrix or str(
        workspace.stage_dir("RNAseq", rna_name, "3_RC", annotation=anno_id)
        / "gene_count_matrix.csv"
    )
    chip_name = workspace.default_name(config, "ChIPseq")
    # The comparison file sits next to the metadata, so a --chip-metadata
    # override needs no second path.
    comparisons_path = os.path.join(os.path.dirname(chip_metadata), workspace.COMPARISON_NAME)
    out = args.out or str(
        workspace.provenance_dir(rna_name, chip_name, annotation=anno_id)
        / "validation_report.tsv"
    )

    problems = []
    try:
        chip = read_chip_metadata(chip_metadata)
        chip_inputs = read_chip_inputs(chip_metadata)
    except (OSError, ValueError) as error:
        problems.append(("chip_metadata", "error", f"cannot read metadata: {error}"))
        chip = []
        chip_inputs = None
    problems += check_chip(chip, comparisons_path, chip_inputs)
    problems += check_rna(read_rna(rna_metadata), count_matrix)
    problems += check_mode(chip, args.mode, comparisons_path)

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w") as handle:
        handle.write("check\tlevel\tmessage\n")
        for check, level, message in problems:
            handle.write(f"{check}\t{level}\t{message}\n")

    errors = [p for p in problems if p[1] == "error"]
    for check, level, message in problems:
        print(f"[{level}] {check}: {message}", file=sys.stderr)
    if errors:
        raise SystemExit(f"{len(errors)} design problem(s); see {out}")
    print(f"[validate] {len(problems)} note(s), no errors -> {out}")


if __name__ == "__main__":
    main()
