"""Metadata handling shared by all APIARC modules.

The pipeline contract is the explicit per-dataset schema in
``workflow/data/<omics>/<dataset>/metadata.csv``.

ChIP-seq uses the role-aware schema::

    sample,role,sample_name,condition,input_sample

``role`` is ``IP`` or ``Input``. Every library holds one row: an Input row
names its own group in ``sample_name`` and its repeat in ``condition`` and
leaves ``input_sample`` empty; an IP row names the Input library it belongs to.
``sample_name`` names the group; rows sharing it are replicates, numbered
``rep1``/``rep2``/``rep3`` in ``condition``. Which group is compared to which is
not hard-coded: the preparation page writes ``comparisons.csv`` next to the
metadata and the pipeline runs exactly those pairs.

Files written before ``role`` existed are still parsed: all their rows are IP
rows and the Input is named only through ``input_sample`` /
``input_sample_name``. The same applies to ``input_sample_name``: the role-aware
table resolves the name from the Input row instead of storing it twice.

RNA-seq keeps ``sample,sample_name,group``.

The older ChIP schema (``sample,mark,condition,replicate,input_id,fastq_prefix``
and the even older ``IP sample,Input,IP_name``) is still parsed so that existing
projects keep running; those samples are keyed by ``mark_condition`` and their
comparisons are enumerated from the mark/condition pairs as before.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

SINGLE_CONDITION = "single"

_FIELD_ALIASES = {
    "sample": {"sample", "ip sample", "ip_sample", "srr", "sample_id"},
    "mark": {"mark", "target", "ip_name"},
    "condition": {"condition", "group", "treatment"},
    "replicate": {"replicate", "rep", "bio_rep"},
    "input_id": {"input_id", "input", "control"},
    "fastq_prefix": {"fastq_prefix", "prefix"},
    "sample_name": {"sample_name"},
    "input_sample": {"input_sample"},
    "input_sample_name": {"input_sample_name"},
    "role": {"role", "library_role"},
}

_REP_SUFFIX = re.compile(r"_rep(\d+)$", re.IGNORECASE)
_REP_ONLY = re.compile(r"^(?:rep)?(\d+)$", re.IGNORECASE)


@dataclass(frozen=True)
class ChipSample:
    sample: str
    mark: str
    condition: str
    replicate: int
    input_id: str
    fastq_prefix: str = ""
    sample_name: str = ""
    input_name: str = ""
    role: str = "IP"

    @property
    def group(self) -> str:
        # The five-column schema names the group directly and has no second
        # condition; the legacy schema still keys groups by mark_condition.
        if self.sample_name and self.mark == self.sample_name:
            return self.sample_name
        return group_key(self.mark, self.condition)

    @property
    def comparison_key(self) -> str:
        """Value compared in ``reference``/``treatment``.

        New schema: the group name. Legacy schema: the condition.
        """
        return self.sample_name if self.mark == self.sample_name else self.condition

    @property
    def is_single_condition(self) -> bool:
        return self.condition == SINGLE_CONDITION


def group_key(mark: str, condition: str) -> str:
    """Stable key for a mark/condition pair, used in output paths."""
    return f"{mark}_{condition}"


def uses_group_names(samples: Iterable[ChipSample]) -> bool:
    """True when the five-column schema named every sample after its group.

    In that schema ``mark`` is filled with the ``sample_name`` rather than a
    histone mark, so code that quotes a mark has to say so explicitly.
    """
    items = list(samples)
    return bool(items) and all(s.sample_name and s.mark == s.sample_name for s in items)


def display_mark(name: str) -> str:
    """Mark part of a ``<mark>_<condition>`` group name.

    Group names carry both halves because the five-column schema has no mark
    column; a name without an underscore is its own mark.
    """
    return name.rpartition("_")[0] or name


def comparison_name(mark: str, reference: str, treatment: str) -> str:
    return f"{mark}_{treatment}_vs_{reference}"


def _canonical_columns(header: Iterable[str]) -> dict[str, str]:
    """Map raw header names to canonical field names."""
    mapping: dict[str, str] = {}
    for raw in header:
        key = raw.strip().strip('"').lower()
        for canonical, aliases in _FIELD_ALIASES.items():
            if key in aliases:
                mapping[raw] = canonical
                break
    return mapping


def _parse_chip_rows(
    path: str,
) -> tuple[list[ChipSample], list[ChipSample], bool]:
    """Split one ChIP metadata file into its IP and Input rows.

    Three layouts are understood. The role-aware table lists every library and
    says which ones are Inputs; the five-column table names Inputs only
    through ``input_sample``; the older mark tables predate both. Only the
    first can hold an Input row, so ``role_aware`` tells the caller whether it
    may demand that every cited Input is declared.
    """
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: empty metadata file")
        mapping = _canonical_columns(reader.fieldnames)
        declared = set(mapping.values())
        role_aware = "role" in declared
        if role_aware:
            required = {"sample", "role", "sample_name", "condition", "input_sample"}
        elif "input_sample" in declared:
            required = {"sample", "sample_name", "condition", "input_sample"}
        else:
            required = {"sample", "mark", "input_id"}
        missing = required - declared
        if missing:
            raise ValueError(
                f"{path}: missing required column(s) {sorted(missing)}; "
                f"found {reader.fieldnames}"
            )
        rows = list(reader)

    legacy = (
        not role_aware
        and "input_sample" not in declared
        and "condition" not in declared
    )
    samples: list[ChipSample] = []
    seen: set[str] = set()
    blank_inputs: list[int] = []

    for row in rows:
        record = {mapping[k]: (v or "").strip() for k, v in row.items() if k in mapping}
        sample = record["sample"]
        if not sample:
            continue
        if sample in seen:
            raise ValueError(f"{path}: duplicated sample id '{sample}'")
        seen.add(sample)

        role = "IP"
        if role_aware:
            raw_role = record.get("role", "").lower()
            if raw_role not in ("ip", "input"):
                raise ValueError(
                    f"{path}: sample '{sample}' has role "
                    f"'{record.get('role', '')}'; expected IP or Input"
                )
            role = "IP" if raw_role == "ip" else "Input"

        if role == "Input":
            sample_name = record.get("sample_name", "")
            if not sample_name:
                raise ValueError(f"{path}: Input '{sample}' has no sample_name")
            if record.get("input_sample", ""):
                raise ValueError(
                    f"{path}: Input '{sample}' names an input_sample; only IP rows do"
                )
            replicate_raw = record.get("condition", "")
            if replicate_raw:
                match = _REP_ONLY.match(replicate_raw)
                if match is None:
                    raise ValueError(
                        f"{path}: Input '{sample}' has condition '{replicate_raw}'; "
                        "expected rep1, rep2 or rep3"
                    )
                replicate_raw = match.group(1)
            else:
                blank_inputs.append(len(samples))
            mark = sample_name
            condition = SINGLE_CONDITION
            input_id = ""
            input_name = ""
        elif role_aware or "input_sample" in declared:
            sample_name = record.get("sample_name", "")
            if not sample_name:
                raise ValueError(f"{path}: sample '{sample}' has no sample_name")
            replicate_raw = record.get("condition", "")
            match = _REP_ONLY.match(replicate_raw)
            if match is None:
                raise ValueError(
                    f"{path}: sample '{sample}' has condition '{replicate_raw}'; "
                    "expected rep1, rep2 or rep3"
                )
            replicate_raw = match.group(1)
            # The group name doubles as the mark namespace; this schema has no
            # second condition.
            mark = sample_name
            condition = SINGLE_CONDITION
            input_id = record.get("input_sample", "")
            input_name = record.get("input_sample_name", "") or input_id
        elif legacy:
            mark = record.get("mark", "")
            match = _REP_SUFFIX.search(mark)
            replicate_raw = match.group(1) if match else "1"
            mark = _REP_SUFFIX.sub("", mark)
            condition = SINGLE_CONDITION
            input_id = record["input_id"]
            sample_name = group_key(mark, condition)
            input_name = input_id
        else:
            mark = record.get("mark", "")
            replicate_raw = record.get("replicate", "")
            condition = record.get("condition", "") or SINGLE_CONDITION
            input_id = record["input_id"]
            sample_name = group_key(mark, condition)
            input_name = input_id

        if not mark:
            raise ValueError(f"{path}: sample '{sample}' has no mark")
        try:
            replicate = int(replicate_raw or 1)
        except ValueError as exc:
            raise ValueError(
                f"{path}: sample '{sample}' has non-integer replicate '{replicate_raw}'"
            ) from exc

        samples.append(
            ChipSample(
                sample=sample,
                mark=mark,
                condition=condition,
                replicate=replicate,
                input_id=input_id,
                fastq_prefix=record.get("fastq_prefix", ""),
                sample_name=sample_name,
                input_name=input_name,
                role=role,
            )
        )

    input_names = {s.sample: s.sample_name for s in samples if s.role == "Input"}
    if blank_inputs:
        # A group may leave the repeat out when it holds one library; a second
        # blank row is numbered instead of colliding with the first.
        blanks = set(blank_inputs)
        taken: dict[str, set[int]] = {}
        for index, item in enumerate(samples):
            if index not in blanks and item.role == "Input" and item.replicate:
                taken.setdefault(item.sample_name, set()).add(item.replicate)
        for index in blank_inputs:
            item = samples[index]
            numbers = taken.setdefault(item.sample_name, set())
            number = 1
            while number in numbers:
                number += 1
            numbers.add(number)
            samples[index] = replace(item, replicate=number)
    if role_aware:
        for item in samples:
            if item.role != "IP":
                continue
            if not item.input_id:
                raise ValueError(f"{path}: sample '{item.sample}' has no input_sample")
            if item.input_id not in input_names:
                raise ValueError(
                    f"{path}: sample '{item.sample}' cites input_sample "
                    f"'{item.input_id}', which is not a row with role Input"
                )
    resolved = []
    for item in samples:
        if item.role == "Input":
            resolved.append(item)
            continue
        name = input_names.get(item.input_id, "")
        resolved.append(replace(item, input_name=name or item.input_name or item.input_id))

    ips = [item for item in resolved if item.role == "IP"]
    inputs = [item for item in resolved if item.role == "Input"]
    _assert_unique_replicates(path, ips)
    _assert_unique_inputs(path, inputs)
    return ips, inputs, role_aware


def read_chip_metadata(path: str) -> list[ChipSample]:
    """IP rows of a ChIP metadata file, in file order."""
    return _parse_chip_rows(path)[0]


def read_chip_inputs(path: str) -> list[ChipSample] | None:
    """Input rows when the file carries a role column.

    ``None`` means an older file, which has no place to declare an Input row;
    the design check must not demand one from those projects.
    """
    _ips, inputs, role_aware = _parse_chip_rows(path)
    return inputs if role_aware else None


def read_comparisons(path: str | Path) -> list[tuple[str, str]]:
    """``(reference, treatment)`` pairs declared for this dataset.

    Both sides are group names from the ``sample_name`` column. Rows that are
    half empty or name the same group twice are dropped: the page refuses them
    too, and silently running them would produce a meaningless contrast.
    """
    path = Path(path)
    if not path.is_file():
        return []
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    pairs: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        reference = (row.get("reference") or "").strip()
        treatment = (row.get("treatment") or "").strip()
        if not reference or not treatment or reference == treatment:
            continue
        if (reference, treatment) in seen:
            continue
        seen.add((reference, treatment))
        pairs.append((reference, treatment))
    return pairs


def _assert_unique_replicates(path: str, samples: list[ChipSample]) -> None:
    seen: dict[tuple[str, str, int], str] = {}
    for s in samples:
        key = (s.mark, s.condition, s.replicate)
        if key in seen:
            raise ValueError(
                f"{path}: replicate {s.replicate} of {s.mark}/{s.condition} "
                f"is defined twice ('{seen[key]}' and '{s.sample}')"
            )
        seen[key] = s.sample


def _assert_unique_inputs(path: str, inputs: list[ChipSample]) -> None:
    """One library per Input group and replicate number."""
    seen: dict[tuple[str, int], str] = {}
    for s in inputs:
        key = (s.sample_name, s.replicate)
        if key in seen:
            raise ValueError(
                f"{path}: Input {s.sample_name} replicate {s.replicate} "
                f"is defined twice ('{seen[key]}' and '{s.sample}')"
            )
        seen[key] = s.sample


def mark_conditions(samples: Iterable[ChipSample]) -> dict[str, list[str]]:
    """mark -> sorted list of conditions."""
    table: dict[str, set[str]] = {}
    for s in samples:
        table.setdefault(s.mark, set()).add(s.condition)
    return {mark: sorted(conds) for mark, conds in table.items()}


def groups(samples: Iterable[ChipSample]) -> list[str]:
    """All mark/condition group keys, used as Snakemake wildcards."""
    return sorted({s.group for s in samples})


def samples_by_group(samples: Iterable[ChipSample]) -> dict[str, list[ChipSample]]:
    table: dict[str, list[ChipSample]] = {}
    for s in samples:
        table.setdefault(s.group, []).append(s)
    return {g: sorted(items, key=lambda x: x.replicate) for g, items in table.items()}


def input_ids(samples: Iterable[ChipSample]) -> list[str]:
    return sorted({s.input_id for s in samples})


def resolve_mode(
    config: dict,
    samples: list[ChipSample],
    comparisons_path: str | Path | None = None,
) -> str:
    """Return 'differential' or 'descriptive'.

    ``integrated.mode`` may force a mode. Otherwise the five-column schema is
    differential as soon as the dataset declares comparison pairs; the legacy
    schema is differential when one mark carries two conditions with enough
    replicates.
    """
    requested = str(config.get("integrated", {}).get("mode", "auto")).lower()
    if requested in {"differential", "descriptive"}:
        return requested

    if uses_group_names(samples):
        if comparisons_path and read_comparisons(comparisons_path):
            return "differential"
        return "descriptive"

    min_reps = int(config.get("chip", {}).get("diff", {}).get("min_replicates", 2))
    by_group = samples_by_group(samples)
    for mark, conditions in mark_conditions(samples).items():
        if len(conditions) < 2:
            continue
        if all(
            len(by_group[group_key(mark, cond)]) >= min_reps for cond in conditions
        ):
            return "differential"
    return "descriptive"


def reference_condition(config: dict, mark: str, conditions: list[str]) -> str:
    """Condition used as the reference for a mark."""
    configured = config.get("chip", {}).get("reference_condition", "auto")
    if configured and str(configured).lower() != "auto":
        if configured not in conditions:
            raise ValueError(
                f"chip.reference_condition='{configured}' is not a condition of "
                f"{mark} (available: {', '.join(conditions)})"
            )
        return str(configured)
    if len(conditions) < 2:
        return conditions[0]
    return sorted(conditions)[0]


def comparisons(
    samples: list[ChipSample],
    config: dict,
    comparisons_path: str | Path | None = None,
) -> list[tuple[str, str, str]]:
    """Return (mark, reference, treatment) for every valid comparison.

    The five-column schema reads the pairs from ``comparisons.csv`` and uses the
    treatment group as the path namespace, so one treatment can be contrasted
    against more than one reference. The legacy schema still enumerates every
    condition against the mark's reference condition.
    """
    min_reps = int(config.get("chip", {}).get("diff", {}).get("min_replicates", 2))
    by_group = samples_by_group(samples)
    if uses_group_names(samples):
        out: list[tuple[str, str, str]] = []
        pairs = read_comparisons(comparisons_path) if comparisons_path else []
        for reference, treatment in pairs:
            if len(by_group.get(reference, [])) < min_reps:
                continue
            if len(by_group.get(treatment, [])) < min_reps:
                continue
            out.append((treatment, reference, treatment))
        return out

    out = []
    for mark, conditions in mark_conditions(samples).items():
        if len(conditions) < 2:
            continue
        ref = reference_condition(config, mark, conditions)
        if len(by_group[group_key(mark, ref)]) < min_reps:
            continue
        for cond in conditions:
            if cond == ref:
                continue
            if len(by_group[group_key(mark, cond)]) < min_reps:
                continue
            out.append((mark, ref, cond))
    return out


def comparison_samples(
    samples: Iterable[ChipSample], reference: str, treatment: str, mark: str = ""
) -> list[ChipSample]:
    """Samples entering one contrast, in a stable order.

    Legacy schema: the mark's libraries whose condition is one of the pair.
    Five-column schema: the two groups; ``mark`` is only a path namespace there.
    """
    items = list(samples)
    named = uses_group_names(items)
    picked = [
        s
        for s in items
        if s.comparison_key in {reference, treatment}
        and (named or not mark or s.mark == mark)
    ]
    # Reference first: count matrices then read left to right in the same order
    # as the R scripts' factor levels.
    rank = {reference: 0, treatment: 1}
    return sorted(picked, key=lambda s: (rank.get(s.comparison_key, 2), s.replicate))


def read_rna_metadata(path: str) -> list[dict]:
    """Read RNA-seq metadata (sample, sample_name, group: T=reference, P=treatment)."""
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    cleaned = []
    for row in rows:
        row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
        if row.get("sample"):
            cleaned.append(row)
    if not cleaned:
        raise ValueError(f"{path}: no RNA-seq samples found")
    return cleaned
