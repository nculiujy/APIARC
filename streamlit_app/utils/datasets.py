from __future__ import annotations

import os
import sys
from pathlib import Path

from . import env_paths, i18n

if str(env_paths.SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(env_paths.SCRIPT_ROOT))

from common import annotation as anno_registry  # noqa: E402
from common import geo as geo_module  # noqa: E402
from common import workspace as ws  # noqa: E402

OMICS = ws.OMICS
OMICS_LABEL = ws.OMICS_LABEL
STAGE_DIRS = ws.STAGE_DIRS
INTEGRATED_STAGES = ws.INTEGRATED_STAGES
CONFIG_SECTION = {"RNAseq": "RNAseq", "ChIPseq": "chip"}
RNA_GROUPS = ws.RNA_GROUPS
CHIP_REPLICATES = ws.CHIP_REPLICATES
human_bytes = ws.human_bytes

DATA_ROOT = env_paths.DATA_ROOT
RESULT_ROOT = env_paths.RESULT_ROOT

FINISH = ws.FINISH_NAME
META = ws.META_NAME
DESCRIPTOR = ws.DESCRIPTOR_NAME

_CONFIG_CACHE: dict[str, tuple[float, dict]] = {}


def columns(omics: str) -> tuple[str, ...]:
    return ws.columns_for(omics)


def current_config() -> dict:
    """Config for the annotation segment; cached on mtime to keep path lookups cheap."""
    from . import runtime

    path = env_paths.CONFIG_FILE
    try:
        stamp = path.stat().st_mtime
    except OSError:
        stamp = 0.0
    cached = _CONFIG_CACHE.get("value")
    if cached and cached[0] == stamp:
        return cached[1]
    config = runtime.load_config(path)
    _CONFIG_CACHE["value"] = (stamp, config)
    return config


def annotation_id(config: dict | None = None) -> str:
    """``<assembly>-<source>`` segment shared by both omics and the joint run."""
    return anno_registry.resolve(config if config is not None else current_config())["id"]


def annotation_info(config: dict | None = None) -> dict:
    return anno_registry.resolve(config if config is not None else current_config())


def _version(config: dict | None, annotation: str | None) -> str:
    """Path segment: the caller's version when it names one, else the config's."""
    return annotation if annotation is not None else annotation_id(config)


def result_versions(omics: str, name: str) -> list[str]:
    """Annotations already present under one dataset; shows what a switch left behind."""
    base = RESULT_ROOT / omics / name
    if not base.is_dir():
        return []
    return sorted(
        entry.name for entry in base.iterdir() if entry.is_dir() and not entry.name.startswith(".")
    )


def result_inventory(omics: str) -> list[tuple[str, str]]:
    """Every ``(dataset, annotation)`` directory that exists under one omics."""
    base = RESULT_ROOT / omics
    if not base.is_dir():
        return []
    pairs: list[tuple[str, str]] = []
    for dataset in sorted(base.iterdir()):
        if not dataset.is_dir() or dataset.name.startswith("."):
            continue
        for version in sorted(dataset.iterdir()):
            if version.is_dir() and not version.name.startswith("."):
                pairs.append((dataset.name, version.name))
    return pairs


def config_section(config: dict, omics: str) -> str:
    for key in ws.CONFIG_SECTIONS.get(omics, (omics,)):
        section = config.get(key) or {}
        if str(section.get("dataset") or "").strip():
            return key
    return CONFIG_SECTION[omics]


def config_dataset(config: dict, omics: str) -> str:
    return ws.default_name(config, omics)


def joint_id(config: dict) -> str:
    return ws.integrated_id(
        ws.default_name(config, "RNAseq"), ws.default_name(config, "ChIPseq")
    )


def directory(omics: str, name: str) -> Path:
    return ws.dataset_dir(omics, name, DATA_ROOT)


def result_directory(
    omics: str, name: str, config: dict | None = None, annotation: str | None = None
) -> Path:
    return ws.result_dir(omics, name, RESULT_ROOT, _version(config, annotation))


def joint_directory(
    config: dict,
    rna: str | None = None,
    chip: str | None = None,
    annotation: str | None = None,
) -> Path:
    return ws.integrated_dir(
        rna or ws.default_name(config, "RNAseq"),
        chip or ws.default_name(config, "ChIPseq"),
        RESULT_ROOT,
        _version(config, annotation),
    )


def provenance_directory(config: dict) -> Path:
    return ws.provenance_dir(
        ws.default_name(config, "RNAseq"),
        ws.default_name(config, "ChIPseq"),
        RESULT_ROOT,
        annotation_id(config),
    )


def stage(
    omics: str,
    name: str,
    key: str,
    config: dict | None = None,
    annotation: str | None = None,
) -> Path:
    return ws.stage_dir(omics, name, key, RESULT_ROOT, _version(config, annotation))


def stage_marker(
    omics: str,
    name: str,
    key: str,
    config: dict | None = None,
    annotation: str | None = None,
) -> Path:
    return ws.marker(omics, name, key, RESULT_ROOT, _version(config, annotation))


def joint_stage(
    config: dict,
    key: str,
    rna: str | None = None,
    chip: str | None = None,
    annotation: str | None = None,
) -> Path:
    return ws.integrated_stage_dir(
        rna or ws.default_name(config, "RNAseq"),
        chip or ws.default_name(config, "ChIPseq"),
        key,
        RESULT_ROOT,
        _version(config, annotation),
    )


def joint_marker(
    config: dict,
    key: str,
    rna: str | None = None,
    chip: str | None = None,
    annotation: str | None = None,
) -> Path:
    return ws.integrated_marker(
        rna or ws.default_name(config, "RNAseq"),
        chip or ws.default_name(config, "ChIPseq"),
        key,
        RESULT_ROOT,
        _version(config, annotation),
    )


def available(omics: str) -> list[str]:
    return ws.list_datasets(omics, DATA_ROOT)


def summary(omics: str, name: str) -> dict:
    return ws.dataset_summary(omics, name, DATA_ROOT)


def descriptors(omics: str, name: str) -> dict:
    return ws.read_descriptor(omics, name, DATA_ROOT)


def save_descriptor(omics: str, name: str, extra: dict) -> Path:
    return ws.write_descriptor(omics, name, extra, DATA_ROOT)


def metadata_rows(omics: str, name: str) -> list[dict]:
    return ws.read_metadata_rows(ws.metadata_path(omics, name, DATA_ROOT))


def save_metadata(omics: str, name: str, rows: list[dict]) -> Path:
    path = ws.metadata_path(omics, name, DATA_ROOT)
    return ws.write_metadata_rows(path, rows, columns(omics))


def validate_name(name: str) -> str:
    text = str(name or "").strip()
    if not ws.safe_name(text):
        raise ValueError(i18n.t("ds.error.name_charset"))
    return text


def validate_metadata(omics: str, rows: list[dict]) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    if not rows:
        return [i18n.t("ds.error.empty_metadata")], warnings

    samples = [(row.get("sample") or "").strip() for row in rows]
    if any(not sample for sample in samples):
        errors.append(i18n.t("ds.error.sample_required"))
    duplicates = {sample for sample in samples if samples.count(sample) > 1}
    if duplicates:
        errors.append(i18n.t("ds.error.duplicate_sample", names=", ".join(sorted(duplicates))))

    if omics == "ChIPseq":
        allowed = set(ws.CHIP_REPLICATES)
        role_aware = any("role" in row for row in rows)
        declared = {
            (row.get("sample") or "").strip()
            for row in rows
            if (row.get("role") or "").strip().lower() == "input"
        }
        seen: dict[tuple[str, str, str], int] = {}
        inputs_without_rep: dict[str, int] = {}
        cited: set[str] = set()
        for row in rows:
            label = row.get("sample", "?")
            sample = (row.get("sample") or "").strip()
            raw_role = (row.get("role") or "").strip().lower()
            if role_aware and raw_role not in ("ip", "input"):
                errors.append(i18n.t("ds.error.role", label=label, value=raw_role))
            role = "input" if raw_role == "input" else "ip"
            sample_name = (row.get("sample_name") or "").strip()
            condition = (row.get("condition") or "").strip()
            ref = (row.get("input_sample") or "").strip()
            if not sample_name:
                errors.append(i18n.t("ds.error.sample_name", label=label))
            if role == "input":
                # An Input library owns a row; it never points at another one.
                if ref:
                    errors.append(i18n.t("ds.error.input_ref", label=label))
                if condition and condition not in allowed:
                    errors.append(
                        i18n.t("ds.error.condition", label=label, value=condition)
                    )
                if condition in allowed:
                    key = ("input", sample_name, condition)
                    seen[key] = seen.get(key, 0) + 1
                elif sample_name:
                    inputs_without_rep[sample_name] = (
                        inputs_without_rep.get(sample_name, 0) + 1
                    )
                continue
            if condition not in allowed:
                errors.append(i18n.t("ds.error.condition", label=label, value=condition))
            if not ref:
                errors.append(i18n.t("ds.error.input_required", label=label))
            elif ref == sample:
                errors.append(i18n.t("ds.error.input_self", label=label))
            elif role_aware and ref not in declared:
                errors.append(i18n.t("ds.error.input_unknown", label=label, ref=ref))
            if ref:
                cited.add(ref)
            key = ("ip", sample_name, condition)
            seen[key] = seen.get(key, 0) + 1
        repeated = [
            f"{group}/{condition}"
            for (_role, group, condition), count in seen.items()
            if count > 1 and group and condition in allowed
        ]
        if repeated:
            errors.append(i18n.t("ds.error.duplicate_pair", names=", ".join(repeated)))
        for group, count in sorted(inputs_without_rep.items()):
            if count > 1:
                warnings.append(
                    i18n.t("ds.warning.input_no_condition", name=group, n=count)
                )
        by_group: dict[str, int] = {}
        for row in rows:
            if (row.get("role") or "").strip().lower() == "input":
                continue
            group = (row.get("sample_name") or "").strip()
            by_group[group] = by_group.get(group, 0) + 1
        for group, count in sorted(by_group.items()):
            if group and count < 2:
                warnings.append(
                    i18n.t("ds.warning.replicates", name=group, n=count)
                )
        if role_aware:
            for ref in sorted(declared - cited):
                warnings.append(i18n.t("ds.warning.input_unused", ref=ref))
    else:
        allowed = set(ws.RNA_GROUPS)
        groups: dict[str, int] = {}
        for row in rows:
            label = row.get("sample", "?")
            group = (row.get("group") or "").strip()
            if not group:
                errors.append(i18n.t("ds.error.group_required", label=label))
            elif group not in allowed:
                errors.append(
                    i18n.t(
                        "ds.error.group_value",
                        label=label,
                        value=group,
                        options=", ".join(sorted(allowed)),
                    )
                )
            if not (row.get("sample_name") or "").strip():
                errors.append(i18n.t("ds.error.sample_name", label=label))
            groups[group] = groups.get(group, 0) + 1
        if "control" not in groups:
            errors.append(i18n.t("ds.error.control_missing"))
        if not set(groups) - {"control", ""}:
            errors.append(i18n.t("ds.error.treat_missing"))
        for group, count in sorted(groups.items()):
            if group and count < 2:
                errors.append(
                    i18n.t("ds.error.group_replicates", group=group, n=count)
                )

    return errors, warnings


def fastq_table(omics: str, name: str) -> dict[str, dict[str, Path]]:
    return ws.index_fastq(directory(omics, name))


def missing_fastq(omics: str, name: str) -> list[str]:
    """Libraries named by the table with no FASTQ in the dataset directory.

    A role-aware table gives the Input its own row, so the IP rows citing it
    would report the same missing library again; each id is listed once.
    """
    table = fastq_table(omics, name)
    rows = metadata_rows(omics, name)
    declared = {(row.get("sample") or "").strip() for row in rows}
    missing = []
    for row in rows:
        key = (row.get("fastq_prefix") or row.get("sample") or "").strip()
        if key and ws.lookup_fastq(table, key) == {} and key not in missing:
            missing.append(key)
        ref = (row.get("input_sample") or row.get("input_id") or "").strip()
        if (
            omics == "ChIPseq"
            and ref
            and ref not in declared
            and ref not in missing
            and ws.lookup_fastq(table, ref) == {}
        ):
            missing.append(ref)
    return missing


def comparison_rows(omics: str, name: str) -> list[dict]:
    return ws.read_metadata_rows(ws.comparisons_path(omics, name, DATA_ROOT))


def save_comparisons(omics: str, name: str, rows: list[dict]) -> Path:
    return ws.write_metadata_rows(
        ws.comparisons_path(omics, name, DATA_ROOT), rows, ws.COMPARISON_COLUMNS
    )


def geo_runs(gse: str, omics: str) -> list[dict]:
    return geo_module.lookup(gse, omics)


def geo_run_info_path(omics: str, name: str) -> Path:
    return directory(omics, name) / "geo" / "SraRunInfo.csv"


def read_run_info(omics: str, name: str) -> list[dict]:
    return geo_module.read_run_info(geo_run_info_path(omics, name))


def write_run_info(omics: str, name: str, rows: list[dict]) -> Path:
    return geo_module.write_run_info(geo_run_info_path(omics, name), rows)


def write_srr_list(omics: str, name: str, runs: list[str]) -> Path:
    path = directory(omics, name) / "geo" / "SRR.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(runs) + ("\n" if runs else ""))
    return path


def srr_list(omics: str, name: str) -> list[str]:
    """Runs picked at step 1, in the order the download will fetch them."""
    path = directory(omics, name) / "geo" / "SRR.txt"
    if not path.is_file():
        return []
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return []
    return [line.strip() for line in text.splitlines() if line.strip()]


def describe(omics: str, name: str) -> list[tuple[str, str]]:
    info = summary(omics, name)
    rows = [
        ("dataset", name),
        ("path", os.path.relpath(info["path"], env_paths.BASE_DIR)),
        ("samples", str(info["samples"])),
        ("paired FASTQ", str(info["fastq_pairs"])),
        ("single FASTQ", str(info["fastq_singles"])),
        ("missing FASTQ", ", ".join(info["missing_fastq"]) or "-"),
    ]
    descriptor = info.get("descriptor") or {}
    for key in ("species", "assembly", "source", "gse", "note"):
        if descriptor.get(key):
            rows.append((key, str(descriptor[key])))
    return rows
