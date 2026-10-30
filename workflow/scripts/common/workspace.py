"""Dataset layout contract shared by the pipeline and the dashboard.

A dataset is a self-contained directory::

    workflow/data/<omics>/<dataset>/
        metadata.csv        library table, schema depends on the omics
        dataset.yaml        descriptor: species, assembly, provenance
        *.fastq.gz          raw reads, single-end or paired-end

Results mirror the input layout, one directory per annotation so that two
reference versions never overwrite each other::

    result/<omics>/<dataset>/<assembly>-<source>/<stage>/
    result/Integrated/<rna>__<chip>/<assembly>-<source>/<stage>/
    result/provenance/<rna>__<chip>/<assembly>-<source>/

The annotation segment is optional in the helpers below: an empty value keeps
the flat layout used before annotations became project-level settings.

Both sides import this module so that a path is spelled in exactly one place.
"""

from __future__ import annotations

import csv
import datetime as _dt
import os
import re
from pathlib import Path

OMICS = ("RNAseq", "ChIPseq")

OMICS_LABEL = {"RNAseq": "RNA-seq", "ChIPseq": "ChIP-seq"}

CONFIG_SECTIONS = {"RNAseq": ("RNAseq",), "ChIPseq": ("chip", "ChIPseq")}

DATA_ROOT = Path("workflow") / "data"
RESULT_ROOT = Path("result")

META_NAME = "metadata.csv"
COMPARISON_NAME = "comparisons.csv"
DESCRIPTOR_NAME = "dataset.yaml"
MANIFEST_NAME = "rawdata_manifest.csv"
FINISH_NAME = "finished.txt"

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

CHIP_COLUMNS = (
    "sample",
    "role",
    "sample_name",
    "condition",
    "input_sample",
)
RNA_COLUMNS = ("sample", "sample_name", "group")
COMPARISON_COLUMNS = ("reference", "treatment")

# RNA groups are a closed set: ``control`` is the DESeq2 reference and each
# treat is contrasted against it once. ChIP replicates name the repeat within
# one sample_name group.
RNA_GROUPS = ("control", "treatA", "treatB", "treatC", "treatD")
CHIP_REPLICATES = ("rep1", "rep2", "rep3")

STAGE_DIRS = {
    "RNAseq": {
        "1_prepare": "1_rawdata",
        "2_QC": "2_cleandata",
        "3_RC": "3_align_quant",
        "4_DEseq": "4_deseq",
    },
    "ChIPseq": {
        "1_prepare": "1_rawdata",
        "2_QC": "2_cleandata",
        "3_CC": "3_align",
        "4_peak_result": "4_peaks",
        "5_diffbind": "5_diffbind",
    },
}

INTEGRATED_STAGES = {
    "0_stats": "0_stats",
    "1_enrichment": "1_enrichment",
    "2_motif": "2_motif",
    "3_enhancer": "3_enhancer",
    "4_common": "4_common",
}

_SAFE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

_READ_TOKENS = ("_R1", "_R2", "_1", "_2")


class WorkspaceError(ValueError):
    """Raised when a dataset name or layout violates the contract."""


def default_name(config: dict, omics: str) -> str:
    """Dataset name declared for an omics, falling back to the experiment id."""
    for key in CONFIG_SECTIONS.get(omics, (omics,)):
        section = config.get(key) or {}
        name = str(section.get("dataset") or "").strip()
        if name:
            return name
    fallback = str(config.get("experiment") or "").strip()
    return fallback or "dataset"


def dataset_dir(omics: str, name: str, root: Path | str = DATA_ROOT) -> Path:
    _check_omics(omics)
    _check_name(name)
    return Path(root) / omics / name


def result_dir(
    omics: str, name: str, root: Path | str = RESULT_ROOT, annotation: str = ""
) -> Path:
    _check_omics(omics)
    _check_name(name)
    return annotate_dir(Path(root) / omics / name, annotation)


def integrated_id(rna_name: str, chip_name: str) -> str:
    _check_name(rna_name)
    _check_name(chip_name)
    return f"{rna_name}__{chip_name}"


def annotate_dir(base: Path, annotation: str = "") -> Path:
    """Insert the ``<assembly>-<source>`` segment that keeps versions apart."""
    if not annotation:
        return base
    _check_name(annotation)
    return base / annotation


def integrated_dir(
    rna_name: str, chip_name: str, root: Path | str = RESULT_ROOT, annotation: str = ""
) -> Path:
    return annotate_dir(Path(root) / "Integrated" / integrated_id(rna_name, chip_name), annotation)


def provenance_dir(
    rna_name: str, chip_name: str, root: Path | str = RESULT_ROOT, annotation: str = ""
) -> Path:
    return annotate_dir(Path(root) / "provenance" / integrated_id(rna_name, chip_name), annotation)


def stage_dir(
    omics: str, name: str, stage: str, root: Path | str = RESULT_ROOT, annotation: str = ""
) -> Path:
    return result_dir(omics, name, root, annotation) / STAGE_DIRS[omics][stage]


def integrated_stage_dir(
    rna_name: str,
    chip_name: str,
    stage: str,
    root: Path | str = RESULT_ROOT,
    annotation: str = "",
) -> Path:
    return integrated_dir(rna_name, chip_name, root, annotation) / INTEGRATED_STAGES[stage]


def metadata_path(omics: str, name: str, root: Path | str = DATA_ROOT) -> Path:
    return dataset_dir(omics, name, root) / META_NAME


def comparisons_path(omics: str, name: str, root: Path | str = DATA_ROOT) -> Path:
    return dataset_dir(omics, name, root) / COMPARISON_NAME


def metadata_from_config(config: dict, omics: str) -> str:
    """Metadata path for a configured dataset, honouring an explicit override."""
    for key in CONFIG_SECTIONS.get(omics, (omics,)):
        section = config.get(key) or {}
        override = str(section.get("metadata") or "").strip()
        if override:
            return override
    return str(metadata_path(omics, default_name(config, omics)))


def descriptor_path(omics: str, name: str, root: Path | str = DATA_ROOT) -> Path:
    return dataset_dir(omics, name, root) / DESCRIPTOR_NAME


def list_datasets(omics: str, root: Path | str = DATA_ROOT) -> list[str]:
    base = Path(root) / omics
    if not base.is_dir():
        return []
    names = []
    for entry in sorted(base.iterdir()):
        if not entry.is_dir() or entry.name.startswith("."):
            continue
        if (entry / META_NAME).is_file() or any(entry.glob("*.fastq*")):
            names.append(entry.name)
    return names


def safe_name(name: str) -> bool:
    return bool(_SAFE_NAME.match(str(name or "").strip()))


def read_descriptor(omics: str, name: str, root: Path | str = DATA_ROOT) -> dict:
    path = descriptor_path(omics, name, root)
    if not path.is_file():
        return {}
    import yaml

    with open(path) as handle:
        data = yaml.safe_load(handle) or {}
    return data if isinstance(data, dict) else {}


def write_descriptor(omics: str, name: str, extra: dict, root: Path | str = DATA_ROOT) -> Path:
    import yaml

    path = descriptor_path(omics, name, root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "omics": omics,
        "dataset": name,
        "updated_at": _now(),
    }
    payload.update({k: v for k, v in (extra or {}).items() if v not in (None, "")})
    with open(path, "w") as handle:
        yaml.safe_dump(payload, handle, sort_keys=True, allow_unicode=True)
    return path


def read_metadata_rows(path: Path | str) -> list[dict]:
    path = Path(path)
    if not path.is_file():
        return []
    with open(path, newline="") as handle:
        return [row for row in csv.DictReader(handle)]


def write_metadata_rows(
    path: Path | str, rows: list[dict], columns: tuple[str, ...] | None = None
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        # ``sample_name`` is common to both schemas; ``input_sample`` is not.
        columns = (
            CHIP_COLUMNS if "input_sample" in (rows[0] if rows else {}) else RNA_COLUMNS
        )
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ("" if row.get(key) is None else row.get(key)) for key in columns})
    return path


def columns_for(omics: str) -> tuple[str, ...]:
    return CHIP_COLUMNS if omics == "ChIPseq" else RNA_COLUMNS


def row_role(row: dict) -> str:
    """``IP`` or ``Input`` for one metadata row.

    Files written before the role column existed only held IP rows, so a
    missing role is an IP row rather than an error.
    """
    return "Input" if (row.get("role") or "").strip().lower() == "input" else "IP"


def declared_inputs(omics: str, rows: list[dict]) -> list[str]:
    """Input library ids declared by the ChIP table, in row order.

    A role-aware file gives each Input its own row; older five-column files
    only named the library in ``input_sample``. Stage 1 links both, otherwise
    peak calling opens a BAM that was never linked. The old schema spelled the
    column ``input_id``.
    """
    if omics != "ChIPseq":
        return []
    seen: dict[str, None] = {}
    for row in rows:
        if row_role(row) == "Input":
            sample = (row.get("sample") or "").strip()
            if sample:
                seen.setdefault(sample, None)
        value = (row.get("input_sample") or row.get("input_id") or "").strip()
        if value:
            seen.setdefault(value, None)
    return list(seen)


def index_fastq(directory: Path | str) -> dict[str, dict[str, Path]]:
    """Map a read prefix to its files: ``{prefix: {"r1": path, "r2": path}}``."""
    directory = Path(directory)
    table: dict[str, dict[str, Path]] = {}
    if not directory.is_dir():
        return table
    for entry in sorted(directory.iterdir()):
        if not entry.is_file() and not entry.is_symlink():
            continue
        name = entry.name
        suffix = next((s for s in FASTQ_SUFFIXES if name.endswith(s)), None)
        if suffix is None:
            continue
        stem = name[: -len(suffix)]
        prefix, read = _split_read(stem)
        slot = table.setdefault(prefix, {})
        if read == 2:
            slot["r2"] = entry
        else:
            slot["r1"] = entry
    return table


def _split_read(stem: str) -> tuple[str, int]:
    for token in _READ_TOKENS:
        if stem.endswith(token):
            return stem[: -len(token)], int(token[-1])
    return stem, 1


def lookup_fastq(table: dict[str, dict[str, Path]], key: str) -> dict[str, Path]:
    if key in table:
        return table[key]
    lowered = {name.lower(): name for name in table}
    match = lowered.get(str(key).lower())
    return table.get(match, {}) if match else {}


def dataset_summary(omics: str, name: str, root: Path | str = DATA_ROOT) -> dict:
    directory = dataset_dir(omics, name, root)
    table = index_fastq(directory)
    rows = read_metadata_rows(metadata_path(omics, name, root))
    ip_rows = [row for row in rows if row_role(row) == "IP"]
    input_rows = [row for row in rows if row_role(row) == "Input"]
    paired = sum(1 for slot in table.values() if slot.get("r2"))
    single = sum(1 for slot in table.values() if not slot.get("r2"))
    missing = []
    seen_missing: set[str] = set()
    for row in rows:
        key = (row.get("fastq_prefix") or row.get("sample") or "").strip()
        if key and key not in table and key not in seen_missing:
            missing.append(key)
            seen_missing.add(key)
        # Older files name the Input only through the IP row; a role-aware
        # table reports the missing library once, on its own row.
        ref = (row.get("input_sample") or row.get("input_id") or "").strip()
        if omics == "ChIPseq" and ref and ref not in table and ref not in seen_missing:
            missing.append(ref)
            seen_missing.add(ref)
    return {
        "omics": omics,
        "dataset": name,
        "path": str(directory),
        "exists": directory.is_dir(),
        "samples": len(ip_rows),
        "inputs": len(input_rows),
        "fastq_pairs": paired,
        "fastq_singles": single,
        "missing_fastq": missing,
        "descriptor": read_descriptor(omics, name, root),
    }


def list_well_formed(omics: str, root: Path | str = DATA_ROOT) -> list[dict]:
    return [dataset_summary(omics, name, root) for name in list_datasets(omics, root)]


def _check_omics(omics: str) -> None:
    if omics not in OMICS:
        raise WorkspaceError(f"unknown omics '{omics}'; expected one of {list(OMICS)}")


def _check_name(name: str) -> None:
    if not safe_name(name):
        raise WorkspaceError(
            f"invalid dataset name '{name}'; use letters, digits, dot, dash or underscore"
        )


def _now() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def ensure_layout(root: Path | str = DATA_ROOT) -> None:
    for omics in OMICS:
        (Path(root) / omics).mkdir(parents=True, exist_ok=True)


def rel(path: Path | str) -> str:
    path = Path(path)
    try:
        return str(path.resolve().relative_to(repo_root()))
    except (ValueError, OSError):
        return str(path)


def human_bytes(value: float) -> str:
    """Size for a caption, in the same decimal units as GEO's ``size_MB``."""
    number = float(value or 0)
    for unit in ("B", "KB", "MB", "GB"):
        if number < 1000:
            return f"{number:.1f} {unit}"
        number /= 1000
    return f"{number:.1f} TB"


def is_empty_dir(path: Path | str) -> bool:
    path = Path(path)
    if not path.is_dir():
        return True
    return not any(entry for entry in path.iterdir())


def touch(path: Path | str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")


def marker(
    omics: str, name: str, stage: str, root: Path | str = RESULT_ROOT, annotation: str = ""
) -> Path:
    return stage_dir(omics, name, stage, root, annotation) / FINISH_NAME


def integrated_marker(
    rna_name: str,
    chip_name: str,
    stage: str,
    root: Path | str = RESULT_ROOT,
    annotation: str = "",
) -> Path:
    return integrated_stage_dir(rna_name, chip_name, stage, root, annotation) / FINISH_NAME


def env_manifest(
    omics: str, name: str, root: Path | str = RESULT_ROOT, annotation: str = ""
) -> Path:
    return stage_dir(omics, name, "1_prepare", root, annotation) / MANIFEST_NAME
