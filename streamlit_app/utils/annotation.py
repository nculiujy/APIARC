from __future__ import annotations

import json
import os
import sys
from pathlib import Path

from . import env_paths, i18n

if str(env_paths.SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(env_paths.SCRIPT_ROOT))

from common import annotation as registry  # noqa: E402

KIND_LABEL = {"rna": "RNA-seq", "chip": "ChIP-seq"}
KIND_ALIGNER = {"rna": "hisat2-build", "chip": "bowtie2-build"}

# Step keys written by scripts/build-annotation.sh. The script owns the order
# and the weights (it is the only side that knows how long each stage takes);
# the labels live here with the rest of the page's wording.
STEP_KEY = {
    "prepare": "annotation.step.prepare",
    "genome": "annotation.step.genome",
    "gtf": "annotation.step.gtf",
    "contigs": "annotation.step.contigs",
    "gene_bed": "annotation.step.gene_bed",
    "splice": "annotation.step.splice",
    "link": "annotation.step.link",
    "index": "annotation.step.index",
}

RNA_INDEX_PARTS = tuple(f"genome.{index}.ht2" for index in range(1, 9))
CHIP_INDEX_PARTS = (
    "1.bt2",
    "2.bt2",
    "3.bt2",
    "4.bt2",
    "rev.1.bt2",
    "rev.2.bt2",
)


def anno_root() -> Path:
    root = Path(registry.ANNO_ROOT)
    return root if root.is_absolute() else env_paths.BASE_DIR / root


def assembly_names() -> list[str]:
    return registry.assembly_names()


def source_names(assembly: str) -> list[str]:
    return registry.sources(assembly)


def assembly_of(species: str) -> str:
    return registry.assembly_of_species(species)


def annotation_id(assembly: str, source: str = "") -> str:
    info = registry.single(assembly, source)
    return info["id"]


def meta(assembly: str, source: str = "") -> dict:
    return registry.single(assembly, source)


def resolve(config: dict) -> dict:
    """Annotation declared by the config, with paths anchored at the repo root."""
    info = registry.resolve(config)
    for key in ("genome", "genome_rna", "gtf", "gene_bed", "rna_index", "chip_index"):
        path = Path(info[key])
        info[key] = path if path.is_absolute() else env_paths.BASE_DIR / path
    for key in ("rna_dir", "chip_dir"):
        path = Path(info[key])
        info[key] = path if path.is_absolute() else env_paths.BASE_DIR / path
    return info


def kind_dir(kind: str, assembly: str) -> Path:
    return anno_root() / registry.KIND_DIRS[kind] / registry.canonical_assembly(assembly)


def index_prefix(kind: str, assembly: str) -> Path:
    """Path stem of the aligner index inside the annotation directory.

    The registry keeps the RNA index down to the ``genome`` stem
    (``.../grcm38/genome``) and the ChIP index only to the assembly-named
    directory (``.../grcm38``), so the two kinds need different arithmetic.
    """
    info = registry.single(assembly)
    folder = kind_dir(kind, assembly)
    path = Path(info["rna_index"] if kind == "rna" else info["chip_index"])
    return folder / path.parent.name / path.name if kind == "rna" else folder / path.name


def expected_files(kind: str, assembly: str, source: str = "") -> list[tuple[str, Path]]:
    info = registry.single(assembly, source)
    folder = kind_dir(kind, assembly)
    genome = Path(info["genome"]).name
    gtf = Path(info["gtf"]).name
    if kind == "rna":
        # Splice sites and exons are named after the index directory
        # (``grcm38``), while the index files carry the ``genome`` stem inside
        # it. Reading one name for both made a finished build look missing.
        index = Path(info["rna_index"])
        items = [
            (i18n.t("annotation.file.genome"), folder / genome),
            (i18n.t("annotation.file.gtf"), folder / gtf),
            (i18n.t("annotation.file.gene_bed"), folder / Path(info["gene_bed"]).name),
            (i18n.t("annotation.file.splice"), folder / f"splicesites_{index.parent.name}.txt"),
            (i18n.t("annotation.file.exons"), folder / f"exons_{index.parent.name}.txt"),
        ]
        items += [
            (i18n.t("annotation.file.hisat2"), folder / index.parent.name / part)
            for part in RNA_INDEX_PARTS
        ]
        return items
    prefix = Path(info["chip_index"]).name
    items = [(i18n.t("annotation.file.genome"), folder / genome)]
    items += [
        (i18n.t("annotation.file.bowtie2"), folder / f"{prefix}.{part}")
        for part in CHIP_INDEX_PARTS
    ]
    return items


def size_text(value: int) -> str:
    return i18n.size(value)


def dir_size(path: Path) -> int:
    if not path.is_dir():
        return 0
    if path.is_symlink():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.path.getsize(os.path.join(root, name))
            except OSError:
                continue
    return total


def status(kind: str, assembly: str, source: str = "") -> dict:
    items = []
    ready = 0
    for label, path in expected_files(kind, assembly, source):
        exists = path.exists()
        if exists:
            ready += 1
        try:
            size = path.stat().st_size if exists else 0
        except OSError:
            size = 0
        items.append(
            {
                "label": label,
                "name": path.name,
                "path": path,
                "ok": exists,
                "size": size,
                "size_text": size_text(size) if exists else "-",
            }
        )
    folder = kind_dir(kind, assembly)
    return {
        "kind": kind,
        "assembly": assembly,
        "dir": folder,
        "dir_size": size_text(dir_size(folder)),
        "items": items,
        "done": ready,
        "total": len(items),
        "ready": ready == len(items),
        "index_prefix": index_prefix(kind, assembly),
    }


def sources(assembly: str, source: str = "") -> list[tuple[str, str]]:
    info = registry.single(assembly, source)
    base = info["base_url"]
    return [
        (i18n.t("annotation.source.genome"), f"{base}/{Path(info['genome']).name}.gz"),
        (i18n.t("annotation.source.gtf"), f"{base}/{Path(info['gtf']).name}.gz"),
    ]


def build_command(
    kind: str,
    assembly: str,
    threads: int = 8,
    force: bool = False,
    progress: Path | str | None = None,
) -> list[str]:
    command = [
        "bash",
        str(env_paths.BASE_DIR / "scripts" / "build-annotation.sh"),
        "--mine",
        kind,
        "--assembly",
        registry.canonical_assembly(assembly),
        "--threads",
        str(max(1, int(threads))),
    ]
    if force:
        command.append("--force")
    if progress:
        command += ["--progress-file", str(progress)]
    return command


def progress_file(kind: str) -> Path:
    return env_paths.RUN_STATE / "progress" / f"annotation-{kind}.json"


def job_scope(kind: str) -> str:
    """Job slot for one index build, named like its progress file.

    HISAT2 and bowtie2 read different files and write different directories,
    so they get a slot each and can build at the same time.
    """
    return progress_file(kind).stem


def begin_progress(kind: str) -> Path:
    """Drop the previous run's file before starting a new build.

    A build writes its own state; one left over from the last run would be
    drawn as this run's progress until the script overwrites it, which reads
    as a bar that starts at 100%.
    """
    path = progress_file(kind)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        path.unlink()
    return path


def read_progress(kind: str) -> dict:
    if not kind:
        return {}
    path = progress_file(kind)
    if not path.is_file():
        return {}
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return state if isinstance(state, dict) else {}


def progress_view(state: dict) -> dict:
    """Bar fraction, headline and stage rows for one progress file.

    Weights come from the script so the bar spends its time where the build
    does; a file without them still renders, just with flat steps.
    """
    steps = [str(name) for name in (state.get("steps") or [])]
    if not steps:
        return {"overall": 0.0, "text": "", "rows": [], "steps": [], "done": False}
    weights = []
    for weight in (state.get("weights") or []):
        try:
            weights.append(max(0.0, float(weight)))
        except (TypeError, ValueError):
            weights.append(1.0)
    if len(weights) != len(steps):
        weights = [1.0] * len(steps)
    try:
        index = int(state.get("step_index", -1))
    except (TypeError, ValueError):
        index = -1
    try:
        fraction = float(state.get("step_progress") or 0.0)
    except (TypeError, ValueError):
        fraction = 0.0
    fraction = max(0.0, min(1.0, fraction))
    index = max(0, min(index, len(steps) - 1))
    done = bool(state.get("done"))

    if done:
        overall = 1.0
    else:
        total = sum(weights) or 1.0
        overall = (sum(weights[:index]) + weights[index] * fraction) / total

    rows = []
    for position, name in enumerate(steps):
        caption = step_text(name)
        if done or position < index:
            rows.append((caption, "", i18n.t("state.done"), "ok"))
        elif position == index:
            detail = str(state.get("detail") or state.get("message") or "")
            rows.append((caption, detail, i18n.t("state.running"), "run"))
        else:
            rows.append((caption, "", i18n.t("state.idle"), "idle"))

    if done:
        headline = f"{i18n.percent(1)} · {state.get('message') or i18n.t('annotation.built')}"
    else:
        headline = f"{i18n.percent(overall)} · {step_text(steps[index])}"
        detail = str(state.get("detail") or "")
        if detail:
            headline = f"{headline} · {detail}"
    return {"overall": overall, "text": headline, "rows": rows, "steps": steps, "done": done}


def step_text(name: str) -> str:
    key = STEP_KEY.get(name)
    return i18n.t(key) if key else name


def label(kind: str) -> str:
    return KIND_LABEL.get(kind, kind)
