from __future__ import annotations

from pathlib import Path

from . import annotation, datasets, env_paths, i18n, runtime

STAGE_KEY = {
    "RNAseq": {
        "1_prepare": "stage.common.1_prepare",
        "2_QC": "stage.common.2_QC",
        "3_RC": "stage.rna.3_RC",
        "4_DEseq": "stage.rna.4_DEseq",
    },
    "ChIPseq": {
        "1_prepare": "stage.common.1_prepare",
        "2_QC": "stage.common.2_QC",
        "3_CC": "stage.chip.3_CC",
        "4_peak_result": "stage.chip.4_peak_result",
        "5_diffbind": "stage.chip.5_diffbind",
    },
    "Integrated": {
        "0_stats": "stage.joint.0_stats",
        "1_enrichment": "stage.joint.1_enrichment",
        "2_motif": "stage.joint.2_motif",
        "3_enhancer": "stage.joint.3_enhancer",
        "4_common": "stage.joint.4_common",
    },
}

# Derived from the workspace list so that a stage added there also appears in
# the launcher's run range, the module switches and the marker lookup. A second
# copy here silently dropped 4_common from all three.
JOINT_STAGES = tuple(datasets.INTEGRATED_STAGES)

# jobs.json files the joint run under "joint" while the views name the scope
# "Integrated"; without the mapping those stages never see a live job.
JOB_SCOPE = {"Integrated": "joint"}


def job_scope(scope: str) -> str:
    return JOB_SCOPE.get(scope, scope)


def snakemake_bin() -> str:
    found = env_paths.resolved_tool("snakemake")
    return found or "snakemake"


def base_arguments(cores: int, config_file: Path | None = None) -> list[str]:
    return [
        snakemake_bin(),
        "--cores",
        str(max(1, int(cores))),
        "--configfile",
        str(config_file or env_paths.CONFIG_FILE),
        "--use-conda",
        "--conda-prefix",
        str(env_paths.stage_envs()),
        "--rerun-incomplete",
        "--printshellcmds",
        "--nolock",
    ]


def rna_targets(config: dict) -> dict[str, str]:
    name = datasets.config_dataset(config, "RNAseq")
    return {
        key: str(datasets.stage_marker("RNAseq", name, key, config))
        for key in datasets.STAGE_DIRS["RNAseq"]
    }


def chip_targets(config: dict) -> dict[str, str]:
    name = datasets.config_dataset(config, "ChIPseq")
    return {
        key: str(datasets.stage_marker("ChIPseq", name, key, config))
        for key in datasets.STAGE_DIRS["ChIPseq"]
    }


def joint_targets(config: dict) -> dict[str, str]:
    return {key: str(datasets.joint_marker(config, key)) for key in JOINT_STAGES}


def stage_targets(scope: str, config: dict) -> dict[str, str]:
    """Finished marker of every stage the scope knows, keyed by stage."""
    if scope == "RNAseq":
        return rna_targets(config)
    if scope == "ChIPseq":
        return chip_targets(config)
    return joint_targets(config)


def target_for(scope: str, stage: str, config: dict) -> str:
    return stage_targets(scope, config)[stage]


def stage_command(
    scope: str,
    stage: str,
    config: dict,
    cores: int = 8,
    config_file: Path | None = None,
) -> list[str]:
    # Same reason as scope_command: the rules spell their outputs relative to
    # the repo root, so an absolute target aborts before the DAG is built.
    target = env_paths.relative(target_for(scope, stage, config))
    return base_arguments(cores, config_file) + [target]


def active_stages(scope: str, config: dict) -> list[str]:
    """Stages the Snakefile actually compiled in, in run order.

    Snakemake builds the DAG from the rule files it included, so asking for a
    stage that was switched off - or for 5_diffbind on a single-condition
    design - aborts with MissingRuleException before any job starts. The
    launcher has to compute the same subset: module flags plus the design mode.
    """
    flags = module_flags(scope, config)
    if scope == "Integrated" and not flags.get("0_stats", True):
        # The enrichment, motif, enhancer and common-gene rules are only
        # included on top of 0_stats, so switching it off takes the rest with it.
        return []
    if scope == "Integrated" and flags.get("4_common"):
        # Common genes intersect differential binding with the DEG list. In
        # descriptive mode every region counts as bound, so the rule file is not
        # included and asking for its marker aborts before the DAG is built.
        flags["4_common"] = resolved_mode(config) == "differential"
    if scope == "ChIPseq" and flags.get("5_diffbind"):
        flags["5_diffbind"] = resolved_mode(config) == "differential"
    return [name for name in flags if flags[name]]


def stage_unavailable(scope: str, stage: str, config: dict) -> str:
    """Why a stage listed in the UI will not run now, or '' when it will."""
    if stage in active_stages(scope, config):
        return ""
    if scope == "Integrated" and stage == "4_common":
        return i18n.t("pipe.unavailable.common")
    if scope == "ChIPseq" and stage == "5_diffbind":
        return i18n.t("pipe.unavailable.comparisons")
    return i18n.t("pipe.unavailable.off")


def scope_command(
    scope: str,
    config: dict,
    cores: int = 8,
    config_file: Path | None = None,
) -> list[str]:
    stages = active_stages(scope, config)
    if not stages:
        raise ValueError(i18n.t("pipe.no_active_stage", scope=scope))
    # Every enabled marker is a target, not only the deepest one. A stage's
    # finished flag is no input of the next stage's first rule - 5_diffbind
    # reads the peak BEDs directly - so a run to 5_diffbind never executed the
    # annotation, coverage and QC rules of stage 4 and never wrote its marker,
    # which the pages then reported as a stage still running. The same holds
    # for the joint stages: enrichment, motif and enhancer each read what
    # 0_stats wrote, none of them feeds another, so naming only the last one
    # built a three-rule DAG and skipped the outputs the results page reads.
    markers = stage_targets(scope, config)
    targets = [env_paths.relative(markers[name]) for name in stages]
    # The rules spell their outputs relative to the repo root; snakemake 7
    # refuses an absolute target against them and aborts with
    # MissingRuleException before the DAG is built.
    return base_arguments(cores, config_file) + targets


def stage_label(scope: str, name: str) -> str:
    key = STAGE_KEY.get(scope, {}).get(name)
    return i18n.t(key) if key else name


def stage_markers(
    scope: str, config: dict, markers: dict[str, Path] | None = None
) -> list[tuple[str, str, Path]]:
    """Stage key, folder detail and finished marker for one scope."""
    if scope == "Integrated":
        folders = datasets.ws.INTEGRATED_STAGES
        # Keyed off the workspace's own stage list: a stage added there shows up
        # here without a second list to keep in step.
        defaults = {key: datasets.joint_marker(config, key) for key in folders}
    else:
        name = datasets.config_dataset(config, scope)
        folders = datasets.STAGE_DIRS[scope]
        defaults = {
            key: datasets.stage_marker(scope, name, key, config) for key in folders
        }
    return [
        (key, detail, (markers or {}).get(key) or defaults[key])
        for key, detail in folders.items()
    ]


def stage_status(
    scope: str,
    config: dict,
    markers: dict[str, Path] | None = None,
    running: bool | None = None,
) -> list[tuple[str, str, str, str]]:
    """One row per stage: label, folder detail, status, status kind.

    ``markers`` overrides where the finished markers live, so a page can show a
    dataset or annotation other than the configured one without touching config.
    ``running`` overrides whether the scope's job is alive; by default it comes
    from the job record, but only for the configured dataset - a page showing an
    older one must not inherit a live run on a different result directory.
    """
    if running is None:
        running = markers is None and runtime.is_running(job_scope(scope))
    return [
        (i18n.t(STAGE_KEY[scope][key]), detail, *_marker_state(marker, running))
        for key, detail, marker in stage_markers(scope, config, markers)
    ]


def _marker_state(marker: Path, running: bool = False) -> tuple[str, str]:
    """Status of one stage, decided by its marker and by the job record.

    Outputs without a marker mean the stage was interrupted or never asked for,
    so only a live process makes it "running"; reading that off the folder left
    every finished run reporting "running" on the stage it skipped.
    """
    if marker.is_file():
        return i18n.t("state.done"), "ok"
    if marker.parent.is_dir() and any(marker.parent.iterdir()):
        if running:
            return i18n.t("state.running"), "run"
        return i18n.t("state.partial"), "warn"
    return i18n.t("state.idle"), "idle"


def progress(
    scope: str,
    config: dict,
    markers: dict[str, Path] | None = None,
    only: list[str] | None = None,
) -> tuple[int, int]:
    """Finished stages out of the stages that can run.

    ``only`` drops stages the Snakefile never compiled in, otherwise a
    single-condition ChIP-seq design would show 4/5 forever.
    """
    rows = stage_markers(scope, config, markers)
    if only is not None:
        allowed = set(only)
        rows = [row for row in rows if row[0] in allowed]
    done = sum(1 for _key, _detail, marker in rows if marker.is_file())
    return done, len(rows)


def mode_label(mode: str) -> str:
    """The integration mode in the reader's language; an unknown one stays as is."""
    key = f"joint.mode.{mode}"
    text = i18n.t(key)
    return mode if text == key else text


def resolved_mode(config: dict) -> str:
    try:
        from common import metadata as meta

        declared = config.get("integrated", {}).get("mode", "auto")
        if declared != "auto":
            return declared
        chip_name = datasets.config_dataset(config, "ChIPseq")
        root = env_paths.DATA_ROOT
        samples = meta.read_chip_metadata(str(datasets.ws.metadata_path("ChIPseq", chip_name, root)))
        # The five-column schema decides from comparisons.csv, the file the
        # Snakefile also hands to resolve_mode. Without it every such dataset
        # looked descriptive here and 5_diffbind dropped out of the launcher.
        pairs = datasets.ws.comparisons_path("ChIPseq", chip_name, root)
        return meta.resolve_mode(config, samples, pairs)
    except Exception:  # noqa: BLE001
        return "unknown"


def module_flags(scope: str, config: dict) -> dict[str, bool]:
    key = f"{scope}_modules"
    section = config.get(key) or {}
    if scope == "RNAseq":
        names = list(datasets.STAGE_DIRS["RNAseq"])
    elif scope == "ChIPseq":
        names = list(datasets.STAGE_DIRS["ChIPseq"])
    else:
        names = list(JOINT_STAGES)
    return {name: bool(section.get(name, True)) for name in names}


def annotation_issues(config: dict) -> list[tuple[str, str, str]]:
    """Reference files the current annotation needs before a run can start."""
    info = annotation.resolve(config)
    issues: list[tuple[str, str, str]] = []
    if not Path(f"{info['rna_index']}.1.ht2").is_file():
        issues.append(
            (i18n.t("pipe.issue.hisat2", id=info["id"]), i18n.t("pipe.issue.hisat2_hint"), "warn")
        )
    if not Path(f"{info['chip_index']}.1.bt2").is_file():
        issues.append(
            (i18n.t("pipe.issue.bowtie2", id=info["id"]), i18n.t("pipe.issue.bowtie2_hint"), "warn")
        )
    if not Path(info["gene_bed"]).is_file():
        issues.append(
            (i18n.t("pipe.issue.gene_bed", id=info["id"]), i18n.t("pipe.issue.gene_bed_hint"), "warn")
        )
    return issues


def preflight(config: dict) -> list[tuple[str, str, str]]:
    issues: list[tuple[str, str, str]] = []
    issues.extend(annotation_issues(config))
    for omics in datasets.OMICS:
        name = datasets.config_dataset(config, omics)
        base = datasets.directory(omics, name)
        if not base.is_dir():
            issues.append(
                (
                    i18n.t("pipe.issue.dataset", omics=omics),
                    f"workflow/data/{omics}/{name}",
                    "fail",
                )
            )
            continue
        if not (base / datasets.META).is_file():
            issues.append(
                (
                    i18n.t("pipe.issue.metadata", omics=omics),
                    f"workflow/data/{omics}/{name}/metadata.csv",
                    "fail",
                )
            )
        missing = datasets.missing_fastq(omics, name)
        if missing:
            issues.append(
                (
                    i18n.t("pipe.issue.fastq", omics=omics),
                    i18n.t("pipe.issue.fastq_hint", n=i18n.number(len(missing))),
                    "warn",
                )
            )
    if not runtime.ready_for_run():
        issues.append(
            (i18n.t("pipe.issue.analysis"), i18n.t("pipe.issue.analysis_hint"), "warn")
        )
    return issues
