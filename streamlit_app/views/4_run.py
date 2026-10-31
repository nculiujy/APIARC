from __future__ import annotations

import os
import sys
from pathlib import Path

import streamlit as st

_APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import datasets, env_paths, i18n, pipeline, runtime  # noqa: E402
from utils.ui_theme import (  # noqa: E402
    log_console,
    page_header,
    section_title,
    stage_rows,
    status_line,
)

SCOPES = ("RNAseq", "ChIPseq")
config = runtime.load_config()
anno = datasets.annotation_info(config)


def scope_chip(scope: str) -> tuple[str, str]:
    job = runtime.job_status(scope)
    if job.get("running"):
        return i18n.t("run.chip.running", scope=scope), "run"
    return i18n.t("run.chip.idle", scope=scope), "idle"


analysis_text, analysis_kind = runtime.analysis_state()
page_header(
    i18n.t("run.title"),
    chips=[
        (i18n.t("run.chip.analysis", state=analysis_text), analysis_kind),
        (i18n.t("run.chip.annotation", id=anno["id"]), "ok"),
        scope_chip("RNAseq"),
        scope_chip("ChIPseq"),
    ],
)

if not env_paths.analysis_ready():
    st.error(i18n.t("run.error.env"))


def scoped_stages(scope: str) -> tuple[list[str], list[tuple[str, str]]]:
    """Stages that can run now, plus the ones the design filtered out."""
    names = list(pipeline.module_flags(scope, config))
    enabled = pipeline.active_stages(scope, config)
    blocked = [
        (name, pipeline.stage_unavailable(scope, name, config))
        for name in names
        if name not in enabled
    ]
    return enabled, blocked


def outcome(job: dict) -> tuple[str, str]:
    """Status text and tone for a job record.

    A dead pid alone reads as "finished", which hides a snakemake that aborted
    two seconds after launch; the log tail decides between done and broken.
    """
    if job.get("running"):
        return (
            i18n.t(
                "run.outcome.running",
                label=job.get("label"),
                at=i18n.stamp(job.get("started_at")),
            ),
            "run",
        )
    if job.get("stopped_at"):
        return (
            i18n.t(
                "run.outcome.stopped",
                label=job.get("label"),
                at=i18n.stamp(job.get("stopped_at")),
            ),
            "warn",
        )
    tail = runtime.tail_log(job.get("log"), 80)
    if any(token in tail for token in ("Traceback", "Exception", "Error:")):
        return i18n.t("run.outcome.failed", label=job.get("label")), "fail"
    return (
        i18n.t("run.outcome.ended", label=job.get("label"), at=i18n.stamp(job.get("started_at"))),
        "idle",
    )


def launch_panel(scope: str, key: str, job: dict) -> None:
    running = bool(job.get("running"))
    names = list(pipeline.module_flags(scope, config))
    enabled, blocked = scoped_stages(scope)

    cols = st.columns([2, 1])
    labels = {
        name: f"{names.index(name) + 1}. {pipeline.stage_label(scope, name)}"
        for name in enabled
    }
    everything = i18n.t("run.scope.all")
    choice = cols[0].selectbox(
        i18n.t("run.label.range"),
        [everything] + list(labels.values()),
        key=f"scope-{key}",
        disabled=not enabled,
    )
    cores = cols[1].number_input(
        i18n.t("run.label.cores"),
        min_value=1,
        max_value=64,
        value=int(config.get("threads", 8) or 8),
        key=f"cores-{key}",
    )

    if not enabled:
        st.error(i18n.t("run.error.none"))
        return

    st.caption(
        i18n.t("run.caption.all_until", stage=pipeline.stage_label(scope, enabled[-1]))
    )
    for name, reason in blocked:
        st.caption(
            i18n.t(
                "run.caption.blocked",
                stage=pipeline.stage_label(scope, name),
                reason=reason,
            )
        )

    if choice == everything:
        command = pipeline.scope_command(scope, config, cores)
        label = f"{scope}-all"
    else:
        stage = next(name for name, text in labels.items() if text == choice)
        command = pipeline.stage_command(scope, stage, config, cores)
        label = f"{scope}-{stage}"

    cols = st.columns([1, 1, 3])
    if cols[0].button(i18n.t("run.button.run"), type="primary", key=f"go-{key}", disabled=running):
        try:
            runtime.start_job(command, label, scope=scope)
        except RuntimeError as error:
            st.error(str(error))
        else:
            st.rerun()

    if cols[1].button(i18n.t("run.button.stop"), key=f"stop-{key}", disabled=not running):
        runtime.stop_job(scope)
        st.session_state[f"stopped-{key}"] = job.get("label", scope)
        st.rerun()


def module_panel(scope: str, key: str) -> None:
    flags = pipeline.module_flags(scope, config)
    with st.expander(i18n.t("run.section.modules"), expanded=False):
        cols = st.columns(len(flags))
        updated = {}
        for column, (name, enabled) in zip(cols, flags.items()):
            updated[name] = column.checkbox(
                pipeline.stage_label(scope, name), value=enabled, key=f"flag-{key}-{name}"
            )
        if updated != flags:
            section = f"{scope}_modules"
            runtime.update_config({section: updated})
            st.rerun()


def omics_tab(scope: str, key: str) -> None:
    job = runtime.job_status(scope)
    stopped = st.session_state.pop(f"stopped-{key}", "")
    if stopped:
        # A stopped job leaves its stage unfinished, not running, so both labels
        # name the stage snakemake resumes from.
        paused = [
            label
            for label, _detail, status, _kind in pipeline.stage_status(scope, config)
            if status in (i18n.t("state.running"), i18n.t("state.partial"))
        ]
        tail = i18n.t("run.stopped.paused", stage=paused[0]) if paused else ""
        st.info(i18n.t("run.stopped.note", label=stopped, tail=tail))

    name = datasets.config_dataset(config, scope)
    enabled = pipeline.active_stages(scope, config)
    done, total = pipeline.progress(scope, config, only=enabled)
    st.caption(
        i18n.t(
            "run.caption.dataset",
            name=name,
            done=i18n.number(done),
            total=i18n.number(total),
            id=anno["id"],
        )
    )
    versions = [item for item in datasets.result_versions(scope, name) if item != anno["id"]]
    if versions:
        st.caption(i18n.t("run.caption.other_versions", names=", ".join(versions)))

    missing = datasets.missing_fastq(scope, name)
    if missing:
        # A run still listed at step 1 is a download to finish in another tab;
        # only a name outside that list is a table mistake.
        runs = datasets.srr_list(scope, name)
        pending = [item for item in missing if item in runs]
        absent = [item for item in missing if item not in pending]
        if pending:
            st.warning(
                i18n.t(
                    "run.warning.pending",
                    n=i18n.number(len(pending)),
                    names=", ".join(pending[:8]),
                )
            )
        if absent:
            st.warning(
                i18n.t(
                    "run.warning.absent",
                    n=i18n.number(len(absent)),
                    names=", ".join(absent[:8]),
                )
            )

    blocked = {
        pipeline.stage_label(scope, name) for name, _reason in scoped_stages(scope)[1]
    }
    not_applicable = i18n.t("run.na")
    stage_rows(
        [
            (label, detail, not_applicable, "warn")
            if label in blocked
            else (label, detail, status, kind)
            for label, detail, status, kind in pipeline.stage_status(scope, config)
        ]
    )
    module_panel(scope, key)
    if job:
        status_line(*outcome(job))
    section_title(i18n.t("run.section.launch"))
    launch_panel(scope, key, job)


rna_tab, chip_tab = st.tabs(["RNA-seq", "ChIP-seq"])
with rna_tab:
    omics_tab("RNAseq", "rna")
with chip_tab:
    omics_tab("ChIPseq", "chip")

RUN_DRAWERS = 5


def verdict(record: dict, live: bool) -> tuple[str, str]:
    """What one run left behind, read from its own log.

    jobs.json only describes the latest run per slot, so an older drawer has to
    look at the log to tell a completed run from an aborted one.
    """
    if live:
        return i18n.t("run.verdict.live", stamp=record["stamp"]), "run"
    tail = runtime.tail_log(record["path"], 80)
    if any(
        token in tail
        for token in ("Traceback", "Error in rule", "Exiting because", "Exception")
    ):
        return i18n.t("run.verdict.error"), "fail"
    if "Nothing to be done" in tail:
        return i18n.t("run.verdict.uptodate"), "idle"
    if "steps (100%) done" in tail:
        return i18n.t("run.verdict.done"), "ok"
    return i18n.t("run.verdict.ended"), "idle"


def log_groups(scope: str, modules: list[Path]) -> dict[str, list[Path]]:
    """Rule logs of one run, grouped by the stage their number points at."""
    stages = list(pipeline.module_flags(scope, config))
    groups: dict[str, list[Path]] = {}
    for path in modules:
        number = path.name.split("_", 1)[0]
        index = int(number) if number.isdigit() else 0
        if 0 < index <= len(stages):
            label = f"{index}. {pipeline.stage_label(scope, stages[index - 1])}"
        else:
            label = i18n.t("run.other")
        groups.setdefault(label, []).append(path)
    return groups


def run_drawer(record: dict) -> None:
    scope = record["scope"]
    job = runtime.job_status(scope)
    live = bool(job.get("running")) and str(job.get("log")) == str(record["path"])
    text, kind = verdict(record, live)
    # The outcome goes in the drawer title: with ten stacked drawers closed, the
    # one that failed was otherwise only found by opening them all.
    title = i18n.t(
        "run.drawer.title",
        n=i18n.number(record["run"]),
        stamp=record["stamp"],
        label=record["label"],
        verdict=i18n.t(f"run.verdict.short.{kind}"),
        size=i18n.size(record["size"]),
    )
    with st.expander(title, expanded=live):
        status_line(text, kind)
        groups = log_groups(
            scope, runtime.module_logs(scope, record["started"], record["until"])
        )
        whole = i18n.t("run.log.whole")
        options = [whole] + list(groups)
        picked = st.selectbox(
            i18n.t("run.log.pick_group"), options, key=f"log-group-{record['path'].stem}"
        )
        if picked == whole:
            path = record["path"]
        else:
            files = groups[picked]
            if len(files) == 1:
                path = files[0]
            else:
                names = [item.name for item in files]
                name = st.selectbox(
                    i18n.t("run.log.pick_module"), names, key=f"log-file-{record['path'].stem}"
                )
                path = files[names.index(name)]
        try:
            written = i18n.stamp(path.stat().st_mtime)
        except OSError:
            written = i18n.t("common.none")
        st.caption(
            i18n.t("run.log.caption", path=env_paths.relative(path), at=written)
        )
        log_console(runtime.tail_log(path, 240))


head, action = st.columns([5, 1], vertical_alignment="center")
with head:
    section_title(i18n.t("run.section.logs"))
with action:
    if st.button(
        i18n.t("run.logs.clear"),
        type="tertiary",
        key="clear-logs",
        help=i18n.t("run.logs.clear_help"),
    ):
        st.session_state["clear-logs-open"] = True

cleared = st.session_state.pop("clear-logs-done", "")
if cleared:
    st.success(cleared)

if st.session_state.get("clear-logs-open"):
    count, size = runtime.log_footprint()
    if count:
        st.warning(
            i18n.t(
                "run.logs.confirm",
                n=i18n.number(count),
                size=i18n.size(size),
            )
        )
        confirm, cancel, _rest = st.columns([1, 1, 4])
        if confirm.button(i18n.t("run.logs.yes"), key="clear-logs-yes"):
            removed, freed = runtime.clear_logs()
            st.session_state.pop("clear-logs-open", None)
            tail = i18n.t("run.logs.freed", size=i18n.size(freed)) if freed else ""
            st.session_state["clear-logs-done"] = i18n.t(
                "run.logs.done", n=i18n.number(removed), tail=tail
            )
            st.rerun()
        if cancel.button(i18n.t("run.logs.cancel"), type="tertiary", key="clear-logs-no"):
            st.session_state.pop("clear-logs-open", None)
            st.rerun()
    else:
        st.caption(i18n.t("run.logs.none"))
        if st.button(i18n.t("run.logs.close"), type="tertiary", key="clear-logs-close"):
            st.session_state.pop("clear-logs-open", None)
            st.rerun()

runs = runtime.run_records(SCOPES)
parallel = [scope for scope in SCOPES if runtime.is_running(scope)]
if len(parallel) > 1:
    status_line(
        i18n.t(
            "run.parallel",
            n=i18n.number(len(parallel)),
            names=", ".join(parallel),
        ),
        "run",
    )
if not runs:
    st.info(i18n.t("run.empty.runs"))
else:
    for scope in SCOPES:
        mine = [record for record in runs if record["scope"] == scope]
        if not mine:
            continue
        head, action = st.columns([4, 1], vertical_alignment="center")
        head.markdown(
            i18n.t("run.group.head", scope=scope, n=i18n.number(len(mine)))
        )
        show_all = (
            action.toggle(
                i18n.t("run.group.all"),
                value=False,
                key=f"runs-all-{scope}",
                help=i18n.t("run.group.all_help", n=i18n.number(RUN_DRAWERS)),
            )
            if len(mine) > RUN_DRAWERS
            else False
        )
        for record in mine if show_all else mine[:RUN_DRAWERS]:
            run_drawer(record)
