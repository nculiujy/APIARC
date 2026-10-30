from __future__ import annotations

import os
import subprocess
import sys

import streamlit as st

_APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import datasets, env_paths, i18n, pipeline, runtime  # noqa: E402
from utils.ui_theme import (  # noqa: E402
    empty_state,
    log_console,
    page_header,
    section_title,
    stage_rows,
)

job = runtime.job_status("joint")
running = bool(job.get("running"))

config = runtime.load_config()
mode = pipeline.resolved_mode(config)
joint = datasets.joint_id(config)
anno = datasets.annotation_info(config)


def mode_text(value: str) -> str:
    return pipeline.mode_label(value)


page_header(
    i18n.t("joint.title"),
    chips=[
        (i18n.t("joint.chip.dir", id=joint), "idle"),
        (i18n.t("joint.chip.annotation", id=anno["id"]), "ok"),
        (
            i18n.t("joint.chip.mode", mode=mode_text(mode)),
            "ok" if mode == "differential" else "warn",
        ),
        (
            i18n.t(
                "joint.chip.process",
                state=i18n.t("joint.job.running") if running else i18n.t("joint.job.idle"),
            ),
            "run" if running else "idle",
        ),
    ],
)

rna_name = datasets.config_dataset(config, "RNAseq")
chip_name = datasets.config_dataset(config, "ChIPseq")
rna_done = datasets.stage_marker("RNAseq", rna_name, "4_DEseq", config).is_file()
chip_done = datasets.stage_marker("ChIPseq", chip_name, "4_peak_result", config).is_file()
diff_done = datasets.stage_marker("ChIPseq", chip_name, "5_diffbind", config).is_file()


def done_text(flag: bool) -> tuple[str, str]:
    return (
        i18n.t("joint.prereq.done") if flag else i18n.t("joint.prereq.todo"),
        "ok" if flag else "warn",
    )


section_title(i18n.t("joint.section.prereq"))
prereq = [
    (
        i18n.t("joint.prereq.rna"),
        f"{rna_name} · {anno['id']} · 4_deseq",
        *done_text(rna_done),
    ),
    (
        i18n.t("joint.prereq.chip"),
        f"{chip_name} · {anno['id']} · 4_peaks",
        *done_text(chip_done),
    ),
]
if mode == "differential":
    prereq.append(
        (
            i18n.t("joint.prereq.diff"),
            f"{chip_name} · {anno['id']} · 5_diffbind",
            *done_text(diff_done),
        )
    )
stage_rows(prereq, numbered=False)

rna_versions = datasets.result_versions("RNAseq", rna_name)
chip_versions = datasets.result_versions("ChIPseq", chip_name)
if set(rna_versions) != set(chip_versions):
    st.warning(
        i18n.t(
            "joint.warning.versions",
            rna=", ".join(rna_versions) or i18n.t("common.none"),
            chip=", ".join(chip_versions) or i18n.t("common.none"),
            id=anno["id"],
        )
    )

if not (rna_done and chip_done):
    st.warning(i18n.t("joint.warning.prereq"))

done, total = pipeline.progress(
    "Integrated", config, only=pipeline.active_stages("Integrated", config)
)
section_title(i18n.t("joint.section.progress"), i18n.ratio(done, total))
stage_rows(pipeline.stage_status("Integrated", config))

section_title(i18n.t("joint.section.params"))
left, right = st.columns(2)
with left:
    window = st.number_input(
        i18n.t("joint.label.window"),
        min_value=200,
        max_value=50000,
        step=100,
        value=int(config.get("integrated", {}).get("promoter_window", 2000)),
        help=i18n.t("joint.help.window"),
    )
    distance = st.number_input(
        i18n.t("joint.label.distance"),
        min_value=1000,
        max_value=5000000,
        step=1000,
        value=int(config.get("integrated", {}).get("enhancer_max_distance", 1000000)),
        help=i18n.t("joint.help.distance"),
    )
    universe = st.selectbox(
        i18n.t("joint.label.universe"),
        ["expressed", "peaks"],
        index=0
        if config.get("integrated", {}).get("universe", "expressed") == "expressed"
        else 1,
        help=i18n.t("joint.help.universe"),
        format_func=lambda value: i18n.t(f"joint.universe.{value}"),
    )
with right:
    permutations = st.number_input(
        i18n.t("joint.label.permutations"),
        min_value=100,
        max_value=1000000,
        step=100,
        value=int(config.get("stats", {}).get("permutations", 10000)),
        help="stats.permutations",
    )
    seed = st.number_input(
        i18n.t("joint.label.seed"),
        min_value=0,
        max_value=10**9,
        step=1,
        value=int(config.get("stats", {}).get("seed", 20261023)),
        help="stats.seed",
    )
    alpha = st.number_input(
        i18n.t("joint.label.alpha"),
        min_value=0.001,
        max_value=0.5,
        step=0.005,
        format="%.3f",
        value=float(config.get("stats", {}).get("alpha", 0.05)),
        help="stats.alpha",
    )

if st.button(i18n.t("joint.button.save")):
    runtime.update_config(
        {
            "integrated": {
                "promoter_window": int(window),
                "enhancer_max_distance": int(distance),
                "universe": universe,
            },
            "stats": {
                "permutations": int(permutations),
                "seed": int(seed),
                "alpha": float(alpha),
            },
        }
    )
    st.success(i18n.t("joint.saved"))
    st.rerun()

with st.expander(i18n.t("joint.section.modules")):
    flags = pipeline.module_flags("Integrated", config)
    cols = st.columns(len(flags))
    updated = {}
    for column, (name, enabled) in zip(cols, flags.items()):
        updated[name] = column.checkbox(
            pipeline.stage_label("Integrated", name),
            value=enabled,
            key=f"flag-joint-{name}",
        )
    if updated != flags:
        runtime.update_config({"Integrated_modules": updated})
        st.rerun()

section_title(i18n.t("joint.section.launch"))
stopped = st.session_state.pop("joint-stopped", "")
if stopped:
    paused = [
        label
        for label, _detail, status, _kind in pipeline.stage_status("Integrated", config)
        if status == i18n.t("state.running")
    ]
    tail = i18n.t("run.stopped.paused", stage=paused[0]) if paused else ""
    st.info(i18n.t("run.stopped.note", label=stopped, tail=tail))
joint_targets = pipeline.joint_targets(config)
names = list(joint_targets)
# A stage the Snakefile did not compile in cannot be a target: snakemake aborts
# with MissingRuleException before building the DAG, so the reason is shown and
# the run button stays off instead of starting a job that fails at once.
blocked = {
    name: pipeline.stage_unavailable("Integrated", name, config) for name in names
}
cols = st.columns([2, 1])
everything = i18n.t("run.scope.all")
options = [everything] + [
    (
        f"{index}. {pipeline.stage_label('Integrated', name)}"
        + (
            i18n.t("joint.option.skipped")
            if blocked[name]
            else ""
        )
    )
    for index, name in enumerate(names, start=1)
]
choice = cols[0].selectbox(i18n.t("run.label.range"), options, key="joint-scope")
cores = cols[1].number_input(
    i18n.t("run.label.cores"),
    min_value=1,
    max_value=64,
    value=int(config.get("threads", 8) or 8),
    key="joint-cores",
)

reason = ""
if choice == everything:
    command = pipeline.scope_command("Integrated", config, cores)
    label = "joint-all"
else:
    stage = names[int(choice.split(".")[0]) - 1]
    reason = blocked[stage]
    command = pipeline.stage_command("Integrated", stage, config, cores)
    label = f"joint-{stage}"

if reason:
    st.warning(i18n.t("joint.warning.blocked", reason=reason))

cols = st.columns([1, 1, 3])
if cols[0].button(
    i18n.t("run.button.run"),
    type="primary",
    key="joint-go",
    disabled=running or bool(reason),
):
    runtime.start_job(command, label, scope="joint")
    st.rerun()

if cols[1].button(i18n.t("run.button.stop"), key="joint-stop", disabled=not running):
    runtime.stop_job("joint")
    st.session_state["joint-stopped"] = job.get("label", "joint")
    st.rerun()

section_title(i18n.t("joint.section.acceptance"))
if st.button(i18n.t("joint.button.acceptance"), disabled=running):
    command = [
        sys.executable,
        str(env_paths.SCRIPT_ROOT / "common" / "check_acceptance.py"),
        "--config",
        str(env_paths.CONFIG_FILE),
        "--result",
        str(env_paths.RESULT_ROOT),
    ]
    completed = subprocess.run(
        command,
        cwd=str(env_paths.BASE_DIR),
        capture_output=True,
        text=True,
        env=env_paths.routed_env(),
    )
    st.session_state["acceptance"] = {
        "code": completed.returncode,
        "output": (completed.stdout or "") + (completed.stderr or ""),
    }

result = st.session_state.get("acceptance")
if result:
    if result["code"] == 0:
        st.success(i18n.t("joint.acceptance.ok"))
    else:
        st.error(i18n.t("joint.acceptance.failed", code=result["code"]))
    log_console(result["output"], 200)
    report = datasets.provenance_directory(config) / "acceptance_report.tsv"
    if report.is_file():
        st.caption(
            i18n.t(
                "joint.acceptance.report",
                path=os.path.relpath(report, env_paths.BASE_DIR),
            )
        )
        st.dataframe(
            [line.split("\t") for line in report.read_text().strip().splitlines()][1:],
            width="stretch",
            hide_index=True,
        )

section_title(i18n.t("joint.section.logs"))
if job:
    st.caption(
        i18n.t(
            "joint.log.caption",
            label=job.get("label"),
            pid=job.get("pid"),
            at=i18n.stamp(job.get("started_at")),
        )
    )
    cols = st.columns([1, 5])
    if cols[0].button(i18n.t("joint.button.refresh"), key="joint-refresh"):
        st.rerun()
    log_console(runtime.tail_log(job.get("log"), 260))
else:
    empty_state(i18n.t("joint.empty.logs"))
