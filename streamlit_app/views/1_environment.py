from __future__ import annotations

import os
import sys

import pandas as pd
import streamlit as st

_APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import env_paths, i18n, runtime  # noqa: E402
from utils.ui_theme import (  # noqa: E402
    empty_state,
    key_values,
    metrics,
    page_header,
    section_title,
    status_line,
)

launcher_text, launcher_kind = runtime.launcher_state()
analysis_text, analysis_kind = runtime.analysis_state()


def dependency_frame(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                i18n.t("env.column.tool"): row["package"],
                i18n.t("env.column.declared"): row["pinned"] or i18n.t("common.none"),
                i18n.t("env.column.installed"): row["version"],
                i18n.t("env.column.status"): runtime.status_text(row["status"]),
                i18n.t("env.column.command"): row["command"],
            }
            for row in rows
        ]
    )


page_header(
    i18n.t("env.title"),
    i18n.t("env.lead"),
    chips=[
        (i18n.t("env.chip.launcher", state=launcher_text), launcher_kind),
        (i18n.t("env.chip.analysis", state=analysis_text), analysis_kind),
        (i18n.t("env.chip.stages", n=i18n.number(runtime.stage_env_count())), "idle"),
    ],
)

section_title(i18n.t("env.section.runtime"))
key_values(
    [
        (i18n.t("env.row.launcher"), env_paths.relative(env_paths.LAUNCH_PREFIX)),
        (i18n.t("env.row.analysis"), env_paths.relative(env_paths.analysis_prefix())),
        (i18n.t("env.row.stages"), env_paths.relative(env_paths.stage_envs())),
        (i18n.t("env.row.definitions"), env_paths.relative(env_paths.BASE_DIR / "environment.yml")),
        (i18n.t("env.row.python"), env_paths.relative(env_paths.launcher_python())),
        ("snakemake", env_paths.resolved_tool("snakemake") or i18n.t("common.none")),
    ]
)

report = runtime.dependency_report()
summary = runtime.dependency_summary(report)
base_rows = report.get("environment.yml", [])
stage_rows = {label: rows for label, rows in report.items() if label != "environment.yml"}

section_title(i18n.t("env.section.dependencies"), i18n.t("env.section.dependencies_hint"))
metrics(
    [
        (i18n.t("env.metric.declared"), i18n.number(summary["base_total"]), "environment.yml"),
        (
            i18n.t("env.metric.readable"),
            i18n.number(summary["base_ready"]),
            i18n.t("env.metric.readable_hint"),
        ),
        (i18n.t("env.metric.absent"), i18n.number(summary["base_absent"]), ""),
        (i18n.t("env.metric.broken"), i18n.number(summary["broken"]), ""),
    ]
)
if summary["base_absent"]:
    status_line(
        i18n.t("env.status.absent", n=i18n.number(summary["base_absent"])), "fail"
    )
else:
    status_line(i18n.t("env.status.ok", n=i18n.number(summary["base_total"])), "ok")

if not base_rows:
    empty_state(i18n.t("env.empty.base"))
else:
    st.dataframe(dependency_frame(base_rows), hide_index=True, width="stretch")


def remember_stage() -> None:
    runtime.save_pref("stage_env", st.session_state.get("env-stage", ""))


section_title(i18n.t("env.section.stages"))
st.caption(
    i18n.t(
        "env.stages.summary",
        total=i18n.number(summary["stage_total"]),
        ready=i18n.number(summary["stage_ready"]),
        pending=i18n.number(summary["stage_pending"]),
    )
)
if stage_rows:
    names = list(stage_rows)
    defaults = [
        name
        for name in names
        if any(row["status"] != runtime.STATUS_READABLE for row in stage_rows[name])
    ]
    default = defaults[0] if defaults else names[0]
    saved = runtime.read_pref("stage_env")
    start = names.index(saved) if saved in names else names.index(default)
    stage = st.selectbox(
        i18n.t("env.stage.picker"),
        names,
        index=start,
        key="env-stage",
        on_change=remember_stage,
    )
    rows = stage_rows[stage]
    readable = sum(1 for row in rows if row["status"] == runtime.STATUS_READABLE)
    status_line(
        i18n.t("env.stage.line", stage=stage, ready=readable, total=len(rows)),
        "ok" if readable == len(rows) else "idle",
    )
    st.dataframe(dependency_frame(rows), hide_index=True, width="stretch")
else:
    empty_state(i18n.t("env.empty.stages"))

section_title(i18n.t("env.section.panel"))
panel_rows = [
    {
        "name": item["name"],
        "version": item["version"],
        "status": item["status"],
        "where": item["where"],
    }
    for item in runtime.panel_report()
]
st.dataframe(
    pd.DataFrame(
        [
            {
                i18n.t("env.column.dependency"): row["name"],
                i18n.t("env.column.version"): row["version"],
                i18n.t("env.column.status"): runtime.status_text(row["status"]),
                i18n.t("env.column.module"): row["where"],
            }
            for row in panel_rows
        ]
    ),
    hide_index=True,
    width="stretch",
)
missing_panel = [
    row["name"] for row in panel_rows if row["status"] != runtime.STATUS_READABLE
]
if missing_panel:
    st.warning(i18n.t("env.panel.missing", names=", ".join(missing_panel)))

section_title(i18n.t("env.section.rescue"))
base_missing = [row for row in base_rows if row["status"] != runtime.STATUS_READABLE]
stage_missing = [
    row
    for rows in stage_rows.values()
    for row in rows
    if row["status"] != runtime.STATUS_READABLE
]
if not base_missing and not stage_missing:
    status_line(i18n.t("env.rescue.ok"), "ok")
if base_missing:
    st.markdown(
        i18n.t("env.rescue.base", names=", ".join(f"`{row['package']}`" for row in base_missing))
    )
    specs = [
        f'"{row["package"]}{("=" + row["pinned"]) if row["pinned"] else ""}"'
        for row in base_missing
    ]
    st.code(
        f"conda install -p {env_paths.relative(env_paths.analysis_prefix())} "
        "-c conda-forge -c bioconda " + " ".join(specs),
        language="bash",
    )
if stage_missing:
    st.markdown(i18n.t("env.rescue.stages", n=i18n.number(len(stage_missing))))
    st.code("apiarc-app-update analysis", language="bash")
st.caption(i18n.t("env.rescue.full"))
