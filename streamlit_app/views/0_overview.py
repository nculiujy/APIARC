from __future__ import annotations

import os
import sys

import streamlit as st

_APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import annotation, datasets, env_paths, i18n, pipeline, prepare, runtime  # noqa: E402
from utils.ui_theme import (  # noqa: E402
    empty_state,
    list_rows,
    log_console,
    metrics,
    page_header,
    section_title,
    stage_rows,
)
from utils.version import APP_FULL_NAME  # noqa: E402

config = runtime.load_config()
running = runtime.running_jobs()
user = runtime.current_user()

rna_name = datasets.config_dataset(config, "RNAseq")
chip_name = datasets.config_dataset(config, "ChIPseq")
joint = datasets.joint_id(config)
mode = pipeline.resolved_mode(config)
anno = datasets.annotation_info(config)
issues = pipeline.preflight(config)
blocking = [item for item in issues if item[2] == "fail"]

head, switch = st.columns([5, 1.2], vertical_alignment="center")
with head:
    page_header(
        i18n.t("overview.title", app=APP_FULL_NAME),
        chips=[
            (i18n.t("overview.chip.rna", name=rna_name), "idle"),
            (i18n.t("overview.chip.chip", name=chip_name), "idle"),
            (i18n.t("overview.chip.annotation", id=anno["id"]), "ok"),
            (
                i18n.t("overview.chip.integration", mode=pipeline.mode_label(mode)),
                "ok" if mode != "unknown" else "warn",
            ),
            (
                i18n.t(
                    "overview.chip.admin" if user["admin"] else "overview.chip.user",
                    name=user["name"],
                    uid=user["uid"],
                ),
                "warn" if user["admin"] else "ok",
            ),
        ],
    )
with switch:
    i18n.language_switch()

rna_done, rna_total = pipeline.progress(
    "RNAseq", config, only=pipeline.active_stages("RNAseq", config)
)
chip_done, chip_total = pipeline.progress(
    "ChIPseq", config, only=pipeline.active_stages("ChIPseq", config)
)
joint_done, joint_total = pipeline.progress("Integrated", config)

metrics(
    [
        (
            i18n.t("overview.metric.rna"),
            i18n.ratio(rna_done, rna_total),
            f"{rna_name} · {anno['id']}",
        ),
        (
            i18n.t("overview.metric.chip"),
            i18n.ratio(chip_done, chip_total),
            f"{chip_name} · {anno['id']}",
        ),
        (
            i18n.t("overview.metric.joint"),
            i18n.ratio(joint_done, joint_total),
            f"{joint} · {anno['id']}",
        ),
        (
            i18n.t("overview.metric.jobs"),
            i18n.number(len(running)),
            " + ".join(job.get("label", "") for job in running) or i18n.t("common.idle"),
        ),
    ]
)

if blocking:
    for title, detail, _kind in blocking:
        st.error(f"{title} — {detail}")
elif issues:
    for title, detail, _kind in issues:
        st.warning(f"{title} — {detail}")
else:
    st.success(i18n.t("overview.ok"))

section_title(i18n.t("overview.section.progress"))
left, right = st.columns(2)
with left:
    st.markdown(f"**{i18n.t('overview.omics.rna')}**")
    stage_rows(pipeline.stage_status("RNAseq", config))
with right:
    st.markdown(f"**{i18n.t('overview.omics.chip')}**")
    stage_rows(pipeline.stage_status("ChIPseq", config))

st.markdown(f"**{i18n.t('overview.omics.joint')}**")
stage_rows(pipeline.stage_status("Integrated", config))

section_title(i18n.t("overview.section.config"))


@st.cache_data(ttl=20, show_spinner=False)
def env_summary() -> dict:
    """Reading every conda-meta file costs most of a second; the overview is
    reloaded far more often than the environment changes."""
    return runtime.dependency_summary()


def env_rows() -> tuple[bool, list[tuple[str, str, str, str]]]:
    summary = env_summary()
    panel = runtime.panel_report()
    panel_bad = [item["name"] for item in panel if item["status"] != runtime.STATUS_READABLE]
    rows = [
        (
            i18n.t("overview.row.env.base"),
            i18n.t(
                "overview.row.env.base_detail",
                path=env_paths.relative(env_paths.analysis_prefix()),
                ready=i18n.number(summary["base_ready"]),
                total=i18n.number(summary["base_total"]),
            ),
            i18n.t("overview.state.configured")
            if not summary["base_absent"]
            else i18n.t("overview.state.missing", n=i18n.number(summary["base_absent"])),
            "ok" if not summary["base_absent"] else "fail",
        ),
        (
            i18n.t("overview.row.env.stages"),
            i18n.t(
                "overview.row.env.stages_detail",
                built=i18n.number(runtime.stage_env_count()),
                declared=i18n.number(summary["stage_total"]),
            ),
            i18n.t("overview.state.configured")
            if not summary["stage_pending"]
            else i18n.t("overview.state.pending", n=i18n.number(summary["stage_pending"])),
            "ok" if not summary["stage_pending"] else "idle",
        ),
        (
            i18n.t("overview.row.env.panel"),
            i18n.t(
                "overview.row.env.panel_detail",
                ready=i18n.number(len(panel) - len(panel_bad)),
                total=i18n.number(len(panel)),
            ),
            i18n.t("overview.state.configured")
            if not panel_bad
            else i18n.t("overview.state.missing", n=i18n.number(len(panel_bad))),
            "ok" if not panel_bad else "warn",
        ),
    ]
    ok = not summary["base_absent"] and not summary["stage_pending"] and not panel_bad
    return ok, rows


def annotation_rows() -> tuple[bool, list[tuple[str, str, str, str]]]:
    rows = [
        (
            i18n.t("overview.row.reference"),
            f"{anno['species']} · {anno['assembly']} · {anno['source']}",
            i18n.t("overview.state.selected"),
            "ok",
        )
    ]
    ok = True
    for kind, label in (("rna", "RNA-seq · HISAT2"), ("chip", "ChIP-seq · bowtie2")):
        state = annotation.status(kind, anno["assembly"], anno["source"])
        ok = ok and state["ready"]
        rows.append(
            (
                label,
                i18n.t(
                    "overview.row.aligner_detail",
                    name=state["dir"].name,
                    done=i18n.number(state["done"]),
                    total=i18n.number(state["total"]),
                ),
                i18n.t("overview.state.deployed" if state["ready"] else "overview.state.undeployed"),
                "ok" if state["ready"] else "warn",
            )
        )
    return ok, rows


def prepare_rows() -> tuple[bool, list[tuple[str, str, str, str]]]:
    rows = [(i18n.t("overview.row.project"), datasets.joint_id(config), "", "idle")]
    ok = True
    for omics, label_key in (
        ("RNAseq", "overview.row.rna_table"),
        ("ChIPseq", "overview.row.chip_table"),
    ):
        name = datasets.config_dataset(config, omics)
        ready = bool(name) and prepare.prepared(omics, name)
        ok = ok and ready
        rows.append(
            (
                i18n.t(label_key),
                name,
                i18n.t("common.configured") if ready else i18n.t("common.not_configured"),
                "ok" if ready else "warn",
            )
        )
    rows.append(
        (
            i18n.t("overview.row.annotation"),
            f"{anno['species']} · {anno['source']}",
            i18n.t("overview.state.selected"),
            "ok",
        )
    )
    return ok, rows


blocks = [
    (
        i18n.t("overview.block.env"),
        *env_rows(),
        i18n.t("overview.block.env.ok"),
        i18n.t("overview.block.env.bad"),
    ),
    (
        i18n.t("overview.block.annotation"),
        *annotation_rows(),
        i18n.t("overview.block.annotation.ok"),
        i18n.t("overview.block.annotation.bad"),
    ),
    (
        i18n.t("overview.block.prepare"),
        *prepare_rows(),
        i18n.t("overview.block.prepare.ok"),
        i18n.t("overview.block.prepare.bad"),
    ),
]
incomplete = [name for name, ok, _rows, _good, _bad in blocks if not ok]
for name, ok, rows, good, bad in blocks:
    with st.expander(
        f"{name} · {good if ok else bad}",
        expanded=bool(incomplete) and name == incomplete[0],
    ):
        stage_rows(rows, numbered=False)

section_title(i18n.t("overview.section.recent"))
logs = runtime.recent_logs(5)
if not logs:
    empty_state(i18n.t("overview.empty.logs"))
else:
    list_rows(
        [
            (
                path.name,
                f"{i18n.stamp(path.stat().st_mtime)} · {i18n.size(path.stat().st_size)}",
                "",
                "idle",
            )
            for path in logs
        ]
    )
    with st.expander(i18n.t("overview.logs.expand")):
        log_console(runtime.tail_log(logs[0], 120))
