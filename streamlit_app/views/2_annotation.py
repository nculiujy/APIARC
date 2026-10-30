from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import pandas as pd
import streamlit as st

_APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import annotation, env_paths, i18n, motif_resources, runtime  # noqa: E402
from utils.ui_theme import (  # noqa: E402
    elapsed_text,
    empty_state,
    key_values,
    log_console,
    metrics,
    page_header,
    section_title,
    stage_rows,
    status_line,
)

KINDS = ("rna", "chip")
KIND_LABEL = {"rna": "HISAT2", "chip": "bowtie2"}


def current_jobs() -> dict[str, dict]:
    """Job per kind, including one started before the slots were split.

    Builds launched by the older page were recorded under the single
    ``annotation`` scope. Adopting those keeps the stop button, the log and the
    disabled build button attached to a build that is still running.
    """
    legacy = runtime.job_status("annotation")
    found: dict[str, dict] = {}
    for kind in KINDS:
        record = runtime.job_status(annotation.job_scope(kind))
        if not record and str(legacy.get("label") or "").startswith(f"annotation-{kind}-"):
            record = legacy
        if record:
            found[kind] = record
    return found


def kind_state(kind: str) -> tuple[str, str]:
    aligner = annotation.KIND_ALIGNER[kind]
    ready = env_paths.resolved_tool(aligner) is not None
    if ready:
        return i18n.t("anno.kind.ready", tool=aligner), "ok"
    return i18n.t("anno.kind.missing", tool=aligner), "warn"


def genome_writer() -> str:
    """The kind currently downloading the shared genome FASTA, if any.

    Both builds read ``RNAseq_anno/<assembly>/`` — the ChIP side only symlinks
    the genome — so the second build must not start while the first is still
    writing it: ``download`` reuses any non-empty file, and gunzipping half an
    archive leaves a truncated FASTA that later builds happily reuse.
    """
    for kind in KINDS:
        record = jobs.get(kind) or {}
        if not record.get("running"):
            continue
        state = annotation.read_progress(kind)
        if state.get("step") == "genome" and not state.get("done"):
            return kind
    return ""


jobs = current_jobs()
running = {kind: bool(jobs.get(kind, {}).get("running")) for kind in KINDS}
busy = any(running.values())
config = runtime.load_config()
current = annotation.resolve(config)
# While a build runs, the bar belongs in the header: the job panel sits below
# both tabs and would otherwise only be visible after scrolling.
build_views = {
    kind: annotation.progress_view(annotation.read_progress(kind)) if running[kind] else {}
    for kind in KINDS
}

chips = [
    (i18n.t("anno.chip.current", id=current["id"]), "ok"),
    kind_state("rna"),
    kind_state("chip"),
    (
        i18n.t(
            "anno.chip.job",
            state=i18n.t("anno.job.running") if busy else i18n.t("anno.job.idle"),
        ),
        "run" if busy else "idle",
    ),
]
for kind in KINDS:
    if build_views[kind].get("steps"):
        chips.append(
            (
                i18n.t(
                    "anno.chip.progress",
                    kind=KIND_LABEL[kind],
                    percent=i18n.percent(build_views[kind]["overall"]),
                ),
                "run",
            )
        )

page_header(i18n.t("anno.title"), chips=chips)

for kind in KINDS:
    if build_views[kind].get("steps"):
        st.progress(
            build_views[kind]["overall"],
            text=f"{KIND_LABEL[kind]} · {build_views[kind]['text']}",
        )

section_title(i18n.t("anno.section.project"))
assemblies = annotation.assembly_names()
left, middle = st.columns([1, 1])
assembly = left.selectbox(
    i18n.t("anno.label.assembly"),
    assemblies,
    index=assemblies.index(current["assembly"]),
    key="anno-assembly",
)
available_sources = annotation.source_names(assembly)
source = middle.selectbox(
    i18n.t("anno.label.source"),
    available_sources,
    index=available_sources.index(current["source"]) if current["source"] in available_sources else 0,
    key="anno-source",
)
st.caption(i18n.t("anno.caption.assemblies"))
if assembly != current["assembly"] or source != current["source"]:
    st.warning(i18n.t("anno.warning.switch", id=annotation.annotation_id(assembly, source)))
    if st.button(i18n.t("anno.button.apply"), type="primary", key="anno-apply"):
        runtime.update_config({"annotation": {"assembly": assembly, "source": source}})
        st.rerun()

# The four metrics below describe the version picked in the two selectboxes,
# not the one already written to config: switching to GRCh38 has to show the
# human TxDb immediately, otherwise the page says mouse while the file list
# underneath already lists the human directory.
info = annotation.resolve(
    runtime.deep_merge(
        config, {"annotation": {"assembly": assembly, "source": source}}
    )
)
with st.expander(i18n.t("anno.expander.identity", id=info["id"])):
    key_values(
        [
            (i18n.t("anno.row.id"), info["id"]),
            (i18n.t("anno.row.assembly"), f"{info['assembly']} · {info['species']}"),
            (i18n.t("anno.row.source"), f"{info['source']} · {info['release']}"),
            (i18n.t("anno.row.peak_db"), f"{info['txdb']} · {info['orgdb']}"),
        ]
    )

section_title(i18n.t("anno.motif.section"), i18n.t("anno.motif.hint"))
note = st.session_state.pop("motif-note", "")
if note:
    st.success(note)
motif = motif_resources.status()
motif_status = i18n.t(
    "anno.motif.status",
    release=motif["release"] or motif_resources.RELEASE["label"],
    n=i18n.number(motif["matrices"]),
)
if motif["downloaded"]:
    motif_status = i18n.t("anno.motif.downloaded", line=motif_status, n=i18n.number(motif["downloaded"]))
status_line(motif_status, "ok" if motif["ready"] else "warn")
# The counts and paths sit in one drawer: on a ready install the page only needs
# the pill and the button, and the four cards pushed the two build blocks down a
# screen for numbers nobody reads twice.
with st.expander(i18n.t("anno.motif.expander"), expanded=not motif["ready"]):
    metrics(
        [
            (
                i18n.t("anno.motif.metric.matrices"),
                i18n.number(motif["matrices"]),
                motif_resources.MATRIX_DIR,
            ),
            (
                i18n.t("anno.motif.metric.size"),
                i18n.size(motif["bytes"]) if motif["bytes"] else i18n.t("common.none"),
                "",
            ),
            (
                i18n.t("anno.motif.metric.map"),
                i18n.t("anno.motif.metric.map_rows", n=i18n.number(motif["map_rows"])),
                motif_resources.MAP_NAME,
            ),
            (
                i18n.t("anno.motif.metric.tf"),
                f"{i18n.number(motif['tf']['homo'])} / {i18n.number(motif['tf']['mm'])}",
                i18n.t("anno.motif.metric.tf_hint"),
            ),
        ]
    )
    key_values(
        [
            (
                i18n.t("anno.motif.row.matrix_dir"),
                str(motif_resources.RESOURCE_DIR / motif_resources.MATRIX_DIR),
            ),
            (
                i18n.t("anno.motif.metric.map"),
                str(motif_resources.RESOURCE_DIR / motif_resources.MAP_NAME),
            ),
            (i18n.t("anno.motif.row.tf_lists"), ", ".join(motif_resources.TF_LISTS.values())),
            (
                i18n.t("anno.motif.row.source"),
                motif["source"] or i18n.t("anno.motif.source.unrecorded"),
            ),
            (
                i18n.t("anno.motif.row.stamp"),
                i18n.stamp(motif["stamp"]) if motif["stamp"] else i18n.t("common.none"),
            ),
        ]
    )
    st.caption(i18n.t("anno.motif.caption"))
    st.markdown(
        i18n.t(
            "anno.motif.bullet.matrix",
            label=motif_resources.RELEASE["label"],
            url=motif_resources.RELEASE["url"],
        )
    )
    st.markdown(i18n.t("anno.motif.bullet.tf"))
if not motif["ready"]:
    st.warning(i18n.t("anno.motif.missing", dir=motif_resources.MATRIX_DIR))
release = motif_resources.RELEASE["label"]
if st.button(
    i18n.t("anno.motif.redownload" if motif["ready"] else "anno.motif.download", release=release),
    key="motif-download",
    type="secondary" if motif["ready"] else "primary",
):
    bar = st.progress(
        0.0,
        text=i18n.t("anno.motif.connect", host=motif_resources.RELEASE["url"].split("/")[2]),
    )

    def tick(fraction: float, done: int = 0, total: int = 0) -> None:
        if total:
            label = i18n.t(
                "anno.motif.downloading",
                done=i18n.size(done),
                total=i18n.size(total),
            )
        elif fraction >= 1:
            label = i18n.t("anno.motif.unpacking")
        else:
            label = i18n.t("anno.motif.retrying")
        bar.progress(fraction, text=label)

    try:
        result = motif_resources.install(progress=tick)
    except motif_resources.DownloadError as error:
        bar.empty()
        st.error(str(error))
    else:
        bar.empty()
        st.session_state["motif-note"] = i18n.t(
            "anno.motif.installed",
            n=i18n.number(result["matrices"]),
            size=i18n.size(result["bytes"]),
        )
        st.rerun()


def render_status(kind: str) -> dict:
    info = annotation.status(kind, assembly, source)
    status_line(
        i18n.t(
            "anno.status.line",
            release=annotation.meta(assembly, source)["release"],
            done=i18n.number(info["done"]),
            total=i18n.number(info["total"]),
        ),
        "ok" if info["ready"] else ("warn" if info["done"] else "idle"),
    )
    metrics(
        [
            (
                i18n.t("anno.metric.dir"),
                str(info["dir"].relative_to(annotation.anno_root())),
                "",
            ),
            (i18n.t("anno.metric.size"), info["dir_size"], ""),
            (
                i18n.t("anno.metric.files"),
                i18n.ratio(info["done"], info["total"]),
                "",
            ),
            (
                i18n.t("anno.metric.prefix"),
                str(info["index_prefix"].relative_to(info["dir"])),
                "",
            ),
        ]
    )
    # HISAT2 ships eight index shards and bowtie2 six: one row each printed
    # "HISAT2 index / built" eight times and pushed the genome and GTF rows off
    # screen, so the shards of one kind share a row.
    groups: dict[str, dict] = {}
    for item in info["items"]:
        entry = groups.setdefault(item["label"], {"files": [], "done": 0, "size": 0})
        entry["files"].append(item["name"])
        entry["done"] += 1 if item["ok"] else 0
        entry["size"] += item["size"]
    frame = pd.DataFrame(
        [
            {
                i18n.t("anno.column.type"): name,
                i18n.t("anno.column.file"): entry["files"][0]
                + (
                    i18n.t("anno.and_more", n=i18n.number(len(entry["files"])))
                    if len(entry["files"]) > 1
                    else ""
                ),
                i18n.t("anno.column.size"): i18n.size(entry["size"])
                if entry["size"]
                else i18n.t("common.none"),
                i18n.t("anno.column.status"): i18n.t("anno.generated")
                if entry["done"] == len(entry["files"])
                else i18n.ratio(entry["done"], len(entry["files"])),
            }
            for name, entry in groups.items()
        ]
    )
    st.dataframe(frame, hide_index=True, width="stretch")
    missing = [item["name"] for item in info["items"] if not item["ok"]]
    if missing:
        st.caption(i18n.t("anno.missing_files", names=", ".join(missing)))
    return info


def render_build(kind: str, key: str) -> None:
    aligner = annotation.KIND_ALIGNER[kind]
    threads = int(config.get("threads", 8) or 8)
    cols = st.columns([1, 1, 2])
    cores = cols[0].number_input(
        i18n.t("anno.build.threads"),
        min_value=1,
        max_value=64,
        value=threads,
        key=f"threads-{key}",
    )
    force = cols[1].toggle(i18n.t("anno.build.force"), value=False, key=f"force-{key}")

    tool_ready = env_paths.resolved_tool(aligner) is not None
    if not tool_ready:
        st.warning(i18n.t("anno.build.tool_missing", tool=aligner))
    holder = genome_writer()
    blocked = bool(holder) and holder != kind
    if blocked:
        st.caption(i18n.t("anno.build.genome_busy", kind=KIND_LABEL[holder]))
    if (assembly, source) != (current["assembly"], current["source"]):
        st.caption(
            i18n.t(
                "anno.build.other_annotation",
                id=annotation.annotation_id(assembly, source),
                current=current["id"],
            )
        )

    if st.button(
        i18n.t("anno.build.start"),
        type="primary",
        key=f"build-{key}",
        disabled=running[kind] or blocked or not tool_ready,
    ):
        progress = annotation.begin_progress(kind)
        runtime.start_job(
            annotation.build_command(kind, assembly, cores, force, progress=progress),
            f"annotation-{kind}-{assembly}",
            scope=annotation.job_scope(kind),
        )
        st.rerun()

    with st.expander(i18n.t("anno.build.sources")):
        for label, url in annotation.sources(assembly, source):
            st.markdown(f"- {label}: `{url}`")


def annotation_tab(kind: str, key: str) -> None:
    render_status(kind)
    render_build(kind, key)


def outcome(record: dict, view: dict) -> tuple[str, str]:
    """Status text and tone once the job is no longer running.

    A dead pid alone reads as "finished", which hides a script that died at the
    download step; the script's own last word, then the log, decide.
    """
    label = record.get("label", "job")
    if record.get("stopped_at"):
        return (
            i18n.t(
                "anno.job.stopped",
                label=label,
                at=i18n.stamp(record.get("stopped_at")),
            ),
            "warn",
        )
    if view.get("done"):
        return (
            i18n.t(
                "anno.job.done",
                label=label,
                elapsed=elapsed_text(record.get("started_at")),
            ),
            "ok",
        )
    tail = runtime.tail_log(record.get("log"), 80)
    if any(token in tail for token in ("[fail ]", "Traceback", "Error:")):
        return i18n.t("anno.job.failed", label=label), "fail"
    return (
        i18n.t("anno.job.ended", label=label, at=i18n.stamp(record.get("started_at"))),
        "idle",
    )


def task_block(kind: str, record: dict) -> bool:
    """One build's status, progress, buttons and log; True while it runs."""
    scope = record.get("scope") or annotation.job_scope(kind)
    state = annotation.read_progress(kind)
    view = annotation.progress_view(state)
    live = bool(record.get("running"))

    if live:
        status_line(
            i18n.t(
                "anno.job.live",
                label=record.get("label"),
                elapsed=elapsed_text(record.get("started_at")),
                pid=record.get("pid"),
            ),
            "run",
        )
    else:
        status_line(*outcome(record, view))

    if view["steps"]:
        st.progress(view["overall"], text=view["text"])
        stage_rows(view["rows"])

    if live:
        if st.button(i18n.t("anno.job.stop"), key=f"anno-stop-{kind}"):
            runtime.stop_job(scope)
            st.rerun()

    with st.expander(
        i18n.t("anno.job.log", kind=KIND_LABEL[kind]), expanded=not view["steps"]
    ):
        log_console(runtime.tail_log(record.get("log"), 300))
    return live


def annotation_logs() -> list[Path]:
    """Logs of annotation builds only.

    The run page owns the pipeline logs; a cleanup control here must not reach
    outside this page's own builds.
    """
    if not env_paths.RUN_ROOT.is_dir():
        return []
    return sorted(env_paths.RUN_ROOT.glob("*annotation*.log"))


def clear_log_panel() -> None:
    flash = st.session_state.pop("anno-logs-flash", "")
    if flash:
        st.success(flash)
    if not st.session_state.get("anno-logs-open"):
        return
    files = annotation_logs()
    if not files:
        st.caption(i18n.t("anno.logs.none"))
        if st.button(i18n.t("anno.logs.close"), type="tertiary", key="anno-logs-close"):
            st.session_state.pop("anno-logs-open", None)
            st.rerun()
        return
    total = sum(path.stat().st_size for path in files)
    st.warning(
        i18n.t("anno.logs.confirm", n=i18n.number(len(files)), size=i18n.size(total))
    )
    confirm, cancel, _rest = st.columns([1, 1, 4])
    if confirm.button(i18n.t("anno.logs.yes"), key="anno-logs-yes"):
        removed = freed = 0
        for path in files:
            try:
                freed += path.stat().st_size
                path.unlink()
            except OSError:
                continue
            removed += 1
        st.session_state.pop("anno-logs-open", None)
        tail = (
            i18n.t("anno.logs.freed_part", size=i18n.size(freed)) if freed else ""
        )
        st.session_state["anno-logs-flash"] = i18n.t(
            "anno.logs.freed", n=i18n.number(removed), tail=tail
        )
        st.rerun()
    if cancel.button(i18n.t("anno.logs.cancel"), type="tertiary", key="anno-logs-no"):
        st.session_state.pop("anno-logs-open", None)
        st.rerun()


def task_panel() -> None:
    head, action = st.columns([5, 1], vertical_alignment="center")
    with head:
        section_title(i18n.t("anno.section.tasks"))
    with action:
        if st.button(i18n.t("anno.logs.clear"), type="tertiary", key="anno-logs-clear"):
            st.session_state["anno-logs-open"] = True
    clear_log_panel()

    records = current_jobs()
    if not records:
        empty_state(i18n.t("anno.empty.tasks"))
        return
    live = False
    for kind in KINDS:
        if records.get(kind):
            live = task_block(kind, records[kind]) or live

    # The reference project refreshes every second; two seconds is enough for
    # a bar driven by 3-second index polls and keeps the page usable.
    if live and st.toggle(i18n.t("anno.auto_refresh"), value=True, key="anno-auto"):
        time.sleep(2)
        st.rerun()


tabs = st.tabs(
    [i18n.t("anno.tab.rna"), i18n.t("anno.tab.chip")]
)
with tabs[0]:
    section_title(i18n.t("anno.section.hisat2"))
    annotation_tab("rna", "rna")
with tabs[1]:
    section_title(i18n.t("anno.section.bowtie2"))
    annotation_tab("chip", "chip")

task_panel()
