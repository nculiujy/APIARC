from __future__ import annotations

import os
import sys

import pandas as pd
import streamlit as st

_APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import datasets, env_paths, i18n, prepare, runtime  # noqa: E402
from utils.ui_theme import (  # noqa: E402
    elapsed_text,
    empty_state,
    log_console,
    metrics,
    page_header,
    pill,
    section_title,
    stage_rows,
    status_line,
)

running = bool(runtime.running_jobs())
config = runtime.load_config()


def read_prep(omics: str) -> dict:
    stored = st.session_state.get(f"prep-{omics}")
    if stored:
        return dict(stored)
    # A browser refresh opens a new Streamlit session and session_state is
    # empty; the dataset being edited sits in ui.json, its table on disk.
    name = runtime.read_pref(f"prepare-{omics}")
    draft = prepare.read_draft(omics, name) if name else {}
    samples = list(draft.get("chosen") or [])
    if not samples:
        return {}
    prep = {"name": name, "samples": samples, "source": draft.get("source") or "local"}
    st.session_state[f"prep-{omics}"] = prep
    return dict(prep)


def forget_rows(omics: str) -> None:
    for key in (f"prep-rows-sig-{omics}", f"prep-rows-{omics}"):
        st.session_state.pop(key, None)


def write_prep(omics: str, name: str, samples: list[str], source: str) -> None:
    st.session_state[f"prep-{omics}"] = {
        "name": name,
        "samples": list(samples),
        "source": source,
    }
    forget_rows(omics)
    runtime.save_pref(f"prepare-{omics}", name)
    # A new pick at step 1 invalidates the rows drafted for the previous one.
    prepare.write_draft(omics, name, {"source": source, "chosen": list(samples), "rows": []})


# What each column of the naming table is for. The table is the only place the
# columns are read from, so the wording stays next to it as a plain list.
COLUMN_HINT_KEY = {
    "RNAseq": (
        ("sample", "prep.hint.sample_fastq"),
        ("sample_name", "prep.hint.rna_name"),
        ("group", "prep.hint.rna_group"),
    ),
    "ChIPseq": (
        ("sample", "prep.hint.sample_fastq"),
        ("role", "prep.hint.role"),
        ("sample_name", "prep.hint.chip_name"),
        ("condition", "prep.hint.condition"),
        ("input_sample", "prep.hint.input_sample"),
    ),
}


def column_hints(omics: str) -> None:
    lines = "\n".join(
        f"| `{name}` | {i18n.t(key)} |" for name, key in COLUMN_HINT_KEY[omics]
    )
    st.markdown(
        i18n.t("prep.hint.table", header=i18n.t("prep.hint.column"), meaning=i18n.t("prep.hint.meaning"), rows=lines)
    )


def clean_rows(frame: pd.DataFrame) -> list[dict]:
    rows = []
    for record in frame.to_dict("records"):
        row = {}
        for key, value in record.items():
            if value is None or (isinstance(value, float) and pd.isna(value)):
                row[key] = ""
            else:
                # Trailing spaces would be stripped again when the pipeline
                # reads the CSV, splitting one group into two.
                row[key] = str(value).strip()
        rows.append(row)
    return prepare.drop_blank_rows(rows)


page_header(
    i18n.t("prep.title"),
    chips=[
        (
            i18n.t("prep.chip.rna", name=datasets.config_dataset(config, "RNAseq")),
            "ok",
        ),
        (
            i18n.t("prep.chip.chip", name=datasets.config_dataset(config, "ChIPseq")),
            "ok",
        ),
        (
            i18n.t(
                "anno.chip.job",
                state=i18n.t("anno.job.running") if running else i18n.t("anno.job.idle"),
            ),
            "run" if running else "idle",
        ),
    ],
)


def render_current(omics: str) -> None:
    names = datasets.available(omics)
    current = datasets.config_dataset(config, omics)
    if not names:
        empty_state(i18n.t("prep.empty.datasets", omics=omics))
        return
    index = names.index(current) if current in names else 0
    cols = st.columns([2, 1], vertical_alignment="bottom")
    picked = cols[0].selectbox(
        i18n.t("prep.label.current"), names, index=index, key=f"pick-{omics}"
    )
    if cols[1].button(
        i18n.t("prep.button.use"), key=f"use-{omics}", disabled=picked == current
    ):
        runtime.update_config({datasets.config_section(config, omics): {"dataset": picked}})
        st.rerun()
    info = prepare.status(omics, picked)
    runs = datasets.srr_list(omics, picked)
    # A run listed at step 1 but absent on disk is a download still to do, not
    # a misspelt sample; the wording has to tell the two apart.
    pending = [item for item in info["missing"] if item in runs]
    cards = [
        (
            i18n.t("prep.metric.samples"),
            i18n.number(info["samples"]),
            i18n.t("prep.metric.samples_hint") if omics == "ChIPseq" else "",
        ),
        (i18n.t("prep.metric.pairs"), i18n.number(info["pairs"]), ""),
        (i18n.t("prep.metric.singles"), i18n.number(info["singles"]), ""),
        (
            i18n.t("prep.metric.missing"),
            i18n.number(len(info["missing"])),
            i18n.t("prep.metric.missing_hint", n=i18n.number(len(pending)))
            if pending
            else "",
        ),
    ]
    if omics == "ChIPseq":
        cards.insert(
            1,
            (
                i18n.t("prep.metric.inputs"),
                i18n.number(info["inputs"]),
                i18n.t("prep.metric.inputs_hint"),
            ),
        )
    if info["ready"]:
        verdict = i18n.t("prep.state.ready")
    elif pending:
        verdict = i18n.t("prep.state.pending", n=i18n.number(len(pending)))
    elif info["missing"]:
        verdict = i18n.t("prep.state.missing", n=i18n.number(len(info["missing"])))
    else:
        verdict = i18n.t("prep.state.incomplete")
    with st.expander(f"{picked} · {verdict}"):
        metrics(cards)


def render_steps(omics: str) -> None:
    prep = read_prep(omics)
    target = prep.get("name", "")
    samples = len(prep.get("samples") or [])
    done = bool(target) and prepare.prepared(omics, target)
    # Only a GEO dataset has raw data still to fetch; a directory imported from
    # disk already holds its FASTQ and would sit on "not downloaded" forever.
    waiting = bool(target) and bool(datasets.srr_list(omics, target))
    loaded = waiting and prepare.status(omics, target)["ready"]
    hint = i18n.t("prep.step2.hint_rna" if omics == "RNAseq" else "prep.step2.hint_chip")
    rows = [
        (
            i18n.t("prep.step1.title"),
            (
                i18n.t("prep.step1.detail", name=target, n=i18n.number(samples))
                if prep
                else i18n.t("prep.step1.empty", omics=omics)
            ),
            i18n.t("prep.step1.ok") if prep else i18n.t("prep.step1.wait"),
            "ok" if prep else "idle",
        ),
        (
            i18n.t("prep.step2.title"),
            i18n.t("prep.step2.done", name=target) if done else hint,
            i18n.t("prep.step2.ok") if done else i18n.t("prep.step2.wait"),
            "ok" if done else "idle",
        ),
    ]
    if waiting:
        rows.append(
            (
                i18n.t("prep.row.raw"),
                (
                    i18n.t("prep.row.raw_ready")
                    if loaded
                    else i18n.t(
                        "prep.row.raw_wait",
                        n=i18n.number(len(datasets.srr_list(omics, target))),
                    )
                ),
                i18n.t("prep.row.ok") if loaded else i18n.t("prep.row.wait"),
                "ok" if loaded else "idle",
            )
        )
    stage_rows(rows)


def save_geo(omics: str, name: str, payload: dict, rows: list[dict], chosen: list[str]) -> str | None:
    try:
        clean = datasets.validate_name(name)
    except ValueError as error:
        st.error(str(error))
        return None
    path = prepare.register_geo(omics, clean, payload["gse"], rows, chosen)
    write_prep(omics, clean, chosen, "geo")
    st.success(
        i18n.t("prep.geo.saved", path=os.path.relpath(path, env_paths.BASE_DIR))
    )
    return clean


def display_title(value: str) -> str:
    text = str(value or "").strip()
    tag = ""
    if "[" in text and "]" in text:
        start = text.index("[")
        end = text.index("]", start)
        tag = text[start + 1 : end].strip()
        text = (text[:start] + " " + text[end + 1 :]).strip()
    if text.upper().startswith("GSM") and ":" in text:
        text = text.split(":", 1)[1].strip()
    text = " ".join(text.split()).replace(" ;", ";").lstrip(";").strip()
    return f"[{tag}] {text}" if tag else text


def render_geo(omics: str) -> None:
    # The dataset is named after the accession; a second box for the same value
    # only invited one of the two to be edited.
    cols = st.columns([3, 1], vertical_alignment="bottom")
    gse = cols[0].text_input(
        i18n.t("prep.geo.label"), key=f"gse-{omics}", placeholder="GSE140552"
    )
    if cols[1].button(
        i18n.t("prep.geo.fetch"), key=f"fetch-{omics}", disabled=not gse.strip()
    ):
        try:
            with st.spinner(i18n.t("prep.geo.reading")):
                rows = prepare.fetch_geo(gse, omics)
        except Exception as error:  # noqa: BLE001
            st.error(i18n.t("prep.geo.failed", error=error))
        else:
            st.session_state[f"geo-{omics}"] = {"gse": gse.strip().upper(), "rows": rows}

    payload = st.session_state.get(f"geo-{omics}") or {}
    rows = payload.get("rows") or []
    if not rows:
        return
    edited = st.data_editor(
        pd.DataFrame(rows).assign(title=lambda frame: frame["title"].map(display_title)).drop(
            columns=["strategy", "spots"], errors="ignore"
        ),
        hide_index=True,
        width="stretch",
        key=f"geo-table-{omics}",
        disabled=["sample", "gsm", "layout", "size_mb", "title"],
        column_config={
            "use": st.column_config.CheckboxColumn(i18n.t("prep.geo.column.use")),
            "sample": st.column_config.TextColumn(i18n.t("prep.geo.column.sample")),
            "gsm": st.column_config.TextColumn("GSM"),
            "layout": st.column_config.TextColumn(i18n.t("prep.geo.column.layout")),
            "size_mb": st.column_config.TextColumn("MB"),
            "title": st.column_config.TextColumn(i18n.t("prep.geo.column.title")),
        },
    )
    chosen = edited.loc[edited["use"], "sample"].tolist()
    st.caption(
        i18n.t("prep.geo.selected", n=i18n.number(len(chosen)), total=i18n.number(len(rows)))
    )

    cols = st.columns([1, 3], vertical_alignment="bottom")
    if cols[0].button(
        i18n.t("prep.geo.save"), type="primary", key=f"save-geo-{omics}", disabled=not chosen
    ):
        if save_geo(omics, payload["gse"], payload, rows, chosen):
            st.rerun()


def render_local(omics: str) -> None:
    names = [name for name in datasets.available(omics) if datasets.directory(omics, name).is_dir()]
    if not names:
        empty_state(i18n.t("prep.local.empty", omics=omics))
        return
    name = st.selectbox(i18n.t("prep.local.label"), names, key=f"local-{omics}")
    rows = prepare.scan_samples(omics, name)
    if not rows:
        empty_state(i18n.t("prep.local.none"))
        return
    edited = st.data_editor(
        pd.DataFrame(rows),
        hide_index=True,
        width="stretch",
        key=f"local-table-{omics}",
        disabled=["sample", "read1", "read2", "layout", "size_gb"],
        column_config={
            "use": st.column_config.CheckboxColumn(i18n.t("prep.geo.column.use")),
            "sample": st.column_config.TextColumn(i18n.t("prep.geo.column.sample")),
            "read1": st.column_config.TextColumn("R1"),
            "read2": st.column_config.TextColumn("R2"),
            "layout": st.column_config.TextColumn(i18n.t("prep.geo.column.layout")),
            "size_gb": st.column_config.NumberColumn("GB", format="%.2f"),
        },
    )
    chosen = edited.loc[edited["use"], "sample"].tolist()
    st.caption(
        i18n.t("prep.geo.selected", n=i18n.number(len(chosen)), total=i18n.number(len(rows)))
    )
    if st.button(
        i18n.t("prep.local.use"), type="primary", key=f"use-local-{omics}", disabled=not chosen
    ):
        write_prep(omics, name, chosen, "local")
        st.rerun()


def render_source(omics: str) -> None:
    options = ["geo", "local"]
    picked = st.radio(
        i18n.t("prep.source.label"),
        options,
        format_func=lambda value: i18n.t(f"prep.source.{value}"),
        key=f"mode-{omics}",
        horizontal=True,
        label_visibility="collapsed",
    )
    if picked == "geo":
        render_geo(omics)
    else:
        render_local(omics)


def editor_columns(omics: str, options: list[str]) -> dict:
    if omics == "ChIPseq":
        # Both columns offer the same names: whether a library is the IP or the
        # control is an experiment design question, not something its name says.
        if options:
            sample_col = st.column_config.SelectboxColumn(
                i18n.t("prep.column.sample"),
                options=["", *options],
                help=i18n.t("prep.help.sample_options"),
            )
            input_col = st.column_config.SelectboxColumn(
                i18n.t("prep.column.input_sample"),
                options=["", *options],
                help=i18n.t("prep.help.input_sample"),
            )
        else:
            sample_col = st.column_config.TextColumn(
                i18n.t("prep.column.sample"), help=i18n.t("prep.help.sample_prefix")
            )
            input_col = st.column_config.TextColumn(
                i18n.t("prep.column.input_sample"),
                help=i18n.t("prep.help.input_sample_short"),
            )
        return {
            "sample": sample_col,
            "role": st.column_config.SelectboxColumn(
                i18n.t("prep.column.role"),
                options=["IP", "Input"],
                help=i18n.t("prep.help.role"),
            ),
            "sample_name": st.column_config.TextColumn(
                i18n.t("prep.column.sample_name"), help=i18n.t("prep.help.group_name")
            ),
            "condition": st.column_config.SelectboxColumn(
                i18n.t("prep.column.condition"),
                options=["", *datasets.CHIP_REPLICATES],
                help=i18n.t("prep.help.condition"),
            ),
            "input_sample": input_col,
        }
    return {
        "sample": st.column_config.TextColumn(
            i18n.t("prep.column.sample"), help=i18n.t("prep.help.sample_prefix")
        ),
        "sample_name": st.column_config.TextColumn(
            i18n.t("prep.column.sample_name"), help=i18n.t("prep.help.rna_column")
        ),
        "group": st.column_config.SelectboxColumn(
            i18n.t("prep.column.group"),
            options=list(datasets.RNA_GROUPS),
            help=i18n.t("prep.help.rna_group"),
        ),
    }


def draft(omics: str, name: str, samples: list[str]) -> tuple[list[dict], int]:
    signature = "|".join([name, *samples])
    if st.session_state.get(f"prep-rows-sig-{omics}") != signature:
        st.session_state[f"prep-rows-sig-{omics}"] = signature
        stored = prepare.read_draft(omics, name)
        kept = None
        if list(stored.get("chosen") or []) == list(samples):
            # A draft written before the role column existed still holds an
            # unnamed row per Input; normalising restores its role instead of
            # asking the user to retype the whole table.
            kept = prepare.normalize_rows(omics, stored.get("rows") or []) or None
        st.session_state[f"prep-rows-{omics}"] = kept or prepare.draft_rows(omics, name, samples)
        st.session_state[f"prep-rows-ver-{omics}"] = int(st.session_state.get(f"prep-rows-ver-{omics}", 0)) + 1
    return st.session_state[f"prep-rows-{omics}"], int(st.session_state.get(f"prep-rows-ver-{omics}", 0))


def remember(omics: str, name: str, samples: list[str], source: str, rows: list[dict]) -> None:
    stored = prepare.read_draft(omics, name)
    if list(stored.get("chosen") or []) == list(samples) and stored.get("rows") == rows:
        return
    prepare.write_draft(omics, name, {"source": source, "chosen": list(samples), "rows": rows})


def render_comparisons(omics: str, name: str, rows: list[dict]) -> None:
    """Write ``comparisons.csv``: which sample_name is compared against which.

    The table names groups but not the biological contrast; two groups alone do
    not say which one is the reference. Every downstream differential step runs
    once per pair listed here.
    """
    if omics != "ChIPseq":
        return
    section_title(i18n.t("prep.comp.title"), i18n.t("prep.comp.hint"))
    counts = {
        key: value for key, value in prepare.group_summary(omics, rows) if key
    }
    usable = sorted(key for key, count in counts.items() if count >= 2)
    if not usable:
        empty_state(i18n.t("prep.comp.empty"))
        return
    columns = list(datasets.ws.COMPARISON_COLUMNS)
    stored = datasets.comparison_rows(omics, name)
    frame = pd.DataFrame(stored or [dict.fromkeys(columns, "")], columns=columns)
    # A saved pair may cite a group the edited table no longer has enough rows
    # for; keeping it in the options shows what is stored rather than silently
    # clearing the cell on the next render.
    written = {str(value).strip() for value in frame.to_numpy().ravel() if str(value).strip()}
    options = sorted(set(usable) | written)
    edited = st.data_editor(
        frame,
        hide_index=True,
        width="stretch",
        num_rows="dynamic",
        key=f"prep-cmp-{omics}-{name}",
        column_config={
            "reference": st.column_config.SelectboxColumn(
                i18n.t("prep.column.reference"),
                options=["", *options],
                help=i18n.t("prep.help.reference"),
            ),
            "treatment": st.column_config.SelectboxColumn(
                i18n.t("prep.column.treatment"),
                options=["", *options],
                help=i18n.t("prep.help.treatment"),
            ),
        },
    )
    pairs = clean_rows(edited)
    problems: list[str] = []
    seen: set[tuple[str, str]] = set()
    for pair in pairs:
        reference = pair.get("reference", "")
        treatment = pair.get("treatment", "")
        if not reference or not treatment:
            problems.append(i18n.t("prep.comp.error.row"))
        elif reference == treatment:
            problems.append(i18n.t("prep.comp.error.same", name=reference))
        elif (reference, treatment) in seen:
            problems.append(
                i18n.t("prep.comp.error.duplicate", source=reference, target=treatment)
            )
        elif counts.get(reference, 0) < 2 or counts.get(treatment, 0) < 2:
            problems.append(
                i18n.t("prep.comp.error.replicates", source=reference, target=treatment)
            )
        seen.add((reference, treatment))
    for problem in dict.fromkeys(problems):
        st.warning(problem)

    cols = st.columns([1, 3], vertical_alignment="bottom")
    if cols[0].button(
        i18n.t("prep.comp.save"),
        key=f"prep-cmp-save-{omics}",
        disabled=not pairs or bool(problems),
    ):
        path = datasets.save_comparisons(omics, name, pairs)
        st.session_state[f"prep-flash-{omics}"] = i18n.t(
            "prep.comp.saved",
            path=os.path.relpath(path, env_paths.BASE_DIR),
            n=i18n.number(len(pairs)),
        )
        st.rerun()
    cols[1].caption(i18n.t("prep.comp.caption"))


def render_naming(omics: str) -> None:
    prep = read_prep(omics)
    if not prep:
        empty_state(i18n.t("prep.naming.empty"))
        return
    flash = st.session_state.pop(f"prep-flash-{omics}", "")
    if flash:
        st.info(flash)
    name = prep["name"]
    samples = prep["samples"]
    rows, version = draft(omics, name, samples)
    columns = list(datasets.columns(omics))
    options = prepare.sample_options(omics, name, samples, rows) if omics == "ChIPseq" else []

    column_hints(omics)

    # Counts belong beside the rules they check, but they can only be computed
    # from the edited rows; the placeholder is filled once the table returns.
    summary = st.empty()

    frame = pd.DataFrame(rows, columns=columns)
    # A column whose cells are all blank comes back as float, and the text and
    # selectbox editors refuse to edit a float column. The table is text.
    frame = frame.fillna("").astype(str)

    edited = st.data_editor(
        frame,
        hide_index=True,
        width="stretch",
        num_rows="dynamic",
        key=f"prep-editor-{omics}-{version}",
        column_config=editor_columns(omics, options),
        # The RNA sample is the step-1 pick; renaming it here would break the
        # link to the FASTQ files and the expression matrix columns.
        disabled=("sample",) if omics == "RNAseq" else (),
    )
    current_rows = prepare.normalize_rows(omics, clean_rows(edited))

    remember(omics, name, samples, prep["source"], current_rows)

    counts = prepare.group_summary(omics, current_rows)
    gaps = prepare.missing_fields(omics, current_rows)
    notes = prepare.group_notes(omics, current_rows)
    heading = " · ".join(
        i18n.t(
            "prep.naming.count",
            group=key or i18n.t("prep.group.unset"),
            n=i18n.number(value),
        )
        for key, value in counts
    ) or i18n.t("prep.naming.no_groups")
    with summary.container():
        st.markdown(
            pill(heading, "ok" if counts and not gaps and not notes else "warn"),
            unsafe_allow_html=True,
        )
        if gaps:
            head = "; ".join(gaps[:6])
            more = (
                i18n.t("prep.naming.gaps_more", n=i18n.number(len(gaps) - 6))
                if len(gaps) > 6
                else ""
            )
            st.caption(i18n.t("prep.naming.gaps", head=head, more=more))
        for note in notes:
            st.caption(note)

    if st.button(i18n.t("prep.naming.save"), type="primary", key=f"commit-{omics}"):
        errors, warnings, path = prepare.complete(
            omics, name, current_rows, prep["source"], chosen=prep["samples"]
        )
        for warning in warnings:
            st.warning(warning)
        if errors:
            for error in errors:
                st.error(error)
        else:
            runtime.update_config({datasets.config_section(config, omics): {"dataset": name}})
            st.session_state[f"flash-{omics}"] = i18n.t(
                "prep.naming.saved",
                path=os.path.relpath(path, env_paths.BASE_DIR),
                name=name,
            )
            st.rerun()

    render_comparisons(omics, name, current_rows)


def download_outcome(job: dict, view: dict) -> tuple[str, str]:
    """Status text and tone once the download is no longer running.

    A dead pid alone reads as "finished", which hides a run that died on a
    failed accession; the progress file has the script's own last word.
    """
    label = job.get("label", "download")
    if job.get("stopped_at"):
        return (
            i18n.t(
                "prep.dl.outcome.stopped",
                label=label,
                at=i18n.stamp(job.get("stopped_at")),
            ),
            "warn",
        )
    if view.get("error"):
        return i18n.t("prep.dl.outcome.failed", label=label), "fail"
    if view.get("done"):
        stamp = i18n.stamp(view.get("updated_at")) if view.get("updated_at") else ""
        tail = i18n.t("prep.dl.outcome.done_tail", at=stamp) if stamp else ""
        return i18n.t("prep.dl.outcome.done", label=label, tail=tail), "ok"
    return (
        i18n.t("prep.dl.outcome.ended", label=label, at=i18n.stamp(job.get("started_at"))),
        "idle",
    )


@st.fragment(run_every=2)
def download_status(omics: str, scope: str) -> None:
    """Status, bar, per-sample table and log for one omics' download.

    Only this block is refreshed on a timer: a full-page rerun every two
    seconds would keep resetting the naming table the user is typing in.
    """
    job = prepare.download_job(omics)
    if not job:
        empty_state(i18n.t("prep.dl.empty"))
        return
    live = bool(job.get("running"))
    view = prepare.download_view(prepare.read_download_progress(omics))
    if not view["rows"]:
        # A run launched before the progress file existed: the FASTQ on disk is
        # the only record, and counting it still beats a dead bar.
        name = read_prep(omics).get("name") or datasets.config_dataset(config, omics)
        runs = datasets.srr_list(omics, name)
        if runs:
            view = prepare.file_view(omics, name, runs)
    if live:
        status_line(
            i18n.t(
                "anno.job.live",
                label=job.get("label"),
                elapsed=elapsed_text(job.get("started_at")),
                pid=job.get("pid"),
            ),
            "run",
        )
    else:
        status_line(*download_outcome(job, view))

    if view["rows"]:
        st.progress(view["overall"], text=view["text"])
        stage_rows(view["rows"])
        if view["runs"]:
            st.dataframe(
                pd.DataFrame(view["runs"]),
                hide_index=True,
                width="stretch",
            )

    cols = st.columns([1, 3], vertical_alignment="bottom")
    if live:
        if cols[0].button(i18n.t("prep.dl.stop"), key=f"dl-stop-{omics}"):
            runtime.stop_job(job.get("scope") or scope)
            st.rerun(scope="app")
    elif cols[0].button(i18n.t("prep.dl.clear"), key=f"dl-clear-{omics}"):
        runtime.clear_state(job.get("scope") or scope)
        st.rerun(scope="app")

    with st.expander(i18n.t("prep.dl.log"), expanded=not view["rows"]):
        log_console(runtime.tail_log(job.get("log"), 200))


def render_download(omics: str) -> None:
    section_title(i18n.t("prep.dl.section"))
    prep = read_prep(omics)
    name = prep.get("name") or datasets.config_dataset(config, omics)
    runs = datasets.srr_list(omics, name)
    if not runs:
        empty_state(i18n.t("prep.dl.no_runs", name=name, omics=omics))
        return
    if not prep:
        st.caption(i18n.t("prep.dl.current", name=name))
    scope = prepare.download_scope(omics)
    job = prepare.download_job(omics)
    live = bool(job.get("running"))
    have, awaiting = prepare.downloaded(omics, name, runs)
    size = prepare.expected_bytes(omics, name, runs)
    with st.expander(
        i18n.t(
            "prep.dl.expander",
            name=name,
            have=i18n.number(have),
            total=i18n.number(len(runs)),
        )
    ):
        metrics(
            [
                (
                    i18n.t("prep.dl.metric.dataset"),
                    name,
                    i18n.t("prep.dl.metric.dataset_hint", n=i18n.number(len(runs))),
                ),
                (
                    i18n.t("prep.dl.metric.have"),
                    i18n.ratio(have, len(runs)),
                    "",
                ),
                (
                    i18n.t("prep.dl.metric.waiting"),
                    i18n.number(len(awaiting)),
                    "",
                ),
                (
                    i18n.t("prep.dl.metric.size"),
                    i18n.size(size) if size else i18n.t("common.none"),
                    i18n.t("prep.dl.metric.size_hint"),
                ),
            ]
        )
        if omics == "ChIPseq":
            ip, controls, untabled = prepare.download_roles(omics, name, runs)
            st.caption(
                i18n.t(
                    "prep.dl.roles",
                    ip=i18n.number(ip),
                    controls=i18n.number(controls),
                    untabled=i18n.number(untabled),
                )
            )

    cols = st.columns([1, 3], vertical_alignment="bottom")
    if cols[0].button(
        i18n.t("prep.dl.start"),
        type="primary",
        key=f"dl-{omics}",
        disabled=live or not awaiting,
    ):
        progress = prepare.begin_download(omics)
        try:
            runtime.start_job(
                prepare.download_command(omics, name, runs, progress),
                f"download-{omics}-{name}",
                scope=scope,
            )
        except RuntimeError as error:
            st.error(str(error))
        else:
            st.rerun()
    if not awaiting:
        cols[1].caption(i18n.t("prep.dl.all_present"))
    else:
        cols[1].caption(i18n.t("prep.dl.skip"))

    download_status(omics, scope)


def render_omics(omics: str) -> None:
    flash = st.session_state.pop(f"flash-{omics}", "")
    if flash:
        st.success(flash)
    render_current(omics)
    render_steps(omics)
    section_title(i18n.t("prep.section.step1"))
    render_source(omics)
    render_download(omics)
    section_title(i18n.t("prep.section.step2"))
    render_naming(omics)


rna_tab, chip_tab = st.tabs(["RNA-seq", "ChIP-seq"])
with rna_tab:
    render_omics("RNAseq")
with chip_tab:
    render_omics("ChIPseq")
