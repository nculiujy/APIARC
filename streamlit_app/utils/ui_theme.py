from __future__ import annotations

import html
import re
from datetime import datetime
from typing import Iterable, Sequence

import streamlit as st

from . import i18n

PAGES: list[tuple[str, str, str]] = []
NAV_GROUPS: list[tuple[str, list[str]]] = []

STATUS_KIND = {
    "ok": "ok",
    "done": "ok",
    "ready": "ok",
    "pass": "ok",
    "warn": "warn",
    "warning": "warn",
    "skipped": "warn",
    "run": "run",
    "running": "run",
    "fail": "fail",
    "error": "fail",
    "idle": "idle",
    "unknown": "idle",
}

CSS = """
<style>
:root {
  --ux-bg: #f4f6f9;
  --ux-surface: #ffffff;
  --ux-border: #e4e9ef;
  --ux-border-strong: #d3dbe4;
  --ux-text: #16202b;
  --ux-muted: #667686;
  --ux-faint: #8a99a8;
  --ux-accent: #1f6feb;
  --ux-accent-soft: #eaf1ff;
  --ux-ok: #12805c;
  --ux-ok-soft: #e6f6ef;
  --ux-warn: #a86412;
  --ux-warn-soft: #fdf3e3;
  --ux-fail: #b42318;
  --ux-fail-soft: #fdecec;
  --ux-run: #1f6feb;
  --ux-run-soft: #eaf1ff;
  --ux-radius: 12px;
  --ux-radius-sm: 9px;
  --ux-shadow: 0 1px 2px rgba(16, 32, 48, .05);
  --ux-mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, monospace;
}

.stApp { background: var(--ux-bg); color: var(--ux-text); }
header[data-testid="stHeader"] { background: transparent; }
#MainMenu, footer { visibility: hidden; }
.block-container { padding: 1.4rem 2.4rem 3.2rem; max-width: 1320px; }
.block-container [data-testid="stVerticalBlock"] { gap: .68rem; }
.block-container [data-testid="stHorizontalBlock"] { gap: 1rem; }

/* Streamlit cancels the trailing paragraph margin of every markdown block with
   a -16px bottom margin. A block that closes with a div has no such margin, so
   the negative one cut 16px into the next element. Zeroing both sides leaves
   the block gap as the only spacing. */
[data-testid="stMarkdownContainer"],
[data-testid="stCaptionContainer"] { margin-bottom: 0 !important; }
[data-testid="stMarkdownContainer"] > p,
[data-testid="stCaptionContainer"] > p { margin-bottom: 0 !important; }
/* Two paragraphs in one markdown block would touch once the margins above are gone. */
[data-testid="stMarkdownContainer"] > p + p { margin-top: .55em !important; }
/* A list keeps its own bottom margin, which the -16px above used to cancel. */
[data-testid="stMarkdownContainer"] > ul,
[data-testid="stMarkdownContainer"] > ol { margin-bottom: 0 !important; padding-left: 1.35rem; }

[data-testid="stSidebar"] {
  background: #0f1b2a;
  border-right: 1px solid #16283c;
}
[data-testid="stSidebar"] [data-testid="stVerticalBlock"] { gap: .45rem; }
[data-testid="stSidebarNav"],
[data-testid="stSidebarNavItems"],
[data-testid="stSidebarNavSeparator"] { display: none; }
[data-testid="stSidebar"] * { color: #dbe5f0; }
[data-testid="stSidebar"] .ux-brand {
  padding: 2px 2px 12px 2px;
  border-bottom: 1px solid rgba(255,255,255,.08);
}
[data-testid="stSidebar"] .ux-brand-name {
  font-size: 1.02rem; font-weight: 650; letter-spacing: .2px; color: #ffffff;
}
[data-testid="stSidebar"] .ux-brand-sub {
  font-size: .72rem; color: #8fa3b8; margin-top: 3px; letter-spacing: .3px;
}
[data-testid="stSidebar"] .ux-side-state {
  display: flex; align-items: center; justify-content: space-between;
  margin: 8px 0 0 0; padding: 5px 9px; border-radius: 8px;
  background: rgba(255,255,255,.04);
}
[data-testid="stSidebar"] .ux-side-state .k {
  font-size: .7rem; color: #8fa3b8; letter-spacing: .6px;
}
[data-testid="stSidebar"] .ux-pill.ux-pill {
  background: rgba(255,255,255,.08); color: #cbd8e5; border-color: rgba(255,255,255,.12);
}
[data-testid="stSidebar"] .ux-pill.ux-pill.run { background: rgba(31,111,235,.26); color: #cfe0ff; }
[data-testid="stSidebar"] .ux-pill.ux-pill.ok { background: rgba(18,128,92,.24); color: #8ce3bf; }
[data-testid="stSidebar"] .ux-nav-group {
  font-size: .7rem; letter-spacing: 1px;
  color: #6d8399; margin: 12px 0 3px 4px;
}
[data-testid="stSidebar"] a[data-testid="stPageLink-NavLink"] {
  border-radius: 8px; padding: 5px 10px; margin-bottom: 1px;
}
[data-testid="stSidebar"] a[data-testid="stPageLink-NavLink"] [data-testid="stIconMaterial"] {
  font-size: 17px; color: #8fa3b8;
}
[data-testid="stSidebar"] a[data-testid="stPageLink-NavLink"]:hover {
  background: rgba(255,255,255,.06);
}
[data-testid="stSidebar"] a[data-testid="stPageLink-NavLink"]:hover [data-testid="stIconMaterial"] {
  color: #cbd8e5;
}
[data-testid="stSidebar"] .ux-side-foot {
  font-size: .7rem; color: #6d8399; border-top: 1px solid rgba(255,255,255,.08);
  margin-top: 14px; padding-top: 9px; line-height: 1.6;
}

.ux-hero {
  background: var(--ux-surface); border: 1px solid var(--ux-border);
  border-radius: var(--ux-radius); box-shadow: var(--ux-shadow);
  padding: 15px 20px 16px;
}
.ux-hero h1 {
  font-size: 1.26rem; margin: 0; font-weight: 640; letter-spacing: -.2px;
}
.ux-hero p { margin: 3px 0 0; color: var(--ux-muted); font-size: .86rem; }
.ux-hero-row { display: flex; gap: 6px; flex-wrap: wrap; margin-top: 10px; }

.ux-pill {
  display: inline-flex; align-items: center; gap: 5px;
  padding: 2px 9px; border-radius: 999px;
  font-size: .74rem; font-weight: 560; letter-spacing: .1px;
  border: 1px solid transparent; white-space: normal; max-width: 100%;
}
.ux-pill .dot { width: 6px; height: 6px; border-radius: 50%; background: currentColor; }
.ux-pill.ok { background: var(--ux-ok-soft); color: var(--ux-ok); border-color: #c4e8d9; }
.ux-pill.warn { background: var(--ux-warn-soft); color: var(--ux-warn); border-color: #f2dfbe; }
.ux-pill.fail { background: var(--ux-fail-soft); color: var(--ux-fail); border-color: #f7d3d0; }
.ux-pill.run { background: var(--ux-run-soft); color: var(--ux-run); border-color: #cfe0ff; }
.ux-pill.idle { background: #f1f4f7; color: #5c6b7a; border-color: #e1e8ee; }

.ux-section {
  display: flex; align-items: baseline; gap: 8px; flex-wrap: wrap;
  margin: 12px 0 0 0; padding-bottom: 5px;
  border-bottom: 1px solid var(--ux-border);
}
.ux-section h2 { font-size: .98rem; margin: 0; font-weight: 620; }
.ux-section span { font-size: .78rem; color: var(--ux-muted); }

.ux-card {
  background: var(--ux-surface); border: 1px solid var(--ux-border);
  border-radius: var(--ux-radius); box-shadow: var(--ux-shadow);
  padding: 14px 16px;
}
.ux-card h3 { font-size: .84rem; margin: 0 0 8px 0; font-weight: 600; color: var(--ux-muted); }

.ux-grid {
  display: grid; gap: 10px;
  grid-template-columns: repeat(auto-fit, minmax(148px, 1fr));
}
.ux-metric {
  background: var(--ux-surface); border: 1px solid var(--ux-border);
  border-radius: var(--ux-radius); box-shadow: var(--ux-shadow);
  padding: 11px 14px; min-width: 0;
}
.ux-metric .k { font-size: .72rem; color: var(--ux-muted); letter-spacing: .2px; }
.ux-metric .v {
  font-size: 1.22rem; font-weight: 620; margin-top: 3px;
  font-variant-numeric: tabular-nums; letter-spacing: -.3px;
  overflow-wrap: break-word;
}
.ux-metric .v.sm { font-size: .88rem; font-weight: 560; letter-spacing: 0; line-height: 1.45; }
.ux-metric .n {
  font-size: .72rem; color: var(--ux-faint); margin-top: 3px;
  overflow-wrap: anywhere;
}

.ux-kv { display: grid; grid-template-columns: minmax(104px, 168px) 1fr; column-gap: 18px; }
.ux-kv .k, .ux-kv .v { padding: 5px 0; border-bottom: 1px solid #f0f3f7; }
.ux-kv .k { color: var(--ux-muted); font-size: .8rem; }
.ux-kv .v {
  font-family: var(--ux-mono); font-size: .78rem; color: #24313d;
  word-break: break-all;
}

/* Column glossaries are plain markdown tables: one row per value column. */
[data-testid="stMarkdownContainer"] table { border-collapse: collapse; font-size: .82rem; }
[data-testid="stMarkdownContainer"] th,
[data-testid="stMarkdownContainer"] td {
  border: none; padding: 3px 20px 3px 0; text-align: left; vertical-align: top;
}
[data-testid="stMarkdownContainer"] thead th {
  color: var(--ux-muted); font-weight: 560;
  border-bottom: 1px solid var(--ux-border); padding-bottom: 4px;
}
[data-testid="stMarkdownContainer"] tbody td:first-child {
  white-space: normal; overflow-wrap: break-word;
}

.ux-log {
  background: #0f1b2a; color: #cbd8e5; border-radius: 10px;
  padding: 12px 14px; max-height: 460px; overflow: auto;
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: .75rem; line-height: 1.55; white-space: pre-wrap;
}
.ux-log .l-err { color: #ff9b8f; }
.ux-log .l-warn { color: #ffca7a; }
.ux-log .l-ok { color: #7fe0b4; }

.ux-rows { display: flex; flex-direction: column; gap: 5px; }

.ux-row {
  display: flex; justify-content: space-between; align-items: center; gap: 12px;
  padding: 7px 12px; border: 1px solid var(--ux-border); border-radius: var(--ux-radius-sm);
  background: var(--ux-surface);
}
.ux-row > div { min-width: 0; }
.ux-row .t { font-size: .84rem; font-weight: 560; line-height: 1.4; overflow-wrap: break-word; }
.ux-row .s { font-size: .75rem; color: var(--ux-muted); line-height: 1.4; overflow-wrap: break-word; }

.ux-empty {
  border: 1px dashed var(--ux-border-strong); border-radius: var(--ux-radius);
  padding: 20px; text-align: center; color: var(--ux-faint); font-size: .84rem;
  background: var(--ux-surface);
}

.ux-stage { display: flex; flex-direction: column; gap: 5px; }
.ux-stage-row {
  display: grid; grid-template-columns: 20px 1fr auto; gap: 10px; align-items: center;
  padding: 6px 12px; border-radius: var(--ux-radius-sm); background: var(--ux-surface);
  border: 1px solid var(--ux-border);
}
.ux-stage-row .idx {
  font-size: .68rem; color: var(--ux-faint); font-variant-numeric: tabular-nums;
  text-align: right;
}
.ux-stage-row .t { font-size: .84rem; font-weight: 560; line-height: 1.4; color: var(--ux-text); }
.ux-stage-row .s {
  font-family: var(--ux-mono); font-size: .74rem; line-height: 1.4; color: var(--ux-muted);
  word-break: break-all;
}

div[data-testid="stButton"] > button,
div[data-testid="stDownloadButton"] > button {
  border-radius: var(--ux-radius-sm); font-size: .83rem; font-weight: 560;
  border: 1px solid var(--ux-border-strong); padding: .34rem .9rem; min-height: 2.05rem;
}
div[data-testid="stButton"] > button[kind="primary"] {
  background: var(--ux-accent); border-color: var(--ux-accent); color: #fff;
}
div[data-testid="stButton"] > button:hover { border-color: var(--ux-accent); color: var(--ux-accent); }
div[data-testid="stButton"] > button[kind="primary"]:hover { color: #fff; filter: brightness(1.06); }

div[data-baseweb="input"] > div, div[data-baseweb="select"] > div,
div[data-baseweb="textarea"] > div { border-radius: var(--ux-radius-sm) !important; }

[data-testid="stTabs"] [role="tablist"] {
  gap: 4px; border-bottom: 1px solid var(--ux-border); padding-bottom: 0;
}
[data-testid="stTab"] { padding: 4px 12px 6px; border-radius: 8px 8px 0 0; }
[data-testid="stTab"] p { font-size: .84rem; font-weight: 560; }
[data-testid="stTab"]:hover p { color: var(--ux-accent); }
[data-testid="stTab"][aria-selected="true"] p { color: var(--ux-accent); font-weight: 640; }
[data-testid="stTab"] .react-aria-SelectionIndicator { background: var(--ux-accent); height: 2px; }

[data-testid="stExpander"] details {
  border: 1px solid var(--ux-border); border-radius: var(--ux-radius-sm);
  background: var(--ux-surface);
}
[data-testid="stExpander"] summary { padding: 2px 6px; border-radius: var(--ux-radius-sm); }
[data-testid="stExpander"] summary:hover { background: #f6f8fb; }
[data-testid="stExpander"] summary p { font-size: .84rem; font-weight: 560; }

div[data-testid="stDataFrame"] {
  border-radius: var(--ux-radius-sm); overflow: hidden; border: 1px solid var(--ux-border);
}

div[data-testid="stAlert"] { border-radius: var(--ux-radius-sm); font-size: .84rem; }

[data-testid="stProgressBarTrack"] {
  height: 6px !important; background: #e9eef5 !important; border-radius: 999px;
}
[data-testid="stProgressBarTrack"] > div { border-radius: 999px; }
[data-testid="stProgress"] [data-testid="stMarkdownContainer"] p {
  font-size: .78rem; color: var(--ux-muted); margin: 0 0 3px 0;
}

[data-testid="stCaptionContainer"] p { font-size: .79rem; color: var(--ux-muted); }
hr { margin: 6px 0; border-color: var(--ux-border); }
</style>
"""


def inject_theme() -> None:
    st.markdown(CSS, unsafe_allow_html=True)


def _esc(text: object) -> str:
    return html.escape(str(text), quote=False)


def _attr(text: object) -> str:
    return html.escape(str(text), quote=True)


def pill(text: str, kind: str = "idle", dot: bool = True) -> str:
    css = STATUS_KIND.get(str(kind).lower(), "idle")
    marker = '<span class="dot"></span>' if dot else ""
    return f'<span class="ux-pill {css}">{marker}{_esc(text)}</span>'


def page_header(title: str, subtitle: str = "", chips: Sequence[tuple[str, str]] = ()) -> None:
    body = [f"<h1>{_esc(title)}</h1>"]
    if subtitle:
        body.append(f"<p>{_esc(subtitle)}</p>")
    if chips:
        row = "".join(pill(text, kind) for text, kind in chips)
        body.append(f'<div class="ux-hero-row">{row}</div>')
    st.markdown(f'<div class="ux-hero">{"".join(body)}</div>', unsafe_allow_html=True)


def section_title(title: str, hint: str = "") -> None:
    tail = f"<span>{_esc(hint)}</span>" if hint else ""
    st.markdown(
        f'<div class="ux-section"><h2>{_esc(title)}</h2>{tail}</div>',
        unsafe_allow_html=True,
    )


def card(title: str, body_html: str) -> None:
    head = f"<h3>{_esc(title)}</h3>" if title else ""
    st.markdown(f'<div class="ux-card">{head}{body_html}</div>', unsafe_allow_html=True)


def metrics(items: Sequence[tuple[str, str, str]]) -> None:
    cells = []
    for key, value, note in items:
        text = str(value)
        css = "v sm" if len(text) > 16 or "/" in text else "v"
        # Notes carry dataset and annotation names that run past the card; the
        # card grows a line instead of cutting the text.
        note_html = f'<div class="n">{_esc(note)}</div>' if note else ""
        # Dotted library names (TxDb.Mmusculus...) break between dots rather
        # than mid-word when the card is too narrow for one line.
        body = _esc(text).replace("/", "/&#8203;")
        if len(text) > 20:
            body = body.replace(".", ".&#8203;")
        cells.append(
            '<div class="ux-metric">'
            f'<div class="k">{_esc(key)}</div>'
            f'<div class="{css}">{body}</div>{note_html}</div>'
        )
    st.markdown(f'<div class="ux-grid">{"".join(cells)}</div>', unsafe_allow_html=True)


def key_values(items: Iterable[tuple[str, str]]) -> None:
    rows = "".join(
        f'<div class="k">{_esc(key)}</div><div class="v">{_esc(value)}</div>' for key, value in items
    )
    st.markdown(f'<div class="ux-kv">{rows}</div>', unsafe_allow_html=True)


def log_console(text: str, max_lines: int = 400) -> None:
    lines = text.splitlines()[-max_lines:]
    rendered = []
    for line in lines:
        lower = line.lower()
        css = ""
        if "error" in lower or lower.startswith("[fail"):
            css = "l-err"
        elif "warn" in lower:
            css = "l-warn"
        elif "done" in lower or lower.startswith("[ok"):
            css = "l-ok"
        rendered.append(f'<span class="{css}">{_esc(line)}</span>' if css else _esc(line))
    st.markdown(f'<div class="ux-log">{"<br>".join(rendered)}</div>', unsafe_allow_html=True)


def empty_state(text: str) -> None:
    st.markdown(f'<div class="ux-empty">{_esc(text)}</div>', unsafe_allow_html=True)


def row(title: str, subtitle: str = "", status: str = "", kind: str = "idle") -> None:
    badge = pill(status, kind) if status else ""
    sub = f'<div class="s">{_esc(subtitle)}</div>' if subtitle else ""
    st.markdown(
        f'<div class="ux-row"><div><div class="t">{_esc(title)}</div>{sub}</div>{badge}</div>',
        unsafe_allow_html=True,
    )


def list_rows(items: Sequence[tuple[str, str, str, str]]) -> None:
    """Stacked rows inside one block, so the gaps stay tight.

    Called once per row, each row became its own markdown element and picked up
    the page-wide gap, which spaced a plain list like a set of cards.
    """
    rows = []
    for title, subtitle, status, kind in items:
        badge = pill(status, kind) if status else ""
        sub = f'<div class="s">{_esc(subtitle)}</div>' if subtitle else ""
        rows.append(f'<div class="ux-row"><div><div class="t">{_esc(title)}</div>{sub}</div>{badge}</div>')
    st.markdown(f'<div class="ux-rows">{"".join(rows)}</div>', unsafe_allow_html=True)


def stage_rows(items: Sequence[tuple[str, str, str, str]], numbered: bool = True) -> None:
    rows = []
    for index, (name, detail, status, kind) in enumerate(items, start=1):
        marker = f"{index:02d}" if numbered else "•"
        rows.append(
            '<div class="ux-stage-row">'
            f'<div class="idx">{marker}</div>'
            f'<div><div class="t">{_esc(name)}</div><div class="s">{_esc(detail)}</div></div>'
            f"{pill(status, kind)}</div>"
        )
    st.markdown(f'<div class="ux-stage">{"".join(rows)}</div>', unsafe_allow_html=True)


def status_line(text: str, kind: str = "idle") -> None:
    st.markdown(pill(text, kind), unsafe_allow_html=True)


def elapsed_text(stamp: str) -> str:
    """How long a job has been running, from its ``started_at`` record."""
    try:
        started = datetime.strptime(str(stamp), "%Y-%m-%d %H:%M:%S")
    except (TypeError, ValueError):
        return "-"
    return i18n.duration((datetime.now() - started).total_seconds())


def render_sidebar(
    groups: dict[str, list], running: str = "", busy: bool = False, footer: str | None = None, active: str = ""
) -> None:
    from .version import APP_FULL_NAME

    with st.sidebar:
        # The current item is marked by its own href: st.page_link renders no
        # aria-current, and the class it uses for the active state is an emotion
        # hash that changes with every Streamlit build.
        st.markdown(
            "<style>"
            f'[data-testid="stSidebar"] a[data-testid="stPageLink-NavLink"][href="{_attr(active)}"] '
            "{ background: rgba(31,111,235,.22); box-shadow: inset 2px 0 0 var(--ux-accent); }"
            f'[data-testid="stSidebar"] a[data-testid="stPageLink-NavLink"][href="{_attr(active)}"] p '
            "{ color: #ffffff; font-weight: 620; }"
            f'[data-testid="stSidebar"] a[data-testid="stPageLink-NavLink"][href="{_attr(active)}"] '
            "[data-testid=\"stIconMaterial\"] { color: #9dc0ff; }"
            "</style>",
            unsafe_allow_html=True,
        )
        st.markdown(
            '<div class="ux-brand">'
            f'<div class="ux-brand-name">{_esc(APP_FULL_NAME)}</div>'
            f'<div class="ux-brand-sub">{_esc(i18n.t("nav.sidebar.subtitle"))}</div>'
            "</div>",
            unsafe_allow_html=True,
        )
        st.markdown(
            f'<div class="ux-side-state"><span class="k">{_esc(i18n.t("nav.sidebar.task"))}</span>'
            f'{pill(running or i18n.t("common.idle"), "run" if busy else "idle")}</div>',
            unsafe_allow_html=True,
        )
        for group, pages in groups.items():
            st.markdown(f'<div class="ux-nav-group">{_esc(group)}</div>', unsafe_allow_html=True)
            for page in pages:
                st.page_link(page, label=page.title)
        foot = footer or i18n.t("nav.sidebar.local")
        st.markdown(f'<div class="ux-side-foot">{_esc(foot)}</div>', unsafe_allow_html=True)


ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def strip_ansi(text: str) -> str:
    return ANSI.sub("", text)
