"""Locale catalogues, the lookup used by the pages, and the formatters beside it.

Every sentence a user can read lives in ``streamlit_app/i18n/<locale>.json`` and
is fetched by a dotted key, so no display text is compiled into the pages. The
formatters follow the same catalogues: English groups thousands and writes dates
as ``Oct 29, 2026``, Chinese leaves digits ungrouped and writes ``2026-10-29``.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import streamlit as st

from . import runtime

LOCALES = ("en", "zh")
DEFAULT_LOCALE = "en"

_CATALOGUE_DIR = Path(__file__).resolve().parents[1] / "i18n"
_CATALOGUES: dict[str, dict[str, str]] = {}
_NAMES: dict[str, str] = {}
_GROUPED = {"en": True, "zh": False}
_DATE = {"en": "%b %d, %Y", "zh": "%Y-%m-%d"}
_DATETIME = {"en": "%b %d, %Y %H:%M", "zh": "%Y-%m-%d %H:%M"}
_UNITS = ("B", "KB", "MB", "GB", "TB")


def catalogue(name: str) -> dict[str, str]:
    if name not in _CATALOGUES:
        path = _CATALOGUE_DIR / f"{name}.json"
        with open(path, encoding="utf-8") as handle:
            raw = json.load(handle)
        _NAMES[name] = str(raw.get("_name") or name)
        _CATALOGUES[name] = {
            key: value for key, value in raw.items() if not key.startswith("_")
        }
    return _CATALOGUES[name]


def native_name(name: str) -> str:
    """The locale's own name, as its catalogue spells it."""
    catalogue(name)
    return _NAMES.get(name, name)


def current() -> str:
    """The active locale, seeded from ``ui.locale`` in config.yaml."""
    try:
        stored = st.session_state.get("locale")
    except Exception:  # noqa: BLE001 - no script context, e.g. a bare import
        stored = None
    if stored in LOCALES:
        return stored
    value = str((runtime.load_config().get("ui") or {}).get("locale") or "")
    value = value if value in LOCALES else DEFAULT_LOCALE
    try:
        st.session_state["locale"] = value
    except Exception:  # noqa: BLE001
        pass
    return value


def set_locale(value: str) -> None:
    if value not in LOCALES:
        return
    st.session_state["locale"] = value
    section = dict(runtime.load_config().get("ui") or {})
    if section.get("locale") != value:
        section["locale"] = value
        runtime.update_config({"ui": section})


def t(key: str, **fields: object) -> str:
    """Catalogue text for ``key``, falling back to English and then to the key."""
    text = catalogue(current()).get(key)
    if text is None:
        text = catalogue(DEFAULT_LOCALE).get(key, key)
    if not fields:
        return text
    try:
        return text.format(**fields)
    except (KeyError, IndexError, ValueError):
        return text


def number(value: object, digits: int = 0) -> str:
    try:
        amount = float(value)
    except (TypeError, ValueError):
        return str(value)
    group = "," if _GROUPED[current()] else ""
    return f"{amount:{group}.{digits}f}"


def size(value: object) -> str:
    try:
        amount = float(value)
    except (TypeError, ValueError):
        return str(value)
    for unit in _UNITS:
        if amount < 1024 or unit == _UNITS[-1]:
            return f"{number(amount, 0 if unit == 'B' else 1)} {unit}"
        amount /= 1024
    return str(value)


def stamp(value: object, with_time: bool = True) -> str:
    moment: datetime | None = None
    if isinstance(value, datetime):
        moment = value
    elif isinstance(value, (int, float)):
        moment = datetime.fromtimestamp(value)
    elif value:
        for pattern in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
            try:
                moment = datetime.strptime(str(value), pattern)
                break
            except ValueError:
                continue
    if moment is None:
        return "-"
    return moment.strftime(_DATETIME[current()] if with_time else _DATE[current()])


def duration(seconds: object) -> str:
    try:
        total = max(0, int(float(seconds)))
    except (TypeError, ValueError):
        return "-"
    hours, rest = divmod(total, 3600)
    minutes, secs = divmod(rest, 60)
    if hours:
        return t("common.duration.hours", hours=hours, minutes=f"{minutes:02d}")
    if minutes:
        return t("common.duration.minutes", minutes=minutes, seconds=f"{secs:02d}")
    return t("common.duration.seconds", seconds=secs)


def percent(value: object, digits: int = 0) -> str:
    try:
        amount = float(value) * 100
    except (TypeError, ValueError):
        return "-"
    return f"{amount:.{digits}f}%"


def ratio(done: object, total: object) -> str:
    return f"{number(done)} / {number(total)}"


def language_switch(key: str = "locale-switch") -> None:
    """Two-state picker; the label is hidden because the native names speak."""
    options = list(LOCALES)
    active = current()
    label = t("common.language")
    if hasattr(st, "segmented_control"):
        picked = st.segmented_control(
            label,
            options,
            default=active,
            format_func=native_name,
            key=key,
            label_visibility="collapsed",
        )
    else:
        picked = st.radio(
            label,
            options,
            index=options.index(active),
            format_func=native_name,
            horizontal=True,
            key=key,
            label_visibility="collapsed",
        )
    if picked and picked != active:
        set_locale(picked)
        st.rerun()
