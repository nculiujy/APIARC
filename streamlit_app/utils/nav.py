from __future__ import annotations

import os

import streamlit as st

from . import i18n

APP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_GROUPS = (
    (
        "nav.group.start",
        (
            ("overview", "nav.overview", "views/0_overview.py", ":material/dashboard:"),
            ("env", "nav.environment", "views/1_environment.py", ":material/tune:"),
            ("annotation", "nav.annotation", "views/2_annotation.py", ":material/dataset:"),
            ("prepare", "nav.prepare", "views/3_prepare.py", ":material/upload_file:"),
        ),
    ),
    (
        "nav.group.run",
        (
            ("run", "nav.run", "views/4_run.py", ":material/play_circle:"),
            ("joint", "nav.joint", "views/5_joint.py", ":material/hub:"),
        ),
    ),
    (
        "nav.group.review",
        (("results", "nav.results", "views/6_results.py", ":material/insights:"),),
    ),
)


def page_groups() -> dict[str, list[st.Page]]:
    groups: dict[str, list[st.Page]] = {}
    for group_key, items in _GROUPS:
        pages = []
        for key, title_key, relative, icon in items:
            pages.append(
                st.Page(
                    os.path.join(APP_DIR, relative),
                    title=i18n.t(title_key),
                    icon=icon,
                    url_path=key,
                    default=key == "overview",
                )
            )
        groups[i18n.t(group_key)] = pages
    return groups


def activate(groups: dict[str, list[st.Page]]) -> st.Page:
    return st.navigation(groups, position="hidden")
