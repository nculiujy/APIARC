from __future__ import annotations

import os
import socket
import sys

import streamlit as st

_APP_DIR = os.path.dirname(os.path.abspath(__file__))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import live_reload  # noqa: E402

# A server that stays up keeps the modules it imported at start, and the shared
# pipeline code lives outside the app directory: without this, a changed
# signature there only took effect after a restart.
live_reload.refresh()

from utils import i18n, nav, runtime  # noqa: E402
from utils.ui_theme import inject_theme, render_sidebar  # noqa: E402
from utils.version import APP_FULL_NAME  # noqa: E402

st.set_page_config(page_title=APP_FULL_NAME, layout="wide", initial_sidebar_state="expanded")
inject_theme()

running = [job.get("label", "pipeline") for job in runtime.running_jobs()]
running_label = i18n.t("nav.sidebar.running", names=" / ".join(running)) if running else ""
try:
    footer = f"{os.environ.get('USER') or 'user'} @ {socket.gethostname()}"
except Exception:  # noqa: BLE001
    footer = None

groups = nav.page_groups()
page = nav.activate(groups)
render_sidebar(
    groups, running=running_label, busy=bool(running), footer=footer, active=page.url_path or ""
)
page.run()
