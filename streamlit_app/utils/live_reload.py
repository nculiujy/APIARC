"""Forget cached modules when their sources change.

Streamlit reruns the entry script on every rerun, but it keeps the modules it
already imported in ``sys.modules``. Its own watcher clears the module entry and
leaves the package attribute of ``utils`` pointing at the old object, so
``from utils import datasets`` kept returning the previous module: adding an
``annotation`` argument to ``datasets.result_directory`` still raised
``TypeError`` on a page that had merely been open. ``workflow/scripts/common``
is not watched at all. Dropping the packages themselves, and their submodules,
lets a browser refresh pick either change up.
"""

from __future__ import annotations

import sys
from pathlib import Path

from . import env_paths

_SHARED_PACKAGES = ("common", "utils")
_SELF = __name__
_STAMP: dict[str, tuple[float, int]] = {}


def _snapshot(root: Path) -> dict[str, tuple[float, int]]:
    stamp: dict[str, tuple[float, int]] = {}
    for path in root.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        try:
            info = path.stat()
        except OSError:
            continue
        stamp[str(path)] = (info.st_mtime, info.st_size)
    return stamp


def _bases() -> list[Path]:
    """Both halves of the code that a page run depends on."""
    return [env_paths.app_dir(), env_paths.SCRIPT_ROOT]


def refresh(
    root: Path | str | None = None, packages: tuple[str, ...] = _SHARED_PACKAGES
) -> list[str]:
    """Drop cached modules of ``packages`` if any watched source changed.

    ``root`` narrows the watch to one directory; the default covers the app and
    the shared ``workflow/scripts`` package. Returns the changed files, so a
    caller can log what triggered the reload.
    """
    global _STAMP
    bases = [Path(root)] if root is not None else _bases()
    current: dict[str, tuple[float, int]] = {}
    for base in bases:
        current.update(_snapshot(base))
    previous, _STAMP = _STAMP, current
    if not previous:
        return []
    changed = sorted(path for path, stamp in current.items() if previous.get(path) != stamp)
    if not changed:
        return []
    for name in list(sys.modules):
        if name == _SELF:
            continue
        # The package module is dropped too: ``from utils import nav`` reads the
        # attribute off the package, which would otherwise still be the old module.
        if name in packages or name.startswith(tuple(pkg + "." for pkg in packages)):
            sys.modules.pop(name, None)
    return changed
