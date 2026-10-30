from __future__ import annotations

import os
import sys
from pathlib import Path


def _env_path(name: str, default: Path) -> Path:
    value = os.environ.get(name, "").strip()
    return Path(os.path.expanduser(value)) if value else default


def app_dir() -> Path:
    return Path(__file__).resolve().parents[1]


def repo_root() -> Path:
    return app_dir().parent


BASE_DIR = repo_root()
APP_DIR = app_dir()

RUNTIME_DIR = _env_path("APIARC_RUNTIME_DIR", BASE_DIR / ".apiarc-runtime")
LAUNCH_PREFIX = _env_path("APIARC_LAUNCH_ENV", RUNTIME_DIR / "launch")
ANALYSIS_PREFIX = _env_path("APIARC_ANALYSIS_ENV", BASE_DIR / "analysis-env")

CONFIG_DIR = BASE_DIR / "config"
CONFIG_FILE = CONFIG_DIR / "config.yaml"
TEST_CONFIG_FILE = CONFIG_DIR / "test_config.yaml"

DATA_ROOT = BASE_DIR / "workflow" / "data"
RESULT_ROOT = BASE_DIR / "result"
LOG_ROOT = BASE_DIR / "logs"

RUN_ROOT = LOG_ROOT / "runs"
RUN_STATE = RUNTIME_DIR / "state"

SCRIPT_ROOT = BASE_DIR / "workflow" / "scripts"


def launcher_python() -> Path:
    return LAUNCH_PREFIX / "bin" / "python"


def launcher_streamlit() -> Path:
    return LAUNCH_PREFIX / "bin" / "streamlit"


def analysis_marker() -> Path:
    return RUN_STATE / "analysis_env"


def analysis_prefix() -> Path:
    marker = analysis_marker()
    try:
        if marker.is_file():
            text = marker.read_text(encoding="utf-8").strip()
            if text:
                return Path(os.path.expanduser(text))
    except OSError:
        pass
    return ANALYSIS_PREFIX


def analysis_bin() -> Path:
    return analysis_prefix() / "bin"


def stage_envs() -> Path:
    return analysis_prefix() / "envs"


def analysis_ready() -> bool:
    return (analysis_bin() / "snakemake").exists()


def launcher_ready() -> bool:
    return launcher_streamlit().exists() and launcher_python().exists()


def analysis_env_marker_write(prefix: Path | str | None = None) -> Path:
    target = Path(os.path.expanduser(str(prefix or ANALYSIS_PREFIX))).resolve()
    RUN_STATE.mkdir(parents=True, exist_ok=True)
    analysis_marker().write_text(str(target) + "\n", encoding="utf-8")
    return target


def tool_path() -> str:
    parts: list[str] = []
    for candidate in (analysis_bin(), Path(sys.prefix) / "bin"):
        text = str(candidate)
        if candidate.is_dir() and text not in parts:
            parts.append(text)
    for item in os.environ.get("PATH", "").split(os.pathsep):
        if item and item not in parts:
            parts.append(item)
    return os.pathsep.join(parts)


def routed_env(base: dict[str, str] | None = None) -> dict[str, str]:
    env = dict(base if base is not None else os.environ)
    env["PATH"] = tool_path()
    env["APIARC_REPO"] = str(BASE_DIR)
    env["APIARC_RUNTIME_DIR"] = str(RUNTIME_DIR)
    env["APIARC_ANALYSIS_ENV"] = str(analysis_prefix())
    env["APIARC_LAUNCH_ENV"] = str(LAUNCH_PREFIX)
    env["CONDA_PKGS_DIRS"] = str(RUNTIME_DIR / "pkgs")
    env["CONDA_ENVS_DIRS"] = str(RUNTIME_DIR / "envs")
    env["PIP_CACHE_DIR"] = str(RUNTIME_DIR / "cache" / "pip")
    # A shell that already has conda base activated exports CONDA_SHLVL, but the
    # launcher never puts that env on the child's PATH. conda activate then
    # rewrites the inherited stack instead of prepending the rule environment,
    # so snakemake's conda jobs silently ran under the launcher interpreter.
    for leftover in (
        "CONDA_SHLVL",
        "CONDA_PREFIX",
        "CONDA_DEFAULT_ENV",
        "CONDA_PYTHON_EXE",
        "CONDA_PROMPT_MODIFIER",
    ):
        env.pop(leftover, None)
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    return env


def resolved_tool(name: str) -> str | None:
    import shutil

    return shutil.which(name, path=tool_path())


def relative(path: Path | str) -> str:
    path = Path(path)
    try:
        return str(path.resolve().relative_to(BASE_DIR))
    except (ValueError, OSError):
        return str(path)


def ensure_dirs() -> None:
    for path in (RUN_ROOT, RUN_STATE, DATA_ROOT / "RNAseq", DATA_ROOT / "ChIPseq"):
        path.mkdir(parents=True, exist_ok=True)
