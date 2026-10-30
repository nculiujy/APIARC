from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import time
from datetime import datetime
from pathlib import Path

import yaml

from . import env_paths

STATUS_READABLE = "readable"
STATUS_ABSENT = "absent"
STATUS_BROKEN = "broken"
STATUS_PENDING = "pending"


def status_text(code: str) -> str:
    """Message for a dependency status code, resolved at render time."""
    from . import i18n

    return i18n.t(f"runtime.status.{code}")


def config_path() -> Path:
    return env_paths.CONFIG_FILE


def load_config(path: Path | str | None = None) -> dict:
    target = Path(path or config_path())
    if not target.is_file():
        return {}
    with open(target) as handle:
        data = yaml.safe_load(handle) or {}
    return data if isinstance(data, dict) else {}


def deep_merge(base: dict, patch: dict) -> dict:
    merged = dict(base)
    for key, value in patch.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def save_config(config: dict, path: Path | str | None = None) -> Path:
    target = Path(path or config_path())
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return target


def update_config(patch: dict) -> dict:
    config = deep_merge(load_config(), patch)
    save_config(config)
    return config


def config_mtime() -> str:
    target = config_path()
    if not target.is_file():
        return "unknown"
    from . import i18n

    return i18n.stamp(datetime.fromtimestamp(target.stat().st_mtime))


def runtime_marker() -> Path:
    return env_paths.RUNTIME_DIR / "runtime"


def runtime_backend() -> str:
    marker = runtime_marker()
    if marker.is_file():
        return marker.read_text().strip() or "conda"
    return "conda"


def launcher_missing() -> list[str]:
    python = env_paths.launcher_python()
    if not python.is_file():
        return ["streamlit", "pyyaml", "pandas", "altair"]
    probes = [("streamlit", "streamlit"), ("pyyaml", "yaml"), ("pandas", "pandas"), ("altair", "altair")]
    missing = []
    for package, module in probes:
        result = subprocess.run(
            [str(python), "-c", f"import {module}"],
            capture_output=True,
        )
        if result.returncode != 0:
            missing.append(package)
    return missing


def launcher_state() -> tuple[str, str]:
    from . import i18n

    prefix = env_paths.LAUNCH_PREFIX
    if not prefix.is_dir():
        return i18n.t("runtime.launcher.absent"), "fail"
    missing = launcher_missing()
    if not missing:
        return i18n.t("runtime.status.ready"), "ok"
    return i18n.t("runtime.launcher.missing", names=", ".join(missing)), "warn"


def analysis_missing() -> list[str]:
    prefix = env_paths.analysis_prefix()
    if not prefix.is_dir():
        return ["analysis-env"]
    return [name for name in ("python", "snakemake") if not (prefix / "bin" / name).is_file()]


def analysis_state() -> tuple[str, str]:
    from . import i18n

    prefix = env_paths.analysis_prefix()
    if not prefix.is_dir():
        return i18n.t("runtime.analysis.absent"), "fail"
    missing = analysis_missing()
    if not missing:
        return i18n.t("runtime.status.ready"), "ok"
    return i18n.t("runtime.analysis.missing", names=", ".join(missing)), "warn"


def tool_state(name: str) -> tuple[str, str]:
    from . import i18n

    found = env_paths.resolved_tool(name)
    text = i18n.t("runtime.status.ready" if found else "runtime.status.missing")
    return text, "ok" if found else "warn"


def env_report() -> list[dict]:
    from . import i18n

    launcher_text, launcher_kind = launcher_state()
    analysis_text, analysis_kind = analysis_state()
    rows = [
        {
            "name": i18n.t("runtime.row.launcher"),
            "path": str(env_paths.LAUNCH_PREFIX),
            "status": launcher_text,
            "kind": launcher_kind,
        },
        {
            "name": i18n.t("runtime.row.analysis"),
            "path": str(env_paths.analysis_prefix()),
            "status": analysis_text,
            "kind": analysis_kind,
        },
    ]
    for tool in ("conda", "snakemake", "prefetch", "fasterq-dump", "pigz"):
        text, kind = tool_state(tool)
        rows.append({"name": tool, "path": env_paths.resolved_tool(tool) or "-", "status": text, "kind": kind})
    return rows


def ready_for_run() -> bool:
    return env_paths.analysis_ready() and env_paths.resolved_tool("snakemake") is not None


def analysis_bin_dir() -> Path:
    return env_paths.analysis_bin()


def stage_definition_root() -> Path:
    return env_paths.BASE_DIR / "workflow" / "envs"


def base_definition() -> Path:
    return env_paths.BASE_DIR / "environment.yml"


def definition_files() -> list[tuple[str, Path]]:
    files = [("environment.yml", base_definition())]
    stages = sorted(stage_definition_root().glob("*.yaml"), key=_stage_order)
    files.extend((path.stem, path) for path in stages)
    return files


def _stage_order(path: Path) -> tuple[int, str]:
    name = path.stem
    for index, prefix in enumerate(("RNA", "ChIP", "Integrated")):
        if name.startswith(prefix):
            return (index, name)
    return (9, name)


def _prefix_packages(prefix: Path) -> dict[str, dict]:
    meta = prefix / "conda-meta"
    packages: dict[str, dict] = {}
    if not meta.is_dir():
        return packages
    for path in meta.glob("*.json"):
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        name = str(data.get("name") or "")
        if not name:
            continue
        packages[name.lower()] = {
            "version": str(data.get("version") or ""),
            "prefix": prefix,
            "files": [str(item) for item in data.get("files") or []],
        }
    return packages


def installed_packages() -> dict[str, dict]:
    roots = [env_paths.analysis_prefix()]
    stage_root = env_paths.stage_envs()
    if stage_root.is_dir():
        roots.extend(sorted(path for path in stage_root.iterdir() if (path / "conda-meta").is_dir()))
    packages: dict[str, dict] = {}
    for prefix in roots:
        for name, info in _prefix_packages(prefix).items():
            packages.setdefault(name, info)
    return packages


def stage_env_count() -> int:
    root = env_paths.stage_envs()
    if not root.is_dir():
        return 0
    return sum(1 for path in root.iterdir() if (path / "conda-meta").is_dir())


def declared_dependencies() -> dict[str, list[dict]]:
    declared: dict[str, list[dict]] = {}
    for label, path in definition_files():
        if not path.is_file():
            continue
        try:
            data = yaml.safe_load(path.read_text()) or {}
        except (OSError, yaml.YAMLError):
            continue
        rows = declared.setdefault(label, [])
        for spec in data.get("dependencies") or []:
            if isinstance(spec, dict):
                continue
            text = str(spec).strip()
            if not text:
                continue
            name = text.split("=")[0].strip()
            pinned = text[len(name):].lstrip("=")
            rows.append({"spec": text, "package": name, "pinned": pinned})
    return declared


def _command_entry(files: list[str]) -> str:
    for entry in files:
        if entry.startswith("bin/") and "/" not in entry[4:] and not entry.endswith(".py"):
            return entry
    return ""


def dependency_report() -> dict[str, list[dict]]:
    packages = installed_packages()
    report: dict[str, list[dict]] = {}
    for label, rows in declared_dependencies().items():
        scope = "base" if label == "environment.yml" else "stage"
        for row in rows:
            info = packages.get(row["package"].lower())
            if info is None:
                if scope == "base":
                    status, kind = STATUS_ABSENT, "fail"
                else:
                    status, kind = STATUS_PENDING, "idle"
                report.setdefault(label, []).append(
                    {
                        **row,
                        "version": "-",
                        "command": "-",
                        "prefix": "-",
                        "status": status,
                        "kind": kind,
                        "scope": scope,
                    }
                )
                continue
            entry = _command_entry(info["files"])
            binary = info["prefix"] / entry if entry else None
            if entry and not (binary and os.access(binary, os.X_OK)):
                status, kind = STATUS_BROKEN, "fail"
            else:
                status, kind = STATUS_READABLE, "ok"
            report.setdefault(label, []).append(
                {
                    **row,
                    "version": info["version"],
                    "command": Path(entry).name if entry else "-",
                    "prefix": env_paths.relative(info["prefix"]),
                    "status": status,
                    "kind": kind,
                    "scope": scope,
                }
            )
    return report


def dependency_summary(report: dict[str, list[dict]] | None = None) -> dict:
    report = report if report is not None else dependency_report()
    base = [row for rows in report.values() for row in rows if row["scope"] == "base"]
    stage = [row for rows in report.values() for row in rows if row["scope"] == "stage"]
    return {
        "base_total": len(base),
        "base_ready": sum(1 for row in base if row["status"] == STATUS_READABLE),
        "base_absent": sum(1 for row in base if row["status"] != STATUS_READABLE),
        "stage_total": len(stage),
        "stage_ready": sum(1 for row in stage if row["status"] == STATUS_READABLE),
        "stage_pending": sum(1 for row in stage if row["status"] != STATUS_READABLE),
        "broken": sum(1 for row in base if row["status"] == STATUS_BROKEN),
    }


def panel_report() -> list[dict]:
    import importlib.util
    from importlib import metadata

    probes = (
        ("streamlit", "streamlit", "streamlit"),
        ("pandas", "pandas", "pandas"),
        ("numpy", "numpy", "numpy"),
        ("pyyaml", "yaml", "PyYAML"),
        ("altair", "altair", "altair"),
        ("pyBigWig", "pyBigWig", "pyBigWig"),
        ("pymupdf", "fitz", "PyMuPDF"),
    )
    rows = []
    for label, module, distribution in probes:
        spec = importlib.util.find_spec(module)
        version = "-"
        if spec:
            try:
                version = metadata.version(distribution)
            except metadata.PackageNotFoundError:
                version = "-"
        origin = spec.origin if spec and spec.origin else "-"
        rows.append(
            {
                "name": label,
                "version": version,
                "status": STATUS_READABLE if spec else STATUS_ABSENT,
                "kind": "ok" if spec else "warn",
                "path": origin,
                "where": _short_origin(origin),
            }
        )
    return rows


def panel_site_dir() -> Path:
    import importlib.util

    spec = importlib.util.find_spec("streamlit")
    origin = spec.origin if spec and spec.origin else ""
    marker = "/site-packages/"
    if marker in origin:
        return Path(origin.split(marker, 1)[0] + marker.rstrip("/"))
    return env_paths.LAUNCH_PREFIX


def _short_origin(origin: str) -> str:
    marker = "/site-packages/"
    if marker in origin:
        return origin.split(marker, 1)[1]
    return origin


def _jobs_path() -> Path:
    return env_paths.RUN_STATE / "jobs.json"


def read_jobs() -> dict[str, dict]:
    path = _jobs_path()
    if path.is_file():
        try:
            payload = json.loads(path.read_text())
        except (OSError, ValueError):
            payload = {}
        if isinstance(payload, dict):
            return payload
    legacy = env_paths.RUN_STATE / "job.json"
    if legacy.is_file():
        try:
            state = json.loads(legacy.read_text())
        except (OSError, ValueError):
            return {}
        if isinstance(state, dict) and state:
            return {"default": state}
    return {}


def write_jobs(jobs: dict[str, dict]) -> None:
    env_paths.RUN_STATE.mkdir(parents=True, exist_ok=True)
    _jobs_path().write_text(json.dumps(jobs, indent=2))
    legacy = env_paths.RUN_STATE / "job.json"
    if legacy.exists():
        legacy.unlink()


def _prefs_path() -> Path:
    return env_paths.RUN_STATE / "ui.json"


def read_pref(name: str, default: str = "") -> str:
    path = _prefs_path()
    if not path.is_file():
        return default
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError):
        return default
    return str(payload.get(name, default))


def save_pref(name: str, value: str) -> None:
    path = _prefs_path()
    payload: dict = {}
    if path.is_file():
        try:
            payload = json.loads(path.read_text())
        except (OSError, ValueError):
            payload = {}
    payload[name] = value
    env_paths.RUN_STATE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def _alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return True
    fields = stat.rsplit(")", 1)[-1].split()
    return bool(fields) and fields[0] != "Z"


def job_status(scope: str | None = None) -> dict:
    jobs = read_jobs()
    for key, state in jobs.items():
        state["scope"] = key
        state["running"] = _alive(int(state.get("pid") or 0))
    if scope is not None:
        return dict(jobs.get(scope) or {})
    if not jobs:
        return {}
    running = [state for state in jobs.values() if state["running"]]
    pool = running or list(jobs.values())
    return dict(max(pool, key=lambda state: str(state.get("started_at", ""))))


def running_jobs() -> list[dict]:
    jobs = read_jobs()
    active = []
    for key, state in jobs.items():
        if _alive(int(state.get("pid") or 0)):
            state["scope"] = key
            state["running"] = True
            active.append(state)
    return active


def is_running(scope: str | None = None) -> bool:
    if scope is None:
        return bool(running_jobs())
    return bool(job_status(scope).get("running"))


def current_user() -> dict:
    uid = os.geteuid() if hasattr(os, "geteuid") else -1
    name = ""
    try:
        import pwd

        name = pwd.getpwuid(uid).pw_name
    except (ImportError, KeyError):
        name = ""
    name = name or os.environ.get("USER") or os.environ.get("LOGNAME") or "-"
    return {"name": name, "uid": uid, "admin": uid == 0}


def new_log_path(label: str) -> Path:
    env_paths.RUN_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in label)[:40]
    return env_paths.RUN_ROOT / f"{stamp}_{safe}.log"


def start_job(
    command: list[str], label: str, scope: str = "default", log_path: Path | None = None
) -> dict:
    from . import i18n

    if is_running(scope):
        raise RuntimeError(i18n.t("runtime.job.busy", scope=scope))
    env_paths.ensure_dirs()
    log_path = Path(log_path or new_log_path(label))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handle = open(log_path, "a")
    handle.write(f"# {datetime.now():%Y-%m-%d %H:%M:%S} {label}\n")
    handle.write("$ " + " ".join(command) + "\n\n")
    handle.flush()
    process = subprocess.Popen(
        command,
        cwd=str(env_paths.BASE_DIR),
        env=env_paths.routed_env(),
        stdout=handle,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    jobs = read_jobs()
    jobs[scope] = {
        "pid": process.pid,
        "label": label,
        "command": command,
        "log": str(log_path),
        "started_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "cwd": str(env_paths.BASE_DIR),
    }
    write_jobs(jobs)
    return job_status(scope)


def _signal_group(pid: int, sig: int) -> None:
    try:
        os.killpg(os.getpgid(pid), sig)
    except OSError:
        try:
            os.kill(pid, sig)
        except OSError:
            pass


def _terminate(pid: int, timeout: float) -> None:
    for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGKILL):
        _signal_group(pid, sig)
        deadline = time.time() + timeout / 3
        while time.time() < deadline:
            if not _alive(pid):
                return
            time.sleep(0.2)


def stop_job(scope: str | None = None, timeout: float = 12.0) -> bool:
    jobs = read_jobs()
    keys = [scope] if scope is not None else list(jobs)
    stopped = False
    for key in keys:
        state = jobs.get(key)
        if state is None:
            continue
        pid = int(state.get("pid") or 0)
        if _alive(pid):
            _terminate(pid, timeout)
        state["running"] = False
        state["stopped_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        stopped = True
    write_jobs(jobs)
    return stopped


def clear_state(scope: str | None = None) -> None:
    jobs = read_jobs()
    if scope is None:
        jobs = {}
    else:
        jobs.pop(scope, None)
    write_jobs(jobs)


def tail_log(path: Path | str | None, lines: int = 200) -> str:
    """Last lines of a log, read from the end.

    Rule logs reach hundreds of kilobytes and the log drawers read one on every
    rerun, so loading the whole file each time is wasted work.
    """
    if not path:
        return ""
    target = Path(path)
    if not target.is_file():
        return ""
    block = max(64 * 1024, int(lines) * 400)
    try:
        size = target.stat().st_size
        with open(target, "rb") as handle:
            handle.seek(max(0, size - block))
            chunk = handle.read()
    except OSError:
        return ""
    text = chunk.decode("utf-8", errors="replace")
    if size > block:
        # The seek lands mid-line; drop what is left of it.
        text = text.split("\n", 1)[-1]
    content = text.splitlines()
    if len(content) < lines and size > block:
        try:
            content = target.read_text(errors="replace").splitlines()
        except OSError:
            pass
    return "\n".join(content[-lines:])


RUN_LOG_PATTERN = re.compile(r"^(\d{8})-(\d{6})_(?P<label>.+)$")


def run_records(scopes: tuple[str, ...] | None = None) -> list[dict]:
    """Runs kept in logs/runs, newest first.

    Each file is named after its start time and its label, so the run history
    survives restarts; jobs.json only remembers the latest run per slot.
    """
    root = env_paths.RUN_ROOT
    if not root.is_dir():
        return []
    records = []
    for path in root.glob("*.log"):
        match = RUN_LOG_PATTERN.match(path.stem)
        if not match:
            continue
        try:
            started = datetime.strptime(match.group(1) + match.group(2), "%Y%m%d%H%M%S")
        except ValueError:
            continue
        label = match.group("label")
        scope = label.split("-", 1)[0]
        if scopes is not None and scope not in scopes:
            continue
        try:
            size = path.stat().st_size
        except OSError:
            size = 0
        records.append(
            {
                "path": path,
                "label": label,
                "scope": scope,
                "started": started,
                "stamp": started.strftime("%Y-%m-%d %H:%M:%S"),
                "size": size,
            }
        )
    records.sort(key=lambda item: item["started"])
    # Runs are numbered per omics, so "#3" means the third RNA-seq run rather
    # than a position in a mixed list.
    counters: dict[str, int] = {}
    for record in records:
        counters[record["scope"]] = counters.get(record["scope"], 0) + 1
        record["run"] = counters[record["scope"]]
    # Rule logs keep only their latest content, so a log belongs to the run
    # whose window it was last written in. The window closes where the next run
    # of the same scope begins.
    for position, record in enumerate(records):
        later = [item for item in records[position + 1:] if item["scope"] == record["scope"]]
        record["until"] = later[0]["started"] if later else datetime.now()
    return list(reversed(records))


def module_logs(scope: str, start: datetime, until: datetime) -> list[Path]:
    """Rule logs written between two runs, in pipeline order."""
    root = env_paths.LOG_ROOT / scope
    if not root.is_dir():
        return []
    picked = []
    for path in root.glob("*.log"):
        try:
            stamp = datetime.fromtimestamp(path.stat().st_mtime)
        except OSError:
            continue
        if start <= stamp <= until:
            picked.append(path)
    picked.sort(key=_log_order)
    return picked


def _log_order(path: Path) -> tuple[list[int], str]:
    numbers: list[int] = []
    for part in path.name.split("_")[:2]:
        if not part.isdigit():
            break
        numbers.append(int(part))
    return (numbers, path.name)


def _log_files() -> list[Path]:
    """Logs the run drawers list, minus the ones a live job still holds.

    Unlinking the log of a running job only detaches the name: the process
    keeps writing to the removed inode, so the drawer goes empty while the run
    continues.
    """
    busy = {Path(str(state["log"])) for state in running_jobs() if state.get("log")}
    # A set, because logs/runs lives under logs and would otherwise be counted
    # once through RUN_ROOT and once through the folder walk.
    files: set[Path] = set()
    if env_paths.RUN_ROOT.is_dir():
        files.update(env_paths.RUN_ROOT.glob("*.log"))
    if env_paths.LOG_ROOT.is_dir():
        for folder in env_paths.LOG_ROOT.iterdir():
            if folder.is_dir() and folder != env_paths.RUN_ROOT:
                files.update(folder.rglob("*.log"))
    return sorted(path for path in files if path.is_file() and path not in busy)


def log_footprint() -> tuple[int, int]:
    """File count and total size of the logs a cleanup would remove."""
    files = _log_files()
    total = 0
    for path in files:
        try:
            total += path.stat().st_size
        except OSError:
            continue
    return len(files), total


def clear_logs() -> tuple[int, int]:
    """Remove every log the drawers show; returns count and bytes freed."""
    removed = freed = 0
    for path in _log_files():
        try:
            size = path.stat().st_size
            path.unlink()
        except OSError:
            continue
        removed += 1
        freed += size
    return removed, freed


def recent_logs(limit: int = 12) -> list[Path]:
    root = env_paths.RUN_ROOT
    if not root.is_dir():
        return []
    files = [path for path in root.iterdir() if path.is_file() and path.suffix == ".log"]
    files.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    return files[:limit]


def result_size(path: Path) -> str:
    if not path.exists():
        return "-"
    try:
        output = subprocess.run(
            ["du", "-sh", str(path)], capture_output=True, text=True, check=False
        )
        return output.stdout.split()[0] if output.stdout else "-"
    except OSError:
        return "-"
