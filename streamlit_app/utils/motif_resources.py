"""JASPAR CORE matrices that the motif stage scans, and the files beside them.

The motif rules read three things from ``workflow/resources/Motif_tf_anno``:
one MEME file per matrix under ``JASPAR_CORE/``, the matrix-id to TF-symbol
table, and a per-species TF list that decides which matrices are worth scanning.
They are committed to the repository, so a bootstrap never fetches them; this
module exists so a missing or outdated collection can be restored from JASPAR
without leaving the dashboard.

The official download is a zip of one MEME file per matrix, which is exactly the
layout the rules expect, so installing it is unpacking plus rebuilding the id
table from the ``MOTIF <id> <name>`` lines.
"""

from __future__ import annotations

import shutil
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

from . import env_paths, i18n

RESOURCE_DIR = Path("workflow") / "resources" / "Motif_tf_anno"
MATRIX_DIR = "JASPAR_CORE"
MAP_NAME = "JASPAR_CORE_ID_to_NAME.txt"
SOURCE_NAME = "JASPAR_SOURCE.txt"
# One row per gene; the Family column ("zf-C2H2", "TF_bZIP") is AnimalTFDB's.
TF_LISTS = {"homo": "Homo_sapiens_TF.txt", "mm": "Mus_musculus_TF.txt"}

RELEASE = {
    "key": "2024",
    "label": "JASPAR CORE 2024",
    "url": (
        "https://jaspar.elixir.no/download/data/2024/CORE/"
        "JASPAR2024_CORE_non-redundant_pfms_meme.zip"
    ),
    "matrices": 2346,
}

CHUNK = 1 << 16


class DownloadError(RuntimeError):
    """Raised with a sentence the page can print as is."""


def resource_dir() -> Path:
    return env_paths.BASE_DIR / RESOURCE_DIR


def matrix_dir() -> Path:
    return resource_dir() / MATRIX_DIR


def read_source() -> dict[str, str]:
    path = resource_dir() / SOURCE_NAME
    if not path.is_file():
        return {}
    fields: dict[str, str] = {}
    with open(path, errors="replace") as handle:
        for line in handle:
            key, _, value = line.rstrip("\n").partition("\t")
            if key and value:
                fields[key] = value
    return fields


def tf_list_rows() -> dict[str, int]:
    """Gene rows per species TF list, zero when the file is absent."""
    counts: dict[str, int] = {}
    for species, name in TF_LISTS.items():
        path = resource_dir() / name
        if not path.is_file():
            counts[species] = 0
            continue
        with open(path, errors="replace") as handle:
            counts[species] = max(0, sum(1 for _ in handle) - 1)
    return counts


def status() -> dict:
    """What is on disk now, without reading the matrices themselves."""
    folder = matrix_dir()
    files = sorted(folder.glob("*.meme")) if folder.is_dir() else []
    total = sum(path.stat().st_size for path in files)
    mapping = resource_dir() / MAP_NAME
    map_rows = 0
    if mapping.is_file():
        with open(mapping, errors="replace") as handle:
            map_rows = sum(1 for line in handle if "\t" in line)
    newest = max((path.stat().st_mtime for path in files), default=0.0)
    source = read_source()
    return {
        "dir": folder,
        "exists": folder.is_dir(),
        "matrices": len(files),
        "bytes": total,
        "map_rows": map_rows,
        "tf": tf_list_rows(),
        "ready": bool(files) and map_rows > 0,
        "stamp": time.strftime("%Y-%m-%d %H:%M", time.localtime(newest)) if newest else "",
        "source": source.get("url", ""),
        "release": source.get("release", ""),
        "downloaded": source.get("downloaded", ""),
    }


def _motif_line(path: Path) -> str:
    """``<id>\\t<tf name>`` from the first MOTIF line of one matrix file."""
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith("MOTIF"):
                fields = line.rstrip("\n").split(None, 2)
                if len(fields) < 2:
                    break
                name = fields[2] if len(fields) > 2 else ""
                return f"{path.stem}\t{name.strip()}"
    return ""


def _stream(url: str, target: Path, progress=None) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "APIARC"})
    expected = done = 0
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            expected = int(response.headers.get("Content-Length") or 0)
            with open(target, "wb") as handle:
                while True:
                    block = response.read(CHUNK)
                    if not block:
                        break
                    handle.write(block)
                    done += len(block)
                    if progress and expected:
                        progress(min(1.0, done / expected), done, expected)
    except urllib.error.URLError as error:
        reason = getattr(error, "reason", error)
        if isinstance(reason, TimeoutError):
            raise DownloadError(i18n.t("motif.error.timeout")) from error
        raise DownloadError(i18n.t("motif.error.download", reason=reason)) from error
    except ValueError as error:
        raise DownloadError(i18n.t("motif.error.url", reason=error)) from error
    except OSError as error:
        raise DownloadError(i18n.t("motif.error.interrupted", reason=error)) from error
    if expected and done != expected:
        raise DownloadError(
            i18n.t(
                "motif.error.partial",
                done=i18n.number(done),
                expected=i18n.number(expected),
            )
        )


def _fetch(url: str, target: Path, progress=None, attempts: int = 3) -> None:
    """Stream the archive, retrying the connection a couple of times.

    The JASPAR host drops TLS handshakes now and then, so one retry turns a
    failed click into a slow one.
    """
    failure: DownloadError | None = None
    for attempt in range(attempts):
        try:
            _stream(url, target, progress)
            return
        except DownloadError as error:
            failure = error
            if attempt + 1 < attempts:
                time.sleep(2 + 3 * attempt)
                if progress:
                    progress(0.0, 0, 0)
    raise failure if failure else DownloadError(i18n.t("motif.error.failed"))


def _unpack(archive: Path, staging: Path) -> list[Path]:
    """Every ``*.meme`` member of the zip, keeping only the file name."""
    written: list[Path] = []
    with zipfile.ZipFile(archive) as bundle:
        for member in bundle.namelist():
            if not member.endswith(".meme"):
                continue
            target = staging / Path(member).name
            with bundle.open(member) as source, open(target, "wb") as sink:
                shutil.copyfileobj(source, sink)
            written.append(target)
    return written


def install(progress=None, target: Path | None = None) -> dict:
    """Download the release and rebuild the matrix folder and the id table.

    ``progress`` is called with (fraction, bytes_done, bytes_expected) while the
    archive streams. Files are written to a staging directory first so a failed
    download cannot leave a half-filled matrix folder behind.
    """
    folder = Path(target) if target is not None else matrix_dir()
    parent = folder.parent
    parent.mkdir(parents=True, exist_ok=True)
    staging = parent / f"{MATRIX_DIR}.download"
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir()
    try:
        with tempfile.TemporaryDirectory() as scratch:
            archive = Path(scratch) / "jaspar.zip"
            _fetch(RELEASE["url"], archive, progress)
            if progress:
                progress(1.0, 0, 0)
            written = _unpack(archive, staging)
            if not written:
                raise DownloadError(i18n.t("motif.error.no_member"))
            rows = [line for line in (_motif_line(path) for path in sorted(written)) if line]
            if len(rows) < len(written):
                raise DownloadError(
                    i18n.t("motif.error.broken", n=len(written) - len(rows))
                )
            mapping = parent / MAP_NAME
            (staging / MAP_NAME).write_text("\n".join(rows) + "\n")

            # Swap in one step: the old folder is kept until the new one is in
            # place, so an interrupted run leaves either the previous collection
            # or none, never a folder holding half of each.
            retired = parent / f"{MATRIX_DIR}.old"
            shutil.rmtree(retired, ignore_errors=True)
            if folder.exists():
                folder.rename(retired)
            staging.rename(folder)
            shutil.rmtree(retired, ignore_errors=True)
            (folder / MAP_NAME).replace(mapping)
    except BaseException:
        # The resource folder is tracked by git, so a failed run must not leave a
        # stray JASPAR_CORE.download/ that later looks like a real collection.
        shutil.rmtree(staging, ignore_errors=True)
        raise

    total = sum(path.stat().st_size for path in folder.glob("*.meme"))
    reference = folder.parent / SOURCE_NAME
    reference.write_text(
        f"release\t{RELEASE['label']}\n"
        f"url\t{RELEASE['url']}\n"
        f"matrices\t{len(rows)}\n"
        f"downloaded\t{time.strftime('%Y-%m-%d %H:%M')}\n"
    )
    return {"matrices": len(rows), "bytes": total, "map": mapping}
