from __future__ import annotations

import csv
import gzip
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SOFT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}/soft/{gse}_family.soft.gz"

RUN_COLUMNS = (
    "Run",
    "Experiment",
    "GSM",
    "Title",
    "spots",
    "bases",
    "size_MB",
    "LibraryStrategy",
    "LibrarySource",
    "LibraryLayout",
)

REQUEST_PAUSE = 0.4
MAX_RETRY = 3


class GeoError(RuntimeError):
    pass


def _user_agent() -> str:
    contact = os.environ.get("GEO_CONTACT_EMAIL", "").strip() or "apiarc@example.org"
    return f"APIARC/2.0 (mailto:{contact})"


def _fetch(url: str, data: bytes | None = None) -> bytes:
    last: Exception | None = None
    for attempt in range(MAX_RETRY):
        try:
            request = urllib.request.Request(url, data=data, headers={"User-Agent": _user_agent()})
            with urllib.request.urlopen(request, timeout=90) as response:
                return response.read()
        except Exception as error:  # noqa: BLE001
            last = error
            time.sleep(1.5 * (attempt + 1))
    raise GeoError(f"request failed: {url} ({last})")


def _eutils(path: str, params: dict, retmode: str = "xml") -> bytes:
    query = urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})
    payload = _fetch(f"{EUTILS}/{path}?{query}&retmode={retmode}")
    time.sleep(REQUEST_PAUSE)
    return payload


def normalise_gse(text: str) -> str:
    value = str(text or "").strip().upper()
    if not value:
        raise GeoError("empty accession")
    if not value.startswith("GSE"):
        value = "GSE" + value
    if not re.fullmatch(r"GSE\d+", value):
        raise GeoError(f"not a GSE accession: {text}")
    return value


def soft_text(gse: str) -> str:
    accession = normalise_gse(gse)
    bucket = f"GSE{int(accession[3:]) // 1000}nnn"
    raw = _fetch(SOFT_URL.format(bucket=bucket, gse=accession))
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw).decode("utf-8", "replace")
    return raw.decode("utf-8", "replace")


def series_title(text: str) -> str:
    match = re.search(r"!Series_title = (.+)", text)
    return match.group(1).strip() if match else ""


def series_summary(text: str) -> str:
    match = re.search(r"!Series_summary = (.+)", text)
    return match.group(1).strip() if match else ""


def bioproject(text: str) -> str | None:
    for line in text.splitlines():
        if "BioProject:" in line:
            match = re.search(r"PRJN[A-Z]?\d+", line)
            if match:
                return match.group(0)
    match = re.search(r"PRJN[A-Z]?\d+", text)
    return match.group(0) if match else None


def _esearch(db: str, term: str, retmax: int = 100000) -> list[str]:
    xml = _eutils("esearch.fcgi", {"db": db, "term": term, "retmax": retmax}).decode()
    return re.findall(r"<Id>(\d+)</Id>", xml)


def _parse_esummary(payload: dict) -> list[dict]:
    result = payload.get("result", {})
    rows = []
    for uid in result.get("uids", []):
        record = result.get(uid, {})
        exp = record.get("expxml", "")
        runs = record.get("runs", "")
        run_match = re.search(
            r'<Run acc="([^"]+)"[^>]*total_spots="([^"]+)"[^>]*total_bases="([^"]+)"', runs
        )
        if not run_match:
            continue
        title_match = re.search(r"<Title>([^<]+)</Title>", exp)
        title = title_match.group(1) if title_match else ""
        rows.append(
            {
                "Run": run_match.group(1),
                "Experiment": _tag(exp, "Experiment", attr=True),
                "GSM": title.split(":")[0].strip() if ":" in title else "",
                "Title": title,
                "spots": run_match.group(2),
                "bases": run_match.group(3),
                "size_MB": _size_mb(exp),
                "LibraryStrategy": _tag(exp, "LIBRARY_STRATEGY"),
                "LibrarySource": _tag(exp, "LIBRARY_SOURCE"),
                "LibraryLayout": _layout(exp),
            }
        )
    return rows


def _tag(text: str, name: str, attr: bool = False) -> str:
    pattern = rf'<{name} acc="([^"]+)"' if attr else rf"<{name}>([^<]+)</{name}>"
    match = re.search(pattern, text)
    return match.group(1) if match else ""


def _layout(text: str) -> str:
    match = re.search(r"<LIBRARY_LAYOUT>\s*<([A-Z]+)/?>", text)
    return match.group(1) if match else ""


def _size_mb(text: str) -> str:
    match = re.search(r'total_size="([^"]+)"', text)
    return str(round(int(match.group(1)) / 1e6, 1)) if match else ""


def runs_for_gse(gse: str) -> tuple[str, list[dict]]:
    accession = normalise_gse(gse)
    text = soft_text(accession)
    project = bioproject(text)
    if not project:
        raise GeoError(
            f"{accession} has no BioProject relation; the series may not be array-based on SRA"
        )
    ids = _esearch("sra", f"{project}[BioProject]")
    if not ids:
        ids = _esearch("sra", project)
    if not ids:
        raise GeoError(f"{project}: no SRA run found")
    rows: list[dict] = []
    for offset in range(0, len(ids), 100):
        batch = ids[offset : offset + 100]
        payload = json.loads(_eutils("esummary.fcgi", {"db": "sra", "id": ",".join(batch)}, "json"))
        rows += _parse_esummary(payload)
    rows.sort(key=lambda row: row["Run"])
    return project, rows


def is_rna(row: dict) -> bool:
    return "RNA" in (row.get("LibraryStrategy") or "").upper()


def is_chip(row: dict) -> bool:
    strategy = (row.get("LibraryStrategy") or "").upper()
    return "CHIP" in strategy


def write_run_info(path: Path | str, rows: list[dict]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(RUN_COLUMNS), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in RUN_COLUMNS})
    return path


def read_run_info(path: Path | str) -> list[dict]:
    path = Path(path)
    if not path.is_file():
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def lookup(gse: str, omics: str) -> list[dict]:
    _project, rows = runs_for_gse(gse)
    if omics == "RNAseq":
        preferred = [row for row in rows if is_rna(row)]
    else:
        preferred = [row for row in rows if is_chip(row)]
    return preferred or rows
