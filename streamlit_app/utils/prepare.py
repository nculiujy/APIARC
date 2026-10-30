from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import yaml

from . import datasets, env_paths, i18n, runtime

DRAFT_NAME = "prepare_state.yaml"

# One download step per stage, weighted by how long it actually takes: the
# network wait owns the wall clock, not the conversion.
DOWNLOAD_KEY = {
    "download": "prep.step.download",
    "convert": "prep.step.convert",
    "compress": "prep.step.compress",
}
RUN_STAGE_KEY = {
    "queued": "prep.run.queued",
    "download": "prep.run.download",
    "convert": "prep.run.convert",
    "compress": "prep.run.compress",
    "done": "prep.run.done",
    "skip": "prep.run.skip",
    "fail": "prep.run.fail",
}


def fetch_geo(gse: str, omics: str) -> list[dict]:
    rows = []
    for row in datasets.geo_runs(gse, omics):
        rows.append(
            {
                "use": True,
                "sample": row.get("Run", ""),
                "gsm": row.get("GSM", ""),
                "strategy": row.get("LibraryStrategy", ""),
                "layout": row.get("LibraryLayout", ""),
                "spots": row.get("spots", ""),
                "size_mb": row.get("size_MB", ""),
                "title": row.get("Title", ""),
            }
        )
    return rows


def geo_run_rows(rows: list[dict], chosen: list[str]) -> list[dict]:
    picked = []
    for row in rows:
        if row.get("sample") not in chosen:
            continue
        picked.append(
            {
                "Run": row.get("sample", ""),
                "GSM": row.get("gsm", ""),
                "Title": row.get("title", ""),
                "spots": row.get("spots", ""),
                "size_MB": row.get("size_mb", ""),
                "LibraryStrategy": row.get("strategy", ""),
                "LibraryLayout": row.get("layout", ""),
            }
        )
    return picked


def register_geo(omics: str, name: str, gse: str, rows: list[dict], chosen: list[str]) -> Path:
    path = datasets.directory(omics, name)
    path.mkdir(parents=True, exist_ok=True)
    datasets.write_run_info(omics, name, geo_run_rows(rows, chosen))
    datasets.write_srr_list(omics, name, chosen)
    descriptor = dict(datasets.descriptors(omics, name))
    descriptor.update(
        {
            "source": "geo",
            "geo_accession": gse,
            "samples": len(chosen),
            "prepared": False,
        }
    )
    return datasets.save_descriptor(omics, name, descriptor)


def scan_samples(omics: str, name: str) -> list[dict]:
    rows = []
    for prefix, slot in sorted(datasets.fastq_table(omics, name).items()):
        r1 = slot.get("r1")
        r2 = slot.get("r2")
        size = sum(item.stat().st_size for item in (r1, r2) if item is not None)
        rows.append(
            {
                "use": True,
                "sample": prefix,
                "read1": r1.name if r1 else "",
                "read2": r2.name if r2 else "",
                "layout": i18n.t("prep.layout.paired" if r2 else "prep.layout.single"),
                "size_gb": round(size / 1e9, 2),
            }
        )
    return rows


def sample_options(omics: str, name: str, samples: list[str], rows: list[dict]) -> list[str]:
    """Names offered by the ChIP ``sample`` and ``input_sample`` columns.

    Candidates are the FASTQ present in the dataset directory, the step-1
    selection (a GEO dataset is described before its download starts) and the
    values already written in either column. Both columns share one pool: which
    library is the IP and which is its control is an experiment design
    question, not something the name can tell.
    """
    found: dict[str, None] = {}
    for prefix in sorted(datasets.fastq_table(omics, name)):
        found.setdefault(prefix, None)
    for item in (
        samples
        + [row.get("sample") for row in rows]
        + [row.get("input_sample") or row.get("input_id") for row in rows]
    ):
        text = str(item or "").strip()
        if text:
            found.setdefault(text, None)
    # Three sources, each already in some order; concatenated they interleave.
    # A dropdown is searched visually, so one sorted list reads better.
    return sorted(found)


def chip_row(row: dict) -> dict:
    """Map one stored row onto the role-aware ChIP columns.

    Old datasets carry mark/condition/replicate; they are still readable, and
    the group name is rebuilt as ``mark_condition`` so the page shows the same
    grouping the analysis used. A stored row without a role is an IP row:
    only the role-aware schema can hold an Input row, and the missing ones are
    added by ``add_input_rows``.
    """
    if "mark" in row:
        mark = (row.get("mark") or "").strip()
        condition = (row.get("condition") or "").strip()
        replicate = (row.get("replicate") or "").strip()
        group = (
            f"{mark}_{condition}" if condition and condition != "single" else mark
        )
        return {
            "sample": row.get("sample", ""),
            "role": "IP",
            "sample_name": group,
            "condition": f"rep{replicate or 1}",
            "input_sample": row.get("input_sample") or row.get("input_id", ""),
        }
    role = (row.get("role") or "").strip().lower()
    return {
        "sample": row.get("sample", ""),
        "role": "Input" if role == "input" else ("IP" if role == "ip" else ""),
        "sample_name": row.get("sample_name", ""),
        "condition": row.get("condition", ""),
        "input_sample": row.get("input_sample", ""),
    }


def input_key(row: dict) -> str:
    """Group an Input row belongs to; the id is the fallback while unnamed."""
    return (row.get("sample_name") or "").strip() or (row.get("sample") or "").strip()


def normalize_rows(omics: str, rows: list[dict]) -> list[dict]:
    """Give every ChIP row a role and every Input a replicate number.

    Drafts and files written before the role column existed hold a row for an
    Input with no group of its own; a row cited as someone's control, or one
    whose group is named ``input*``, is that Input. Two libraries of one Input
    are told apart by ``condition``, so a blank one is numbered instead of
    being saved as a second rep1.
    """
    if omics != "ChIPseq":
        return rows
    cited = {
        (row.get("input_sample") or row.get("input_id") or "").strip() for row in rows
    } - {""}
    out = []
    for row in rows:
        clean = {key: str(row.get(key, "") or "") for key in datasets.columns(omics)}
        # Read the role from the row itself: a running app can still hold the
        # column list from before the role column existed.
        clean["role"] = str(row.get("role", "") or "")
        raw = clean["role"].strip().lower()
        if raw in ("ip", "input"):
            clean["role"] = "IP" if raw == "ip" else "Input"
        else:
            clean["role"] = (
                "Input" if clean["sample"] in cited and not clean["sample_name"] else "IP"
            )
        # The old table had no role column and no room for an Input row, so a
        # control was written as an IP row whose group is named ``input*``.
        if clean["role"] == "IP" and clean["sample_name"].strip().lower().startswith("input"):
            clean["role"] = "Input"
        out.append(clean)

    taken: dict[str, set[str]] = {}
    for row in out:
        if row["role"] == "Input" and row["condition"]:
            taken.setdefault(input_key(row), set()).add(row["condition"].lower())
    for row in out:
        if row["role"] != "Input" or row["condition"]:
            continue
        used = taken.setdefault(input_key(row), set())
        for number in range(1, 10):
            candidate = f"rep{number}"
            if candidate not in used:
                row["condition"] = candidate
                used.add(candidate)
                break
    return out


def add_input_rows(
    rows: list[dict], samples: list[str], names: dict[str, str] | None = None
) -> list[dict]:
    """A row for every Input cited by an IP row but absent from the table.

    Five-column tables named the Input only in ``input_sample``, so reopening
    one would otherwise fail the check that every cited Input has its own row.
    Only a library the step-1 selection holds gets a row; a hand-typed id with
    no library behind it is left for validation to report.
    """
    present = {(row.get("sample") or "").strip() for row in rows}
    rank = {sample: index for index, sample in enumerate(samples)}
    out = list(rows)
    for row in list(rows):
        ref = (row.get("input_sample") or "").strip()
        if not ref or ref in present or ref not in rank:
            continue
        present.add(ref)
        added = {
            "sample": ref,
            "role": "Input",
            "sample_name": (names or {}).get(ref) or ref,
            "condition": "",
            "input_sample": "",
        }
        position = len(out)
        for index, other in enumerate(out):
            if rank.get((other.get("sample") or "").strip(), len(rank)) > rank[ref]:
                position = index
                break
        out.insert(position, added)
    return normalize_rows("ChIPseq", out)


def draft_rows(omics: str, name: str, samples: list[str]) -> list[dict]:
    stored = {row.get("sample", ""): row for row in datasets.metadata_rows(omics, name)}
    # Older tables kept the Input name on the IP row; the role-aware schema
    # reads it from the Input row, so the name has to travel with the new row.
    legacy_names: dict[str, str] = {}
    cited: set[str] = set()
    for row in stored.values():
        ref = (row.get("input_sample") or row.get("input_id") or "").strip()
        value = (row.get("input_sample_name") or "").strip()
        if not ref:
            continue
        cited.add(ref)
        if value:
            legacy_names.setdefault(ref, value)
    # After a save, the table itself decides which rows exist: a sample the
    # user deleted (chosen at step 1 but absent from the saved rows) must not
    # reappear on the next visit. Pick a different set of samples at step 1 and
    # every row is drafted again.
    table = datasets.descriptors(omics, name).get("table") or {}
    saved = {str(item) for item in (table.get("rows") or [])}
    chosen = {str(item) for item in (table.get("chosen") or [])}
    if saved and set(samples) == chosen:
        picked = [sample for sample in samples if sample in saved]
    else:
        picked = list(samples)
    rows = []
    for sample in picked:
        base = stored.get(sample) or {}
        if omics == "ChIPseq":
            if base:
                row = chip_row(base)
            elif sample in cited:
                # An older file named its controls only on the IP rows; the
                # library still holds a row of its own in this schema.
                row = {
                    "sample": sample,
                    "role": "Input",
                    "sample_name": legacy_names.get(sample, ""),
                    "condition": "",
                    "input_sample": "",
                }
            else:
                row = {
                    "sample": sample,
                    "role": "IP",
                    "sample_name": "",
                    "condition": "",
                    "input_sample": "",
                }
            row["sample"] = sample
            rows.append({key: row.get(key, "") for key in datasets.columns(omics)})
        else:
            rows.append(
                {
                    "sample": sample,
                    "sample_name": base.get("sample_name", ""),
                    "group": base.get("group", ""),
                }
            )
    if omics == "ChIPseq":
        rows = add_input_rows(normalize_rows(omics, rows), samples, legacy_names)
    return rows


def drop_blank_rows(rows: list[dict]) -> list[dict]:
    """Rows with no cell filled in.

    The editor's add-row button appends an empty row; saving it would fail with
    "every row needs a sample" before the user ever noticed it.
    """
    return [row for row in rows if any(str(value or "").strip() for value in row.values())]


def draft_path(omics: str, name: str) -> Path:
    return datasets.directory(omics, name) / DRAFT_NAME


def read_draft(omics: str, name: str) -> dict:
    path = draft_path(omics, name)
    if not path.is_file():
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError):
        return {}
    return data if isinstance(data, dict) else {}


def write_draft(omics: str, name: str, payload: dict) -> Path:
    """Keep the in-progress selection and table on disk.

    A browser refresh opens a new Streamlit session and session_state is gone,
    so without this file every reload asks the user to redo step 1 and step 2.
    """
    path = draft_path(omics, name)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = dict(payload)
    body["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    path.write_text(
        yaml.safe_dump(body, sort_keys=False, allow_unicode=True, default_flow_style=False),
        encoding="utf-8",
    )
    return path


def clear_draft(omics: str, name: str) -> None:
    path = draft_path(omics, name)
    if path.is_file():
        path.unlink()


def group_summary(omics: str, rows: list[dict]) -> list[tuple[str, int]]:
    counts: dict[str, int] = {}
    for row in rows:
        if omics == "ChIPseq":
            # An Input row names its own control group, which no comparison
            # can cite; counting it here would offer it as a contrast member.
            if (row.get("role") or "").strip().lower() == "input":
                continue
            # The group is the sample_name; a blank one belongs to no group and
            # ``missing_fields`` already names the gap, so counting it as a
            # group of one would report the same gap twice.
            key = (row.get("sample_name") or "").strip()
        else:
            key = (row.get("group") or "").strip()
        counts[key] = counts.get(key, 0) + 1
    return sorted(counts.items())


def group_notes(omics: str, rows: list[dict]) -> list[str]:
    """Table problems that save-time validation would reject.

    Mark, condition, group and replicate are typed by hand. A duplicate
    replicate or a name differing only by case splits one group into two and
    only shows up in the analysis, so the page reports it while it is typed.
    """
    notes: list[str] = []
    for label, count in group_summary(omics, rows):
        if not label or count >= 2:
            continue
        notes.append(i18n.t("prep.note.single", label=label, n=i18n.number(count)))

    folded: dict[str, dict[str, None]] = {}
    replicates: dict[tuple[str, str, str], list[str]] = {}
    for row in rows:
        sample = (row.get("sample") or "").strip() or "?"
        if omics == "ChIPseq":
            # An Input group shares the sample_name column with the IP groups
            # but is never a contrast; the two name spaces are checked apart.
            role = (row.get("role") or "").strip().lower() or "ip"
            group = (row.get("sample_name") or "").strip()
            condition = (row.get("condition") or "").strip()
            if group:
                folded.setdefault(f"{role}:{group.casefold()}", {})[group] = None
            if group and condition in datasets.CHIP_REPLICATES:
                replicates.setdefault((role, group, condition), []).append(sample)
            elif condition:
                notes.append(i18n.t("prep.note.condition", sample=sample, value=condition))
        else:
            group = (row.get("group") or "").strip()
            if group:
                folded.setdefault(group.casefold(), {})[group] = None

    for variants in folded.values():
        if len(variants) > 1:
            notes.append(
                i18n.t("prep.note.case", names=" / ".join(sorted(variants)))
            )
    for (_role, group, condition), names in sorted(replicates.items()):
        if len(names) > 1:
            notes.append(
                i18n.t(
                    "prep.note.replicate",
                    group=group,
                    condition=condition,
                    names=", ".join(names),
                )
            )
    return notes


def missing_fields(omics: str, rows: list[dict]) -> list[str]:
    hints = []
    for row in rows:
        sample = (row.get("sample") or "").strip() or "?"
        if omics == "RNAseq":
            fields = ("sample_name", "group")
        elif (row.get("role") or "").strip().lower() == "input":
            # An Input names its own group; it is nobody's control library.
            fields = ("sample_name",)
        else:
            fields = ("sample_name", "condition", "input_sample")
        gaps = [field for field in fields if not (row.get(field) or "").strip()]
        if gaps:
            hints.append(i18n.t("prep.gap", sample=sample, fields=", ".join(gaps)))
    return hints


def unmatched_fastq(omics: str, name: str, rows: list[dict]) -> list[str]:
    """Rows whose sample or input_sample has no FASTQ in the dataset directory.

    The sample column is editable in the naming table, so a mistyped id would
    otherwise save cleanly and only fail later in 1_prepare. A GEO dataset is
    described before its download starts, so with no FASTQ at all there is
    nothing to compare the names against yet.
    """
    table = datasets.fastq_table(omics, name)
    if not table:
        return []
    # A cited Input that holds its own row is checked by that row, so the IP
    # row does not repeat the same gap a second time.
    samples = {(row.get("sample") or "").strip() for row in rows}
    missing = []
    for row in rows:
        sample = (row.get("sample") or "").strip()
        if not sample:
            continue
        key = (row.get("fastq_prefix") or "").strip() or sample
        if not datasets.ws.lookup_fastq(table, key):
            missing.append(sample)
            continue
        if omics != "ChIPseq":
            continue
        ref = (row.get("input_sample") or row.get("input_id") or "").strip()
        if ref and ref not in samples and not datasets.ws.lookup_fastq(table, ref):
            missing.append(i18n.t("prep.unmatched.input", sample=sample, ref=ref))
    return missing


def download_roles(omics: str, name: str, runs: list[str]) -> tuple[int, int, int]:
    """IP, Input and untabled runs among the step-1 selection.

    The download takes every run written to ``geo/SRR.txt``; the table only
    decides how the FASTQ is analysed. Counting the roles here makes that
    difference visible before the download starts.
    """
    rows = datasets.metadata_rows(omics, name)
    if not rows:
        # The download can start right after step 1, before step 2 is saved; the
        # draft is what the page shows, so count its roles rather than calling
        # every run untabled.
        rows = normalize_rows(omics, read_draft(omics, name).get("rows") or [])
    roles: dict[str, str] = {}
    for row in rows:
        sample = (row.get("sample") or "").strip()
        if not sample:
            continue
        role = "input" if (row.get("role") or "").strip().lower() == "input" else "ip"
        roles[sample] = role
    ip = sum(1 for run in runs if roles.get(run) == "ip")
    controls = sum(1 for run in runs if roles.get(run) == "input")
    return ip, controls, len(runs) - ip - controls


def status(omics: str, name: str) -> dict:
    info = datasets.summary(omics, name)
    rows = datasets.metadata_rows(omics, name)
    errors, warnings = (
        datasets.validate_metadata(omics, rows)
        if rows
        else ([i18n.t("prep.error.empty_metadata")], [])
    )
    return {
        "samples": info["samples"],
        "inputs": info.get("inputs", 0),
        "pairs": info["fastq_pairs"],
        "singles": info["fastq_singles"],
        "missing": info["missing_fastq"],
        "errors": errors,
        "warnings": warnings,
        "ready": bool(rows) and not errors and not info["missing_fastq"],
    }


def complete(
    omics: str,
    name: str,
    rows: list[dict],
    source: str,
    chosen: list[str] | None = None,
) -> tuple[list[str], list[str], Path | None]:
    rows = normalize_rows(omics, rows)
    errors, warnings = datasets.validate_metadata(omics, rows)
    missing = unmatched_fastq(omics, name, rows)
    if missing:
        # Raw data is fetched after step 1, so a missing FASTQ is a download
        # still to do, not a reason to reject the table; the Run page checks it
        # again before anything is submitted.
        pending = [item for item in missing if item in datasets.srr_list(omics, name)]
        absent = [item for item in missing if item not in pending]
        if pending:
            warnings.append(i18n.t("prep.warning.pending", names=", ".join(pending[:8])))
        if absent:
            warnings.append(i18n.t("prep.warning.absent", names=", ".join(absent[:8])))
    if errors:
        return errors, warnings, None
    path = datasets.save_metadata(omics, name, rows)
    descriptor = dict(datasets.descriptors(omics, name))
    descriptor.setdefault("source", source)
    controls = (
        sum(1 for row in rows if (row.get("role") or "").strip().lower() == "input")
        if omics == "ChIPseq"
        else 0
    )
    descriptor.update(
        {
            "samples": len(rows) - controls,
            "inputs": controls,
            "prepared": True,
        }
    )
    descriptor["table"] = {
        "chosen": list(chosen if chosen is not None else (row.get("sample") for row in rows)),
        "rows": [row.get("sample") for row in rows],
    }
    datasets.save_descriptor(omics, name, descriptor)
    return [], warnings, path


def prepared(omics: str, name: str) -> bool:
    """The naming table is saved and accepted.

    Raw data plays no part: a GEO dataset is described first and downloaded
    afterwards, so requiring the FASTQ here would make the download button
    wait for the download.
    """
    rows = datasets.metadata_rows(omics, name)
    if not rows:
        return False
    errors, _ = datasets.validate_metadata(omics, rows)
    return not errors


def download_scope(omics: str) -> str:
    """Job slot for one omics' download; the progress file shares the name."""
    return f"download-{omics}"


def download_progress_file(omics: str) -> Path:
    return env_paths.RUN_STATE / "progress" / f"{download_scope(omics)}.json"


def begin_download(omics: str) -> Path:
    """Drop the previous run's file before starting a new download.

    A file left over from the last run would be drawn as this run's progress
    until the script writes its own, which reads as a bar that starts full.
    """
    path = download_progress_file(omics)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        path.unlink()
    return path


def read_download_progress(omics: str) -> dict:
    path = download_progress_file(omics)
    if not path.is_file():
        return {}
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return state if isinstance(state, dict) else {}


def download_job(omics: str) -> dict:
    """Record for this omics' download, including the shared slot older runs used."""
    record = runtime.job_status(download_scope(omics))
    if record:
        return record
    legacy = runtime.job_status("download")
    if str(legacy.get("label") or "").startswith(f"download-{omics}-"):
        return legacy
    return {}


def download_command(omics: str, name: str, runs: list[str], progress: Path) -> list[str]:
    return [
        sys.executable,
        str(env_paths.SCRIPT_ROOT / "common" / "fetch_sra.py"),
        "--omics",
        omics,
        "--dataset",
        name,
        "--srr",
        *runs,
        "--progress-file",
        str(progress),
    ]


def downloaded(omics: str, name: str, runs: list[str]) -> tuple[int, list[str]]:
    """How many of the listed runs already hold FASTQ, and which do not."""
    table = datasets.fastq_table(omics, name)
    missing = [run for run in runs if not datasets.ws.lookup_fastq(table, run)]
    return len(runs) - len(missing), missing


def expected_bytes(omics: str, name: str, runs: list[str]) -> int:
    """Archive size of the listed runs, from the GEO run table of step 1."""
    sizes: dict[str, int] = {}
    for row in datasets.read_run_info(omics, name):
        run = (row.get("Run") or "").strip()
        try:
            megabytes = float(row.get("size_MB") or 0)
        except ValueError:
            continue
        if run and megabytes > 0:
            sizes[run] = int(megabytes * 1000 * 1000)
    return sum(sizes.get(run, 0) for run in runs)


def _number(value: object) -> float:
    """One numeric field of the progress file; a garbled one must not hide the bar."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def download_view(state: dict) -> dict:
    """Bar fraction, headline, stage rows and per-run states for one download.

    Same file contract as ``annotation.progress_view``; the labels and the
    byte counts differ.
    """
    steps = [str(step) for step in (state.get("steps") or [])]
    if not steps:
        return {
            "overall": 0.0,
            "text": "",
            "rows": [],
            "runs": [],
            "done": False,
            "finished": 0,
            "total": 0,
            "downloaded": 0,
            "expected": 0,
            "updated_at": "",
            "error": "",
        }
    weights = [1.0] * len(steps)
    for position, value in enumerate(state.get("weights") or []):
        if position >= len(weights):
            break
        weights[position] = max(0.0, _number(value))
    done = bool(state.get("done"))
    error = str(state.get("error") or "")
    index = int(max(0, min(_number(state.get("step_index")), len(steps) - 1)))
    fraction = max(0.0, min(1.0, _number(state.get("step_progress"))))
    finished = int(_number(state.get("finished")))
    total = int(_number(state.get("total")))

    if done:
        # A download that ended with failures is not a full bar; the status
        # line above reports what went wrong.
        overall = 1.0 if not error else (finished / total if total else 0.0)
    else:
        total_weight = sum(weights) or 1.0
        overall = (sum(weights[:index]) + weights[index] * fraction) / total_weight

    rows = []
    for position, step in enumerate(steps):
        label = i18n.t(DOWNLOAD_KEY[step]) if step in DOWNLOAD_KEY else step
        # ``done`` means the run is over; the last stage need not be the one the
        # script stopped on when every sample was already in place.
        if done:
            badge = (
                (i18n.t("state.done"), "ok")
                if not error
                else (i18n.t("prep.state.ended"), "warn")
            )
            rows.append((label, "", *badge))
        elif position < index:
            rows.append((label, "", i18n.t("state.done"), "ok"))
        elif position == index:
            rows.append(
                (label, str(state.get("detail") or ""), i18n.t("state.running"), "run")
            )
        else:
            rows.append((label, "", i18n.t("state.idle"), "idle"))

    if done and error:
        text = i18n.t(
            "prep.text.partial",
            finished=i18n.number(finished),
            total=i18n.number(total),
        )
    elif done:
        text = i18n.t("prep.text.done", detail=state.get("message") or "")
    else:
        text = f"{i18n.percent(overall)} · {i18n.t(DOWNLOAD_KEY[steps[index]]) if steps[index] in DOWNLOAD_KEY else steps[index]}"
        detail = str(state.get("detail") or "")
        if detail:
            text = f"{text} · {detail}"

    runs = []
    for item in state.get("runs") or []:
        if not isinstance(item, dict):
            continue
        stage = str(item.get("stage") or "")
        size = int(_number(item.get("bytes")))
        runs.append(
            {
                i18n.t("prep.column.sample"): str(item.get("run") or ""),
                i18n.t("prep.column.run_state"): i18n.t(RUN_STAGE_KEY[stage])
                if stage in RUN_STAGE_KEY
                else stage,
                i18n.t("prep.column.downloaded"): i18n.size(size) if size else "",
                i18n.t("prep.column.error"): str(item.get("error") or ""),
            }
        )

    return {
        "overall": overall,
        "text": text,
        "rows": rows,
        "runs": runs,
        "done": done,
        "finished": finished,
        "total": total,
        "downloaded": int(_number(state.get("downloaded"))),
        "expected": int(_number(state.get("expected"))),
        "updated_at": str(state.get("updated_at") or ""),
        "error": error,
    }


def file_view(omics: str, name: str, runs: list[str]) -> dict:
    """Coarse stand-in for a download started before the progress file existed.

    Only the FASTQ on disk can be counted, so the bar moves per sample instead
    of per byte; a job launched by an older page still shows movement.
    """
    total = len(runs)
    have, awaiting = downloaded(omics, name, runs)
    waiting = set(awaiting)
    done = not awaiting
    return {
        "overall": (have / total) if total else 0.0,
        "text": i18n.t(
            "prep.text.fallback",
            finished=i18n.number(have),
            total=i18n.number(total),
            mode=i18n.t("prep.text.fallback_done" if done else "prep.text.fallback_estimate"),
        ),
        "rows": [
            (
                i18n.t("prep.row.raw"),
                i18n.t(
                    "prep.row.raw_detail",
                    finished=i18n.number(have),
                    total=i18n.number(total),
                ),
                i18n.t("state.done") if done else i18n.t("state.running"),
                "ok" if done else "run",
            )
        ],
        "runs": [
            {
                i18n.t("prep.column.sample"): run,
                i18n.t("prep.column.run_state"): i18n.t(
                    "prep.run.download" if run in waiting else "prep.run.done"
                ),
                i18n.t("prep.column.downloaded"): "",
                i18n.t("prep.column.error"): "",
            }
            for run in runs
        ],
        "done": done,
        "finished": have,
        "total": total,
        "downloaded": 0,
        "expected": expected_bytes(omics, name, runs),
        "updated_at": "",
        "error": "",
    }
