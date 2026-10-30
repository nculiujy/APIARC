#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common import geo, workspace

SRA_SUFFIX = ".sra"
# Direct object URL for one run. The bucket serves the same archive prefetch
# fetches, without the resolver in front of it.
SRA_ODP = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"

# The three steps a run goes through. Weights, not step counts, drive the bar:
# the download is the network wait and owns most of the wall clock.
STEPS = ("download", "convert", "compress")
WEIGHTS = (6, 3, 1)
# Where a run sits on that scale; ``skip``/``fail`` are terminal like ``done``.
STAGE_INDEX = {
    "queued": 0,
    "download": 0,
    "convert": 1,
    "compress": 2,
    "done": 3,
    "skip": 3,
    "fail": 3,
}
POLL_SECONDS = 1.5


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download SRA runs into a dataset directory")
    parser.add_argument("--omics", required=True, choices=list(workspace.OMICS))
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--srr", nargs="*", default=[])
    parser.add_argument("--srr-file", default="")
    parser.add_argument("--data-root", default=str(workspace.DATA_ROOT))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--keep-sra", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--progress-file", default="")
    return parser.parse_args(argv)


def collect(args: argparse.Namespace) -> list[str]:
    runs = [item.strip() for item in args.srr if item.strip()]
    if args.srr_file:
        text = Path(args.srr_file).read_text().split()
        runs += [item.strip() for item in text if item.strip()]
    seen: dict[str, None] = {}
    for run in runs:
        seen.setdefault(run.upper(), None)
    return list(seen)


def expected_bytes(target: Path, runs: list[str]) -> dict[str, int]:
    """SRA size per run, from the GEO run table the panel wrote at step 1.

    ``size_MB`` is the archive size, which is what prefetch downloads, so the
    bar can show real bytes instead of a spinner. A dataset imported from disk
    has no run table and reports no sizes.
    """
    rows = geo.read_run_info(target / "geo" / "SraRunInfo.csv")
    sizes: dict[str, int] = {}
    for row in rows:
        run = (row.get("Run") or "").strip()
        try:
            megabytes = float(row.get("size_MB") or 0)
        except ValueError:
            continue
        if run and megabytes > 0:
            sizes[run] = int(megabytes * 1000 * 1000)
    return {run: sizes[run] for run in runs if run in sizes}


class Progress:
    """Stage and byte state the dashboard polls every couple of seconds.

    ``prefetch`` writes its own progress lines but they only appear in the log
    after the fact, so the script reports what it knows itself: which stage the
    slowest run is in, how many bytes reached the SRA cache, how many runs are
    finished. One JSON file, replaced atomically - half a document reads as
    "no progress at all".
    """

    def __init__(
        self,
        path: Path | str | None,
        omics: str,
        dataset: str,
        runs: list[str],
        expected: dict[str, int],
        cache: Path,
    ) -> None:
        self.path = Path(path) if path else None
        self.omics = omics
        self.dataset = dataset
        self.expected = dict(expected)
        self.cache = Path(cache)
        self.lock = threading.Lock()
        self.runs = {run: {"stage": "queued", "error": "", "bytes": 0} for run in runs}
        self.started_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.done = False
        self.error = ""
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def enter(self, run: str, stage: str) -> None:
        with self.lock:
            if run in self.runs:
                self.runs[run]["stage"] = stage
        print(f"[stage] {run}: {stage}", flush=True)

    def prefetched(self, run: str, size: int) -> None:
        with self.lock:
            if run in self.runs:
                self.runs[run]["bytes"] = int(max(0, size))
                self.runs[run]["stage"] = "convert"

    def finish(self, run: str, skipped: bool = False) -> None:
        with self.lock:
            if run in self.runs:
                self.runs[run]["stage"] = "skip" if skipped else "done"

    def fail(self, run: str, message: str) -> None:
        with self.lock:
            if run in self.runs:
                self.runs[run]["stage"] = "fail"
                self.runs[run]["error"] = message

    def live_size(self, run: str) -> int:
        # prefetch writes .sra.tmp and the HTTPS path writes .sra.part, so the
        # finished name alone would report 0 B for the whole download.
        folder = self.cache / run
        size = 0
        for suffix in ("", ".tmp", ".part"):
            try:
                size = max(size, (folder / f"{run}{SRA_SUFFIX}{suffix}").stat().st_size)
            except OSError:
                continue
        return size

    def snapshot(self) -> dict:
        with self.lock:
            runs = {run: dict(state) for run, state in self.runs.items()}
        total = len(runs) or 1
        downloaded = 0
        converted = 0
        compressed = 0
        finished = 0
        open_indexes = []
        for run, state in runs.items():
            index = STAGE_INDEX.get(state["stage"], 0)
            downloaded += max(state["bytes"], self.live_size(run))
            if index >= 1:
                converted += 1
            if index >= 2 and state["stage"] != "fail":
                compressed += 1
            # A failed run also sits on the last index; counting it as finished
            # would draw a full bar for a download that fetched nothing.
            if state["stage"] in ("done", "skip"):
                finished += 1
            if index < 3:
                open_indexes.append(index)

        if self.done:
            index = len(STEPS) - 1
            fraction = 1.0
        elif open_indexes:
            index = min(open_indexes)
            fraction = {
                0: min(0.99, downloaded / sum(self.expected.values()))
                if self.expected
                else (total - open_indexes.count(0)) / total,
                1: converted / total,
                2: compressed / total,
            }.get(index, 0.0)
        else:
            index = len(STEPS) - 1
            fraction = 1.0

        expected_total = sum(self.expected.values())
        if STEPS[index] == "download" and expected_total:
            detail = f"{workspace.human_bytes(downloaded)} / {workspace.human_bytes(expected_total)}"
        elif STEPS[index] == "convert":
            detail = f"下载完成 · 转换 {converted}/{total}"
        else:
            detail = f"压缩 {compressed}/{total}"
        if finished:
            detail += f" · 完成 {finished}/{total}"

        return {
            "omics": self.omics,
            "dataset": self.dataset,
            "steps": list(STEPS),
            "weights": list(WEIGHTS),
            "step": STEPS[index],
            "step_index": index,
            "step_progress": round(max(0.0, min(1.0, float(fraction))), 4),
            "message": f"{self.dataset} · {total} 个样本",
            "detail": detail,
            "downloaded": int(downloaded),
            "expected": int(expected_total),
            "finished": finished,
            "total": total,
            "started_at": self.started_at,
            "updated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "done": self.done,
            "error": self.error,
            "runs": [
                {
                    "run": run,
                    "stage": state["stage"],
                    "bytes": int(state["bytes"]),
                    "error": state["error"],
                }
                for run, state in runs.items()
            ],
        }

    def write(self) -> None:
        if self.path is None:
            return
        payload = self.snapshot()
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            temporary = self.path.with_suffix(self.path.suffix + ".tmp")
            temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
            os.replace(temporary, self.path)
        except OSError:
            return

    def start(self) -> None:
        self.write()
        self._thread = threading.Thread(target=self._poll, daemon=True)
        self._thread.start()

    def _poll(self) -> None:
        while not self._stop.wait(POLL_SECONDS):
            self.write()

    def close(self, error: str = "") -> None:
        self.error = error
        self.done = True
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5)
        self.write()


def existing_fastq(target: Path, run: str) -> bool:
    return any(target.glob(f"{run}*.fastq.gz")) or any(target.glob(f"{run}*.fq.gz"))


def https_sra(cache: Path, run: str) -> Path | None:
    """Fetch the archive from the ODP bucket, or return None to fall back.

    prefetch goes through sra-download.ncbi.nlm.nih.gov, which for a single run
    can sit at tens of KB/s while the same object comes from the bucket at
    megabytes per second on the same network. prefetch is still the fallback:
    it resolves runs the bucket does not list.
    """
    dest = cache / run / f"{run}{SRA_SUFFIX}"
    if dest.is_file():
        return dest
    if not shutil.which("curl"):
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    command = ["curl", "-fsSL", "--retry", "3", "--retry-delay", "5"]
    command += ["--connect-timeout", "30", "--max-time", "7200"]
    if part.is_file() and part.stat().st_size:
        # A rerun continues where the last attempt stopped instead of starting
        # the whole archive again.
        command += ["-C", "-"]
    command += ["-o", str(part), SRA_ODP.format(run=run)]
    try:
        subprocess.run(command, check=True)
    except (subprocess.CalledProcessError, OSError) as error:
        print(f"[warn] {run}: direct fetch failed ({error}); trying prefetch", flush=True)
        part.unlink(missing_ok=True)
        return None
    part.replace(dest)
    return dest


def prefetch(target: Path, cache: Path, run: str) -> Path:
    expected = cache / run / f"{run}{SRA_SUFFIX}"
    if expected.is_file():
        return expected
    fetched = https_sra(cache, run)
    if fetched is not None:
        return fetched
    command = ["prefetch", "--max-size", "100G", "-O", str(cache), run]
    subprocess.run(command, check=True)
    if expected.is_file():
        return expected
    flat = cache / f"{run}{SRA_SUFFIX}"
    if flat.is_file():
        return flat
    raise RuntimeError(f"prefetch finished without {run}.sra")


def dump(target: Path, sra: Path, run: str, threads: int) -> list[Path]:
    command = [
        "fasterq-dump",
        "--split-files",
        "--skip-technical",
        "--threads",
        str(max(1, threads)),
        "-O",
        str(target),
        str(sra),
    ]
    subprocess.run(command, check=True)
    produced = sorted(target.glob(f"{run}*.fastq"))
    if not produced:
        produced = sorted(target.glob(f"{run}*.fq"))
    return produced


def compress(paths: list[Path], threads: int) -> list[Path]:
    done = []
    for path in paths:
        command = ["pigz", "-p", str(max(1, threads))]
        if not shutil.which("pigz"):
            command = ["gzip"]
        subprocess.run(command + [str(path)], check=True)
        done.append(path.with_suffix(path.suffix + ".gz"))
    return done


def handle(
    run: str,
    target: Path,
    cache: Path,
    threads: int,
    keep_sra: bool,
    progress: Progress,
) -> str:
    if existing_fastq(target, run):
        progress.finish(run, skipped=True)
        return f"[skip] {run}: fastq already present"
    progress.enter(run, "download")
    sra = prefetch(target, cache, run)
    try:
        progress.prefetched(run, sra.stat().st_size)
    except OSError:
        progress.prefetched(run, 0)
    progress.enter(run, "convert")
    produced = dump(target, sra, run, threads)
    if not produced:
        raise RuntimeError(f"fasterq-dump produced no FASTQ for {run}")
    progress.enter(run, "compress")
    compress(produced, threads)
    if not keep_sra:
        sra.unlink(missing_ok=True)
        shutil.rmtree(sra.parent, ignore_errors=True)
    progress.finish(run)
    return f"[done] {run}: {len(produced)} file(s)"


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    runs = collect(args)
    if not runs:
        print("[ERROR] no accession given", file=sys.stderr)
        return 2

    target = workspace.dataset_dir(args.omics, args.dataset, args.data_root)
    target.mkdir(parents=True, exist_ok=True)
    cache = target / ".cache" / "sra"
    cache.mkdir(parents=True, exist_ok=True)

    progress = Progress(
        args.progress_file, args.omics, args.dataset, runs, expected_bytes(target, runs), cache
    )
    progress.start()

    failures: list[str] = []
    try:
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = {
                pool.submit(
                    handle, run, target, cache, args.threads, args.keep_sra, progress
                ): run
                for run in runs
            }
            # As they finish, not in submission order: a slow first run must not
            # keep the other runs' results out of the log.
            for future in as_completed(futures):
                run = futures[future]
                try:
                    print(future.result(), flush=True)
                except Exception as error:  # noqa: BLE001
                    failures.append(f"{run}: {error}")
                    progress.fail(run, str(error))
                    print(f"[ERROR] {run}: {error}", file=sys.stderr, flush=True)
    finally:
        if not args.keep_sra:
            shutil.rmtree(cache.parent, ignore_errors=True)
        progress.close(error="; ".join(failures))

    print(f"[Info] dataset={target} requested={len(runs)} failed={len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
