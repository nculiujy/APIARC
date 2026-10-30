#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


class SourceCache:
    def open(self, file, mode="rb"):
        return open(file.get_path_or_uri(), mode)


class Persistence:
    def __init__(self, archive: str) -> None:
        self.conda_env_archive_path = archive


class StubWorkflow:
    def __init__(self, archive: str) -> None:
        self.conda_frontend = "conda"
        self.singularity_args: list[str] = []
        self.sourcecache = SourceCache()
        self.persistence = Persistence(archive)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--env-root", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--archive", default="")
    args = parser.parse_args()

    from snakemake.deployment.conda import Env

    env_root = Path(args.env_root).resolve()
    prefix = Path(args.prefix).resolve()
    archive = args.archive or str(prefix.parent / "env-archive")
    prefix.mkdir(parents=True, exist_ok=True)
    workflow = StubWorkflow(archive)

    specs = sorted(env_root.glob("*.yaml"))
    if not specs:
        print(f"no stage definitions under {env_root}", file=sys.stderr)
        return 1

    created = 0
    present = 0
    failed: list[str] = []
    for spec in specs:
        env = Env(workflow, env_file=str(spec), env_dir=str(prefix))
        if os.path.exists(env.address):
            present += 1
            print(f"present  {spec.name}", flush=True)
            continue
        print(f"create   {spec.name}", flush=True)
        try:
            env.create()
        except Exception as error:  # noqa: BLE001
            failed.append(spec.name)
            print(f"failed   {spec.name}: {error}", file=sys.stderr, flush=True)
            continue
        created += 1

    print(f"done: {created} created, {present} present, {len(failed)} failed")
    if failed:
        print("failed: " + ", ".join(failed), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
