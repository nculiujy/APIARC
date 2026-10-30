#!/usr/bin/env python3
"""Check that a GTF and a genome FASTA use the same contig names.

A GTF that says ``chr1`` while the FASTA says ``1`` is not an error anywhere in
the pipeline: HISAT2 aligns nothing, featureCounts counts nothing, and the run
"completes" with an empty matrix. GENCODE lift37 annotation and the GRCh37
primary assembly are built by different pipelines, so the two conventions can
meet in one directory.

The check runs once per annotation build. When one uniform rename resolves every
mismatch (add or strip a ``chr`` prefix, ``MT``/``chrM``) the GTF is rewritten in
place; when it does not, the build stops with the offending names instead of
producing an annotation nobody can align to.
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys

CHR_ALIASES = {"MT": "chrM", "chrM": "MT"}


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def fasta_contigs(path: str) -> set[str]:
    names: set[str] = set()
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith(">"):
                names.add(line[1:].split()[0].strip())
    return names


def gtf_contigs(path: str) -> list[str]:
    """Every distinct seqname, in first-seen order."""
    seen: dict[str, None] = {}
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            name = line.split("\t", 1)[0]
            if name and name not in seen:
                seen[name] = None
    return list(seen)


def candidate_names(name: str) -> list[str]:
    options = [name]
    if name.startswith("chr"):
        options.append(name[3:])
    else:
        options.append("chr" + name)
    if name in CHR_ALIASES:
        options.append(CHR_ALIASES[name])
    return options


def plan(gtf_names: list[str], fasta_names: set[str]) -> tuple[dict[str, str], list[str]]:
    rename: dict[str, str] = {}
    missing: list[str] = []
    for name in gtf_names:
        if name in fasta_names:
            continue
        match = next((option for option in candidate_names(name) if option in fasta_names), None)
        if match is None:
            missing.append(name)
        else:
            rename[name] = match
    return rename, missing


def rewrite(path: str, rename: dict[str, str]) -> None:
    temporary = path + ".aligned"
    with open_maybe_gzip(path) as source, open(temporary, "w") as target:
        for line in source:
            if not line or line.startswith("#"):
                target.write(line)
                continue
            name, _, rest = line.partition("\t")
            target.write(f"{rename.get(name, name)}\t{rest}")
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--dry-run", action="store_true", help="report without rewriting")
    args = parser.parse_args()

    fasta_names = fasta_contigs(args.fasta)
    if not fasta_names:
        print(f"[contigs] {args.fasta} has no FASTA headers", file=sys.stderr)
        return 2

    names = gtf_contigs(args.gtf)
    rename, missing = plan(names, fasta_names)
    print(
        f"[contigs] gtf={len(names)} fasta={len(fasta_names)} "
        f"rename={len(rename)} missing={len(missing)}",
        flush=True,
    )
    if missing:
        print(
            "[contigs] contigs in the GTF but not in the FASTA: "
            + ", ".join(missing[:10])
            + (" ..." if len(missing) > 10 else ""),
            file=sys.stderr,
        )
        return 2
    if rename and not args.dry_run:
        rewrite(args.gtf, rename)
        pairs = ", ".join(f"{old}->{new}" for old, new in list(rename.items())[:5])
        print(f"[contigs] rewrote GTF contig names: {pairs}", flush=True)
    elif rename:
        print("[contigs] dry run: GTF would be rewritten", flush=True)
    else:
        print("[contigs] GTF and FASTA already agree", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
