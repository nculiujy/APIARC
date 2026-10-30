#!/usr/bin/env python3
"""Create empty reference-data placeholders so the DAG can be validated.

The reference annotation is deliberately not in version control: it is a few
gigabytes of GENCODE and JASPAR files. That makes ``snakemake -n`` impossible on
a fresh clone, which in turn means the only workflow check CI could run was a
Python syntax check.

This script fills exactly the paths that Snakemake resolves while *building* the
job graph (declared rule inputs), so ``snakemake -n`` works without the real
data. Files are created empty and are never valid for a real run.

Paths that are only used as rule ``params`` (bowtie2/HISAT2 index prefixes, the
picard directory) are not created: Snakemake does not check their existence at
parse time, and pretending they exist would hide a real misconfiguration.

Usage::

    python workflow/scripts/common/stub_resources.py [--force]

The script refuses to touch a path that already has content unless --force is
given, so it cannot destroy a downloaded reference set by accident.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from common import annotation, resources

ASSEMBLIES = annotation.assembly_names()

MOTIF_ROOT = os.path.join("workflow", "resources", "Motif_tf_anno")
MOTIF_MAP = os.path.join(MOTIF_ROOT, "JASPAR_CORE_ID_to_NAME.txt")
MOTIF_DIR = os.path.join(MOTIF_ROOT, "JASPAR_CORE")

STUB_MOTIF = os.path.join(MOTIF_DIR, "MA0001.1.meme")

STUB_MEME = """MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF MA0001.1 placeholder
letter-probability matrix: alength= 4 w= 1 nsites= 1 E= 0
 0.25  0.25  0.25  0.25
"""


def targets() -> list[tuple[str, str]]:
    """(path, content) for every reference file the DAG resolves."""
    items: list[tuple[str, str]] = []
    for assembly in ASSEMBLIES:
        items.append((resources.gene_bed(assembly), ""))
        items.append((resources.gtf(assembly), ""))
        items.append((resources.genome_fasta(assembly), ""))
    items.append((MOTIF_MAP, "MA0001.1\tplaceholder\n"))
    if not _has_real_motifs():
        items.append((STUB_MOTIF, STUB_MEME))
    return items


def _has_real_motifs() -> bool:
    if not os.path.isdir(MOTIF_DIR):
        return False
    for name in os.listdir(MOTIF_DIR):
        if name.endswith(".meme") and name != os.path.basename(STUB_MOTIF):
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="overwrite non-empty files")
    args = parser.parse_args()

    created = 0
    kept = 0
    for path, content in targets():
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path) and os.path.getsize(path) > 0 and not args.force:
            kept += 1
            continue
        with open(path, "w") as handle:
            handle.write(content)
        created += 1
        print(f"[stub] {path}", flush=True)

    print(
        f"[stub] created {created}, kept {kept}. "
        "These files are empty on purpose: they validate the graph, not the data.",
        flush=True,
    )
    if created:
        print(
            "[stub] replace them with real references before running the pipeline",
            flush=True,
        )


if __name__ == "__main__":
    main()
