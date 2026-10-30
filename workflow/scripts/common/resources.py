"""Reference resource paths keyed by annotation.

These helpers used to live inside individual rule files. That meant a rule file
could only use them if another rule file happened to be included first: running
``RNAseq_modules.3_RC: False`` plus ``4_DEseq: True`` failed with
``NameError: get_gtf_file is not defined``, a failure caused purely by module
ordering. Reference lookups are a property of the project, not of a rule, so they
live here and the Snakefile imports them once.

Every function accepts an assembly (``GRCh37``, ``GRCh38``, ``GRCm38``) or the
legacy species alias (``homo``, ``mm``); names and index prefixes come from
``common.annotation`` so a rule cannot drift from what the build page produced.
"""

from __future__ import annotations

from . import annotation

__all__ = [
    "chip_index",
    "rna_index",
    "gtf",
    "gene_bed",
    "genome_fasta",
    "assembly",
    "txdb",
]


def _info(name: str) -> dict:
    return annotation.single(name)


def chip_index(name: str) -> str:
    """bowtie2 index prefix for ChIP-seq alignment."""
    return _info(name)["chip_index"]


def rna_index(name: str) -> str:
    """HISAT2 index prefix for RNA-seq alignment."""
    return _info(name)["rna_index"]


def gtf(name: str) -> str:
    return _info(name)["gtf"]


def gene_bed(name: str) -> str:
    return _info(name)["gene_bed"]


def genome_fasta(name: str) -> str:
    return _info(name)["genome"]


def assembly(name: str) -> str:
    return _info(name)["assembly"]


def txdb(name: str) -> str:
    return _info(name)["txdb"]
