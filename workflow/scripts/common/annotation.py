"""Reference annotation registry: assembly (version) and source.

One project runs both omics against one annotation. Version and source are
project-level settings, so RNA-seq and ChIP-seq cannot drift apart; every result
path carries the same ``<assembly>-<source>`` segment.

The dashboard imports this module too, so the build page and the pipeline agree
on file names, index prefixes and the TxDb used for peak annotation.
"""

from __future__ import annotations

import os

# The mock build in scripts/mock-hg19.sh writes to a scratch root so that a
# synthetic reference can never be mistaken for a downloaded one by the dashboard.
ANNO_ROOT = os.environ.get("APIARC_ANNO_ROOT") or os.path.join(
    "workflow", "resources", "anno"
)

KIND_DIRS = {"rna": "RNAseq_anno", "chip": "ChIPseq_anno"}

ASSEMBLIES = {
    "GRCh37": {
        "species": "homo",
        "taxon": "human",
        "orgdb": "org.Hs.eg.db",
        "sources": {
            "gencode": {
                "genome": "GRCh37.primary_assembly.genome.fa",
                "gtf": "gencode.v44lift37.annotation.gtf",
                "index": "grch37",
                "release": "GENCODE human release 44 (GRCh37 lift37)",
                "base_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping",
                "txdb": "TxDb.Hsapiens.UCSC.hg19.knownGene",
            }
        },
    },
    "GRCh38": {
        "species": "homo",
        "taxon": "human",
        "orgdb": "org.Hs.eg.db",
        "sources": {
            "gencode": {
                "genome": "GRCh38.primary_assembly.genome.fa",
                "gtf": "gencode.v44.annotation.gtf",
                "index": "grch38",
                "release": "GENCODE human release 44",
                "base_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44",
                "txdb": "TxDb.Hsapiens.UCSC.hg38.knownGene",
            }
        },
    },
    "GRCm38": {
        "species": "mm",
        "taxon": "mouse",
        "orgdb": "org.Mm.eg.db",
        "sources": {
            "gencode": {
                "genome": "GRCm38.primary_assembly.genome.fa",
                "gtf": "gencode.vM25.annotation.gtf",
                "index": "grcm38",
                "release": "GENCODE mouse release M25",
                "base_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25",
                "txdb": "TxDb.Mmusculus.UCSC.mm10.knownGene",
            }
        },
    },
}

SPECIES_ASSEMBLY = {
    "mm": "GRCm38",
    "homo": "GRCh38",
    "mouse": "GRCm38",
    "human": "GRCh38",
}

ALIASES = {
    "hg19": "GRCh37",
    "hg38": "GRCh38",
    "mm10": "GRCm38",
    "grch37": "GRCh37",
    "grch38": "GRCh38",
    "grcm38": "GRCm38",
    "mm": "GRCm38",
    "homo": "GRCh38",
}


def assembly_names() -> list[str]:
    return list(ASSEMBLIES)


def canonical_assembly(name: str) -> str:
    text = str(name or "").strip()
    key = text if text in ASSEMBLIES else ALIASES.get(text.lower(), "")
    if not key:
        raise ValueError(
            f"unknown assembly '{name}'; supported: {', '.join(assembly_names())}"
        )
    return key


def sources(assembly: str) -> list[str]:
    return list(ASSEMBLIES[canonical_assembly(assembly)]["sources"])


def default_source(assembly: str) -> str:
    return sources(assembly)[0]


def annotation_id(assembly: str, source: str) -> str:
    return f"{canonical_assembly(assembly)}-{source}"


def assembly_of_species(species: str) -> str:
    return SPECIES_ASSEMBLY.get(str(species or "").strip(), "GRCm38")


def single(assembly: str, source: str = "") -> dict:
    """Resolve one annotation without reading a config; accepts species aliases."""
    return resolve({"annotation": {"assembly": assembly, "source": source}})


def _resolve_source(assembly: str, source: str) -> str:
    available = sources(assembly)
    text = str(source or "").strip() or available[0]
    if text not in available:
        raise ValueError(
            f"assembly {assembly} has no source '{source}'; available: {', '.join(available)}"
        )
    return text


def declared(config: dict) -> tuple[str, str]:
    """Assembly and source configured for the project (legacy ``species`` accepted)."""
    block = config.get("annotation") or {}
    raw_assembly = str(block.get("assembly") or "").strip()
    if raw_assembly:
        assembly = canonical_assembly(raw_assembly)
    else:
        assembly = assembly_of_species(str(config.get("species") or "").strip())
    return assembly, _resolve_source(assembly, str(block.get("source") or ""))


def gene_bed_name(gtf_name: str) -> str:
    return gtf_name[:-4] + ".bed" if gtf_name.endswith(".gtf") else gtf_name + ".bed"


def source_dir(kind: str, assembly: str) -> str:
    return os.path.join(ANNO_ROOT, KIND_DIRS[kind], canonical_assembly(assembly))


def resolve(config: dict) -> dict:
    """Everything the pipeline and the dashboard need for one annotation."""
    assembly, source = declared(config)
    meta = ASSEMBLIES[assembly]
    spec = meta["sources"][source]
    rna_dir = source_dir("rna", assembly)
    chip_dir = source_dir("chip", assembly)
    genome = os.path.join(chip_dir, spec["genome"])
    return {
        "id": annotation_id(assembly, source),
        "assembly": assembly,
        "source": source,
        "species": meta["species"],
        "taxon": meta["taxon"],
        "orgdb": meta["orgdb"],
        "release": spec["release"],
        "base_url": spec["base_url"],
        "index": spec["index"],
        "txdb": spec["txdb"],
        "rna_dir": rna_dir,
        "chip_dir": chip_dir,
        "genome": genome,
        "genome_rna": os.path.join(rna_dir, spec["genome"]),
        "gtf": os.path.join(rna_dir, spec["gtf"]),
        "gene_bed": os.path.join(rna_dir, gene_bed_name(spec["gtf"])),
        "rna_index": os.path.join(rna_dir, spec["index"], "genome"),
        "chip_index": os.path.join(chip_dir, spec["index"]),
    }


def manifest(config: dict) -> dict:
    """Provenance record written next to the results."""
    info = resolve(config)
    return {
        "id": info["id"],
        "assembly": info["assembly"],
        "source": info["source"],
        "species": info["species"],
        "release": info["release"],
        "txdb": info["txdb"],
        "orgdb": info["orgdb"],
        "genome": info["genome"],
        "gtf": info["gtf"],
        "gene_bed": info["gene_bed"],
        "rna_index": info["rna_index"],
        "chip_index": info["chip_index"],
    }
