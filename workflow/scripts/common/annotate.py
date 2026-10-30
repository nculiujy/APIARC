"""Region-to-gene assignment shared by the integrated analysis.

A region is linked to the gene whose transcription start site is closest, and
the reported distance is the distance to that TSS. The previous implementation
measured "enhancer distance" against a file of promoter peaks, so the number it
reported was the distance to another binding site, not to a gene.

The gene model is indexed once per chromosome. Rebuilding the coordinate list for
every region turned assignment into an O(regions x genes) scan, which is why the
first version was unusable on a real peak set.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass

PROMOTER = "promoter"
DISTAL = "distal"
UNASSIGNED = "unassigned"


@dataclass(frozen=True)
class Gene:
    chrom: str
    start: int
    end: int
    symbol: str
    gene_id: str
    strand: str

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end - 1


class GeneIndex:
    """Chromosome -> genes sorted by TSS, with the TSS list precomputed."""

    def __init__(self, genes: dict[str, list[Gene]]) -> None:
        self.genes = genes
        self.tss: dict[str, list[int]] = {
            chrom: [gene.tss for gene in items] for chrom, items in genes.items()
        }

    def nearest(self, chrom: str, start: int, end: int) -> tuple[Gene | None, int]:
        items = self.genes.get(chrom)
        if not items:
            return None, -1
        positions = self.tss[chrom]
        index = bisect.bisect_left(positions, start)
        window = items[max(0, index - 5): min(len(items), index + 6)]
        best = None
        best_distance = None
        for gene in window:
            distance = interval_distance(start, end, gene.tss, gene.tss)
            if best_distance is None or distance < best_distance:
                best, best_distance = gene, distance
        return best, (-1 if best_distance is None else best_distance)


def load_genes(bed_path: str) -> GeneIndex:
    """Read a GENCODE BED (chrom, start, end, gene_id, score, strand, symbol)."""
    genes: dict[str, list[Gene]] = {}
    with open(bed_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            gene = Gene(
                chrom=fields[0],
                start=int(fields[1]),
                end=int(fields[2]),
                gene_id=fields[3],
                strand=fields[5] if fields[5] in {"+", "-"} else "+",
                symbol=fields[6] or fields[3],
            )
            genes.setdefault(gene.chrom, []).append(gene)
    for chrom in genes:
        genes[chrom].sort(key=lambda g: g.tss)
    return GeneIndex(genes)


def interval_distance(start: int, end: int, other_start: int, other_end: int) -> int:
    if end < other_start:
        return other_start - end
    if start > other_end:
        return start - other_end
    return 0


def assign_gene(
    index: GeneIndex,
    chrom: str,
    start: int,
    end: int,
    promoter_window: int,
    max_distance: int,
) -> tuple[str, int, str]:
    """Return (symbol, distance_to_tss, region_type) for one region."""
    gene, distance = index.nearest(chrom, start, end)
    if gene is None or distance < 0 or distance > max_distance:
        return ("", -1, UNASSIGNED)
    region_type = PROMOTER if distance <= promoter_window else DISTAL
    return (gene.symbol, distance, region_type)
