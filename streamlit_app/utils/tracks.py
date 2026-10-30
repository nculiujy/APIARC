"""Read bigWig signal into rows and stack one row per track.

The first version of the signal tab drew every selected track as an overlaid
area. Two filled areas cover each other, and a track with a tenth of the signal
still fills its own band from zero, so the question the tab exists for - which
sample is stronger at this locus - could not be read off the picture. One row
per track with a y range shared by default (the layout of the ChIP-seq panel in
the other project) makes the comparison direct.

The figure is plotly, not altair: the tab is read by pointing at it, and
wheel zoom, drag to pan, a spike line across all rows and one hover label per
position are plotly's defaults. Zoom only magnifies the bins already loaded,
so a window that is not on screen is typed as a gene name or ``chr:start-end``
and re-read from the bigWig files - one ``stats()`` call per track, the same
bargain a genome browser makes. The box takes a gene name as well, looked up in
the GENCODE GTF the pipeline aligned against, so a locus can be named the way
the experiment is described (``Nppa``) instead of by coordinates.

Values follow the same convention as the pipeline's deeptools stage: one mean
per equal-width bin, ``nan`` (no coverage) as 0, and smoothing applied to the
display only.

The stack ends in a gene model row read from the same GTF: without it the
signal rows give no way to tell where the gene sits relative to the peak under
it, which is the question a ChIP-seq locus view is opened to answer.
"""

from __future__ import annotations

import gzip
import math
import os
import re
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from . import env_paths, i18n

try:  # Optional launcher dependency; the tab degrades to a hint without it.
    import pyBigWig
except ImportError:  # pragma: no cover - depends on the deployment
    pyBigWig = None

try:  # The chart is drawn with plotly so it can be zoomed in the browser.
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:  # pragma: no cover - depends on the deployment
    go = None
    make_subplots = None

DEFAULT_BINS = 800
PEAK_BAND = "#f17eb2"
FOCUS_BAND = "#a8cdec"
SERIES_COLORS = ("#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b")
INPUT_COLOR = "#8a8a8a"
GRID = "#ececec"
GENE_FLANK = 5_000
WIDE_FLANK = 1_000_000
GENE_COLUMNS = ["chrom", "start", "end", "strand", "name", "gid", "biotype"]
EXON_COLUMNS = ["chrom", "start", "end", "gid"]
GENE_DTYPES = {
    "chrom": str,
    "strand": str,
    "name": str,
    "gid": str,
    "biotype": str,
}
EXON_DTYPES = {"chrom": str, "gid": str}
GENE_PLUS = "#2f5d8c"
GENE_MINUS = "#8a5a2b"


def have_pybigwig() -> bool:
    return pyBigWig is not None


def have_plotly() -> bool:
    return go is not None


def _mtime(path: str | os.PathLike) -> float:
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0.0


@lru_cache(maxsize=64)
def _sizes_cached(path: str, _mtime: float) -> tuple[tuple[str, int], ...]:
    handle = pyBigWig.open(path)
    try:
        return tuple(sorted((handle.chroms() or {}).items()))
    finally:
        handle.close()


def chrom_sizes(path: str | os.PathLike) -> dict[str, int]:
    """Contig lengths of one track, used to keep an interval inside the file."""
    if pyBigWig is None or not os.path.isfile(path):
        return {}
    return dict(_sizes_cached(str(path), _mtime(path)))


@lru_cache(maxsize=512)
def _bins_cached(path: str, chrom: str, start: int, end: int, bins: int, aggregate: str,
                 _mtime: float) -> tuple[float, ...]:
    try:
        handle = pyBigWig.open(path)
    except Exception:  # noqa: BLE001 - a broken file must not take the page down
        return (0.0,) * bins
    try:
        values = handle.stats(chrom, start, end, type=aggregate, nBins=bins)
    except Exception:  # noqa: BLE001 - unknown contig or bad bounds
        return (0.0,) * bins
    finally:
        handle.close()
    if not values:
        return (0.0,) * bins
    return tuple(0.0 if value is None or not math.isfinite(value) else float(value) for value in values)


def read_bins(
    path: str | os.PathLike,
    chrom: str,
    start: int,
    end: int,
    bins: int,
    aggregate: str = "mean",
) -> np.ndarray:
    """Mean signal of ``bins`` equal-width windows; empty windows count as 0.

    ``aggregate="max"`` is for views where one bin is hundreds of kilobases:
    averaging that far flattens every peak into the background, while the peak
    height survives in a maximum.
    """
    bins = max(1, int(bins))
    if pyBigWig is None or not os.path.isfile(path):
        return np.zeros(bins, dtype=float)
    return np.asarray(
        _bins_cached(
            str(path), chrom, int(start), int(end), bins, str(aggregate), _mtime(path)
        ),
        dtype=float,
    )


def smooth(values: np.ndarray, bin_bp: float, smooth_bp: float) -> np.ndarray:
    """Moving average sized in base pairs, not in points.

    Dividing by the count of windows that actually overlap keeps the two ends
    from sagging the way a zero-padded convolution makes them look.
    """
    kernel = int(round(max(0.0, smooth_bp) / max(bin_bp, 1e-9)))
    if kernel <= 1 or values.size < 2:
        return values
    kernel = min(kernel, values.size)
    total = np.convolve(values, np.ones(kernel), mode="same")
    count = np.convolve(np.ones_like(values), np.ones(kernel), mode="same")
    return total / np.maximum(count, 1e-9)


def nice_max(value: float) -> float:
    """Round a maximum up to 1/2/2.5/5 x 10^n so the axis reads cleanly."""
    if not math.isfinite(value) or value <= 0:
        return 1.0
    exponent = math.floor(math.log10(value))
    base = 10.0 ** exponent
    for step in (1.0, 2.0, 2.5, 5.0, 10.0):
        if value <= step * base:
            return step * base
    return 10.0 * base


def unit_for(span: int) -> tuple[str, float]:
    """Display unit and divisor for a window; coordinates stay bp on disk."""
    if span >= 200_000:
        return "Mb", 1_000_000.0
    if span >= 5_000:
        return "kb", 1_000.0
    return "bp", 1.0


def parse_region(text: str, chroms: Mapping[str, int]) -> dict:
    """Turn ``chr1:40000-60000`` into a clipped interval, or explain the format."""
    cleaned = (text or "").replace(",", "").replace(" ", "")
    if not cleaned:
        raise ValueError(i18n.t("track.error.empty_region"))
    if ":" not in cleaned:
        chrom, span = cleaned, ""
    else:
        chrom, span = cleaned.split(":", 1)
    if chrom not in chroms:
        known = ", ".join(list(chroms)[:5]) or i18n.t("track.unknown_chrom_none")
        raise ValueError(i18n.t("track.error.unknown_chrom", chrom=chrom, known=known))
    length = int(chroms[chrom])
    if not span:
        start, end = 0, length
    else:
        if "-" not in span:
            raise ValueError(i18n.t("track.error.format"))
        left, right = span.split("-", 1)
        try:
            start, end = int(left), int(right)
        except ValueError as error:
            raise ValueError(i18n.t("track.error.format")) from error
    start, end = max(0, start), min(end, length)
    if end <= start:
        raise ValueError(
            i18n.t("track.error.range", chrom=chrom, length=i18n.number(length))
        )
    return {"chrom": chrom, "start": start, "end": end, "gene": None}


def resolve_chrom(name: str, chroms: Mapping[str, int]) -> str:
    """
    Map ``1`` / ``chr1`` / ``Chr1`` onto the name the bigWig uses.

    The pipeline writes Ensembl contig names into the bigWig and GENCODE's
    GTF uses the UCSC spelling, so the two only differ in the prefix here.
    """
    if name in chroms:
        return name
    lowered = {chrom.lower(): chrom for chrom in chroms}
    for candidate in (
        name.lower(),
        "chr" + name.lower().lstrip("chr"),
        name.lower().replace("chr", "", 1),
    ):
        if candidate in lowered:
            return lowered[candidate]
    return ""


def _gtf_cache(path: Path, kind: str) -> Path:
    return (
        env_paths.RUNTIME_DIR
        / "cache"
        / "genes"
        / f"{path.name}.{kind}.{int(_mtime(path))}.tsv"
    )


def _read_cached(
    path: Path,
    kind: str,
    build: "Callable[[Path], pd.DataFrame]",
    dtypes: dict[str, str],
) -> pd.DataFrame:
    """One GTF extract, read from the runtime cache or built and stored there."""
    cache = _gtf_cache(path, kind)
    if cache.is_file():
        try:
            return pd.read_csv(cache, sep="\t", dtype=dtypes)
        except (OSError, ValueError, pd.errors.ParserError):
            pass
    frame = build(path)
    try:
        cache.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(cache, sep="\t", index=False)
        # Only this extract's older copies: the gene and exon tables share the
        # directory and must not delete each other.
        for stale in cache.parent.glob(f"{path.name}.{kind}.*.tsv"):
            if stale != cache:
                stale.unlink()
    except OSError:
        pass
    return frame


def _attr(text: str, key: str) -> str:
    found = re.search(rf'{key} "([^"]*)"', text)
    return found.group(1) if found else ""


def _read_gtf_genes(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    rows: list[tuple] = []
    with opener(path, "rt", errors="replace") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = fields[8]
            try:
                start, end = int(fields[3]) - 1, int(fields[4])
            except ValueError:
                continue
            rows.append(
                (
                    fields[0],
                    start,
                    end,
                    fields[6],
                    _attr(attributes, "gene_name") or _attr(attributes, "gene_id"),
                    _attr(attributes, "gene_id"),
                    _attr(attributes, "gene_type"),
                )
            )
    return pd.DataFrame(rows, columns=GENE_COLUMNS)


def _read_gtf_exons(path: Path) -> pd.DataFrame:
    """Exon spans per gene, deduplicated across transcripts.

    The annotation row only needs blocks to draw; which transcript an exon
    belongs to would triple the cache for no visible difference.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    rows: list[tuple] = []
    with opener(path, "rt", errors="replace") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            try:
                start, end = int(fields[3]) - 1, int(fields[4])
            except ValueError:
                continue
            rows.append((fields[0], start, end, _attr(fields[8], "gene_id")))
    return pd.DataFrame(rows, columns=EXON_COLUMNS).drop_duplicates()


def gene_index(gtf: str | os.PathLike) -> pd.DataFrame:
    """Gene-level rows of the GENCODE GTF, cached under the runtime directory.

    A mouse GTF is 1.9M lines and the tab asks for the index on every rerun;
    the gene rows alone are small, so they are written once and re-read by
    later sessions. The cache name carries the GTF's mtime, which also drops
    it after a rebuild. The parsed frame is also kept in memory for the
    lifetime of the process, since reading 55k rows still costs 70 ms.
    """
    path = Path(gtf)
    if not path.is_file():
        return pd.DataFrame(columns=GENE_COLUMNS)
    return _gene_index_cached(str(path), _mtime(path))


@lru_cache(maxsize=4)
def _gene_index_cached(key: str, _mtime_value: float) -> pd.DataFrame:
    return _read_cached(Path(key), "genes", _read_gtf_genes, GENE_DTYPES)


def exon_index(gtf: str | os.PathLike) -> pd.DataFrame:
    """Exon blocks per gene, for the annotation row under the tracks.

    Same cache scheme as the gene index: a first pass over the 1.9M-line GTF
    is kept on disk, and parsed frames stay in memory per process.
    """
    path = Path(gtf)
    if not path.is_file():
        return pd.DataFrame(columns=EXON_COLUMNS)
    return _exon_index_cached(str(path), _mtime(path))


@lru_cache(maxsize=4)
def _exon_index_cached(key: str, _mtime_value: float) -> pd.DataFrame:
    return _read_cached(Path(key), "exons", _read_gtf_exons, EXON_DTYPES)


def locus_genes(
    genes: pd.DataFrame,
    exons: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    limit: int = 12,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Genes overlapping a window plus their exon blocks, for the annotation row.

    In a 1 Mb view a mouse locus can hold dozens of genes; the widest ``limit``
    of them are drawn, which keeps the row readable and the payload small.
    """
    empty = genes.iloc[0:0], exons.iloc[0:0]
    if genes.empty:
        return empty
    hit = genes[(genes["chrom"] == chrom) & (genes["start"] < end) & (genes["end"] > start)]
    if hit.empty:
        return hit, exons.iloc[0:0]
    hit = (
        hit.assign(_span=hit["end"] - hit["start"])
        .nlargest(limit, "_span")
        .sort_values("start")
        .drop(columns="_span")
    )
    if exons.empty:
        return hit, exons
    blocks = exons[
        exons["gid"].isin(set(hit["gid"]))
        & (exons["chrom"] == chrom)
        & (exons["start"] < end)
        & (exons["end"] > start)
    ].sort_values(["gid", "start"])
    # Transcript variants repeat an exon a dozen times over; overlapping spans
    # of one gene are one block on screen, so they are merged here instead of
    # being drawn on top of each other.
    merged: list[dict] = []
    for row in blocks.itertuples(index=False):
        if merged and merged[-1]["gid"] == row.gid and row.start <= merged[-1]["end"]:
            merged[-1]["end"] = max(merged[-1]["end"], int(row.end))
        else:
            merged.append(
                {
                    "chrom": row.chrom,
                    "start": int(row.start),
                    "end": int(row.end),
                    "gid": row.gid,
                }
            )
    return hit, pd.DataFrame(merged, columns=EXON_COLUMNS)


def find_gene(genes: pd.DataFrame, name: str) -> pd.DataFrame:
    """Exact gene name, then gene ID with or without the version suffix."""
    if genes.empty:
        return genes
    lowered = name.strip().lower()
    hit = genes[genes["name"].str.lower() == lowered]
    if not hit.empty:
        return hit
    bare = lowered.split(".")[0]
    return genes[genes["gid"].str.lower().str.split(".").str[0] == bare]


def parse_locus(
    query: str,
    chroms: Mapping[str, int],
    genes: pd.DataFrame | None = None,
    flank: int = GENE_FLANK,
) -> dict:
    """
    ``chr2:180,214,794-180,224,794``, ``chr2:1500000`` or ``Nppa`` → interval.

    A gene name resolves through the GENCODE gene rows to the gene body plus
    ``flank`` either side, which is the window the peak annotation in the
    ChIP-seq result was assigned over.
    """
    text = (query or "").strip()
    if not text:
        raise ValueError(i18n.t("track.error.query_empty"))
    cleaned = text.replace(",", "").replace(" ", "")
    if ":" in cleaned:
        chrom, span = cleaned.split(":", 1)
        if "-" in span:
            return parse_region(cleaned, chroms)
        if not span.isdigit():
            raise ValueError(i18n.t("track.error.format"))
        centre = int(span)
        return parse_region(f"{chrom}:{max(0, centre - 10_000)}-{centre + 10_000}", chroms)
    if genes is None or genes.empty:
        raise ValueError(i18n.t("track.error.need_gtf"))
    hits = find_gene(genes, cleaned)
    if hits.empty:
        raise ValueError(i18n.t("track.error.no_gene", name=cleaned))
    row = hits.iloc[0]
    chrom = resolve_chrom(str(row["chrom"]), chroms)
    if not chrom:
        raise ValueError(
            i18n.t(
                "track.error.gene_off_chrom", name=row["name"], chrom=row["chrom"]
            )
        )
    start = max(0, int(row["start"]) - flank)
    end = min(int(chroms[chrom]), int(row["end"]) + flank)
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "gene": {
            "name": str(row["name"]),
            "chrom": chrom,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "strand": str(row["strand"]),
            "biotype": str(row["biotype"]),
            "hits": int(len(hits)),
        },
    }


def comparison_defaults(
    groups: Mapping[str, Sequence[str]],
    pairs: Iterable[Sequence[str]],
    limit: int = 4,
) -> list[str]:
    """
    Opening tracks: the first two groups of the ChIP-seq comparison table.

    One comparison is a reference group against a treatment group, two
    libraries each, so two groups fill the ``limit`` with four tracks and the
    difference the comparison was designed to show is on screen as soon as the
    tab opens. Input libraries come last inside a group and are skipped while
    an IP library is still available.
    """
    ordered: list[str] = []
    for pair in pairs:
        for group in pair:
            if group in groups and group not in ordered:
                ordered.append(group)
        if len(ordered) >= 2:
            break
    per_group = max(1, limit // 2)
    chosen: list[str] = []
    for group in ordered[:2]:
        chosen.extend(list(groups[group])[:per_group])
    if not chosen:
        chosen = [name for names in groups.values() for name in names][:limit]
    return chosen[:limit]



def read_peaks(paths: Iterable[str | os.PathLike]) -> list[dict]:
    """Peaks from narrowPeak / merged BED files, scored by -log10(q) when present."""
    peaks: list[dict] = []
    for path in paths:
        try:
            with open(path, "r", errors="replace") as handle:
                lines = handle.readlines()
        except OSError:
            continue
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[0].startswith(("#", "track", "browser")):
                continue
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            except ValueError:
                continue
            peaks.append(
                {
                    "sample": Path(path).parent.name,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "score": _peak_score(fields),
                }
            )
    return peaks


def _peak_score(fields: Sequence[str]) -> float:
    # narrowPeak is 1-based: column 9 is -log10(q), column 10 is where the
    # summit sits inside the peak — a distance in bp, not a strength, so it
    # must never be read as the score. The merged BED only carries a constant
    # score in column 5, hence the fallbacks.
    for index in (8, 7, 4):
        if len(fields) > index:
            try:
                return float(fields[index] or 0)
            except ValueError:
                continue
    return 0.0


def peaks_in(
    peaks: Sequence[dict], chrom: str, start: int, end: int, limit: int = 40
) -> list[dict]:
    """Overlapping peaks merged by span, so one locus paints one band."""
    hit = [dict(peak) for peak in peaks if peak["chrom"] == chrom and peak["start"] < end and peak["end"] > start]
    hit.sort(key=lambda peak: peak["start"])
    merged: list[dict] = []
    for peak in hit:
        if merged and peak["start"] <= merged[-1]["end"]:
            merged[-1]["end"] = max(merged[-1]["end"], peak["end"])
            merged[-1]["score"] = max(merged[-1]["score"], peak["score"])
        else:
            merged.append(peak)
    merged.sort(key=lambda peak: -peak["score"])
    return merged[:limit]


def top_peaks(peaks: Sequence[dict], limit: int = 20) -> list[dict]:
    return sorted(peaks, key=lambda peak: -peak["score"])[:limit]


def spread_peaks(peaks: Sequence[dict], limit: int) -> list[dict]:
    """At most ``limit`` bands spread over the interval, strongest of each stretch.

    The score-sorted head piles every band onto the same few hotspots and
    leaves the rest of the chromosome bare, which is the opposite of what an
    overview is for.
    """
    limit = max(1, int(limit))
    if len(peaks) <= limit:
        return list(peaks)
    ordered = sorted(peaks, key=lambda peak: int(peak["start"]))
    origin = int(ordered[0]["start"])
    width = max(1, (int(ordered[-1]["end"]) - origin) // limit)
    buckets: dict[int, dict] = {}
    for peak in ordered:
        key = (int(peak["start"]) - origin) // width
        best = buckets.get(key)
        if best is None or float(peak["score"]) > float(best["score"]):
            buckets[key] = peak
    return [buckets[key] for key in sorted(buckets)][:limit]


def colors_for(names: Sequence[str]) -> dict[str, str]:
    """Grey for Input tracks, the series palette for everything else."""
    palette = iter(SERIES_COLORS)
    assigned: dict[str, str] = {}
    for name in names:
        if "input" in name.lower():
            assigned[name] = INPUT_COLOR
        else:
            assigned[name] = next(palette, INPUT_COLOR)
    return assigned


def read_tracks(
    tracks: Mapping[str, str | os.PathLike],
    *,
    chrom: str,
    start: int,
    end: int,
    bins: int = DEFAULT_BINS,
    smooth_bp: int = 0,
    aggregate: str = "mean",
) -> tuple[pd.DataFrame, list[str]]:
    """Long frame of ``track/position/signal`` plus notes about skipped or clipped tracks."""
    rows: list[dict] = []
    notes: list[str] = []
    for name, path in tracks.items():
        sizes = chrom_sizes(path)
        if chrom not in sizes:
            notes.append(i18n.t("track.note.missing_chrom", name=name, chrom=chrom))
            continue
        stop = min(end, sizes[chrom])
        begin = max(0, min(start, stop - 1))
        values = read_bins(path, chrom, begin, stop, bins, aggregate=aggregate)
        width = (stop - begin) / max(1, bins)
        if smooth_bp:
            values = smooth(values, width, smooth_bp)
        for index, value in enumerate(values):
            # The bin stands for its span, so it is placed at the middle of that
            # span; an edge would shift every peak by half a bin.
            rows.append(
                {
                    "track": name,
                    "position": begin + (index + 0.5) * width,
                    "signal": float(value),
                }
            )
        if (begin, stop) != (start, end):
            notes.append(
                i18n.t(
                    "track.note.clipped",
                    name=name,
                    length=i18n.number(sizes[chrom]),
                    start=i18n.number(begin),
                    stop=i18n.number(stop),
                )
            )
    frame = pd.DataFrame(rows, columns=["track", "position", "signal"])
    return frame, notes


def merge_resolution(
    context: pd.DataFrame, detail: pd.DataFrame, *, start: int, end: int
) -> pd.DataFrame:
    """One curve per track: the window at bin resolution, the rest at chromosome
    resolution.

    Reading only the window leaves the chromosome ends with nothing to draw, and
    zooming out of the window then shows a line that stops mid-air. Keeping the
    coarse bins outside the window means the curve always reaches the ends of
    the chromosome, and the bins the window covers are replaced by the finer one.
    """
    if detail.empty or context.empty:
        return context if detail.empty else detail
    outside = context[(context["position"] < start) | (context["position"] > end)]
    merged = pd.concat([outside, detail], ignore_index=True)
    return merged.sort_values(["track", "position"], ignore_index=True)


def shared_top(frame: pd.DataFrame) -> float:
    if frame.empty:
        return 1.0
    return nice_max(float(frame["signal"].max()))


def peak_draw_cap(rows: int) -> int:
    """How many peak bands fit before the picture stops being readable."""
    return max(12, 240 // max(1, rows))


def _rgba(color: str, alpha: float) -> str:
    digits = color.lstrip("#")
    red, green, blue = (int(digits[index : index + 2], 16) for index in (0, 2, 4))
    return f"rgba({red},{green},{blue},{alpha})"


def _gene_lanes(genes: pd.DataFrame, lanes: int = 2) -> list[tuple[dict, int]]:
    """Pack genes into a couple of vertical lanes so neighbouring bodies don't merge.

    Greedy by start: a gene takes the first lane whose last gene has ended. A
    gene that finds no free lane shares the last one, which reads as one dense
    band instead of a missing gene.
    """
    packed: list[tuple[dict, int]] = []
    ends = [-1] * lanes
    for row in genes.sort_values("start").itertuples(index=False):
        first, last = int(row.start), int(row.end)
        lane = next((index for index, end in enumerate(ends) if end <= first), lanes - 1)
        ends[lane] = max(ends[lane], last)
        packed.append(
            (
                {
                    "start": first,
                    "end": last,
                    "name": str(row.name),
                    "strand": str(row.strand),
                    "gid": str(row.gid),
                },
                lane,
            )
        )
    return packed


def _draw_gene_row(
    figure,
    genes: pd.DataFrame,
    exons: pd.DataFrame | None,
    *,
    start: int,
    end: int,
    divisor: float,
    row: int,
) -> None:
    """Gene bodies, exon blocks, strand tips and names on the annotation row.

    ``start``/``end`` is the window the row is read at, not the span the bins
    cover: the label threshold and the clipping are relative to what is on
    screen, so a 10 kb view labels a 5 kb gene even though the data behind it
    reaches the chromosome ends.
    """
    span = max(1, end - start)
    packed = _gene_lanes(genes)
    lanes = max((lane for _, lane in packed), default=0) + 1
    step = 1.0 / (lanes + 1)
    blocks: dict[str, list[tuple[int, int]]] = {}
    if exons is not None and not exons.empty:
        for part in exons.itertuples(index=False):
            blocks.setdefault(str(part.gid), []).append((int(part.start), int(part.end)))
    tips: dict[str, list[float]] = {"+": [], "-": []}
    for gene, lane in packed:
        color = GENE_MINUS if gene["strand"] == "-" else GENE_PLUS
        base = 0.62 - lane * step
        left, right = max(gene["start"], start), min(gene["end"], end)
        figure.add_shape(
            type="line",
            x0=left / divisor,
            x1=right / divisor,
            y0=base,
            y1=base,
            line=dict(color=color, width=1.6),
            row=row,
            col=1,
        )
        for block_start, block_end in blocks.get(gene["gid"], [])[:40]:
            figure.add_shape(
                type="rect",
                x0=max(block_start, start) / divisor,
                x1=min(block_end, end) / divisor,
                y0=base - step * 0.3,
                y1=base + step * 0.3,
                fillcolor=color,
                line_width=0,
                row=row,
                col=1,
            )
        negative = gene["strand"] == "-"
        tips[gene["strand"]].append((left if negative else right) / divisor)
        if right - left > span * 0.012:
            figure.add_annotation(
                x=left / divisor,
                y=base + step * 0.32,
                text=gene["name"],
                showarrow=False,
                xanchor="left",
                yanchor="bottom",
                font=dict(size=9, color=color),
                row=row,
                col=1,
            )
    for strand, points in tips.items():
        if not points:
            continue
        negative = strand == "-"
        figure.add_trace(
            go.Scatter(
                x=points,
                y=[0.62 - lane * step for gene, lane in packed if gene["strand"] == strand],
                mode="markers",
                hoverinfo="skip",
                showlegend=False,
                marker=dict(
                    symbol="triangle-left" if negative else "triangle-right",
                    size=6,
                    color=GENE_MINUS if negative else GENE_PLUS,
                ),
            ),
            row=row,
            col=1,
        )


def build_chart(
    frame: pd.DataFrame,
    *,
    chrom: str,
    start: int,
    end: int,
    order: Sequence[str],
    peaks: Sequence[dict] = (),
    focus: Sequence[int] | None = None,
    genes: pd.DataFrame | None = None,
    exons: pd.DataFrame | None = None,
    unified: bool = True,
    row_height: int = 74,
    view: Sequence[int] | None = None,
    uid: str = "",
):
    """One subplot per track, with the gene model underneath when given.

    Drawn with plotly because the tab is read by pointing at it: wheel zoom,
    drag to pan, one spike line across every row and a single hover label per
    locus are built in, where altair would only render a still picture. Zoom
    magnifies the bins already loaded; a different window is typed above and
    re-read from the bigWig files, which is what a genome browser does too.

    ``start``/``end`` is the span the bins cover, ``view`` the part on screen:
    they differ once a locus is opened inside a whole-chromosome curve. ``uid``
    is plotly's uirevision: a rerun carrying the same value keeps whatever zoom
    the reader dragged, a new value (a freshly requested window) draws the new
    window as asked.

    ``peaks`` shade called peaks, ``focus`` shades the interval that was asked
    for, ``genes``/``exons`` fill the annotation row.
    """
    view_from, view_to = (int(view[0]), int(view[1])) if view else (start, end)
    span = max(1, view_to - view_from)
    unit, divisor = unit_for(span)
    colors = colors_for(order)
    present = [name for name in order if not frame.loc[frame["track"] == name].empty]
    annotate = genes is not None and not genes.empty
    nrows = len(present) + (1 if annotate else 0)
    if not nrows:
        return go.Figure()

    # Total height is solved backwards from a fixed row height: handing plotly a
    # total instead makes every row shrink as tracks are added, and with four
    # tracks the row names end up on top of the curves.
    spacing = 0.06
    weights = [1.0] * len(present) + ([1.45] if annotate else [])
    figure = make_subplots(
        rows=nrows,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=spacing,
        row_heights=weights,
    )
    margin_top, margin_bottom, margin_left = 44, 36, 62
    plot_height = int(round(row_height * sum(weights) / (1 - spacing * (nrows - 1))))

    for index, text in enumerate(
        present + ([i18n.t("track.gene_model")] if annotate else []), start=1
    ):
        # Track name in the gap above its row, left aligned. yref="y domain" is
        # the row's own 0-1 box, so 1.0 is that row's top edge, not the figure's.
        figure.add_annotation(
            xref="x domain",
            yref="y domain",
            x=0.0,
            y=1.0,
            xanchor="left",
            yanchor="bottom",
            text=f"<b>{text}</b>",
            showarrow=False,
            font=dict(size=11, color="#333"),
            bgcolor="rgba(255,255,255,0.85)",
            borderpad=2,
            row=index,
            col=1,
        )

    if focus:
        focus_from = max(int(focus[0]), start)
        focus_to = min(int(focus[1]), end)
        figure.add_vrect(
            x0=focus_from / divisor,
            x1=focus_to / divisor,
            fillcolor=_rgba(FOCUS_BAND, 0.18),
            line_width=0,
            layer="below",
            row="all",
            col=1,
            exclude_empty_subplots=False,
        )
        if (focus_to - focus_from) < span * 0.004:
            # A 10 kb locus on a whole chromosome is a fraction of a pixel, so
            # the band alone would leave the reader hunting for where they are.
            figure.add_vline(
                x=(focus_from + focus_to) / 2 / divisor,
                line=dict(color=FOCUS_BAND, width=1.2, dash="dot"),
                row="all",
                col=1,
                exclude_empty_subplots=False,
            )

    # Peaks shade the signal rows only: the annotation row has its own marks and
    # a band under them hides the exons. Shapes are per row, so the cap is
    # shared out; peaks arrive sorted by score, so a cap drops the weakest.
    peak_cap = peak_draw_cap(len(present))
    for row_index in range(1, len(present) + 1):
        for peak in list(peaks)[:peak_cap]:
            figure.add_vrect(
                x0=max(int(peak["start"]), start) / divisor,
                x1=min(int(peak["end"]), end) / divisor,
                fillcolor=_rgba(PEAK_BAND, 0.26),
                line_width=0,
                layer="below",
                row=row_index,
                col=1,
                exclude_empty_subplots=False,
            )

    shared = shared_top(frame)
    for index, name in enumerate(present, start=1):
        values = frame[frame["track"] == name]
        positions = values["position"].to_numpy(dtype=float)
        signal = values["signal"].to_numpy(dtype=float)
        figure.add_trace(
            go.Scatter(
                x=positions / divisor,
                y=signal,
                name=name,
                mode="lines",
                line=dict(color=colors[name], width=1.0),
                fill="tozeroy",
                fillcolor=_rgba(colors[name], 0.30),
                hovertemplate=f"%{{y:.3f}}<extra>{name}</extra>",
            ),
            row=index,
            col=1,
        )
        top = shared if unified else nice_max(float(signal.max()))
        figure.update_yaxes(
            range=[0, top * 1.02],
            row=index,
            col=1,
            showgrid=True,
            gridcolor=GRID,
            gridwidth=0.5,
            tickfont=dict(size=9),
            nticks=3,
            # Vertical zoom would pull the rows apart, and comparing them is the
            # point of the stack; dragging therefore only spans the x axis.
            fixedrange=True,
        )

    if annotate:
        _draw_gene_row(
            figure, genes, exons, start=view_from, end=view_to, divisor=divisor, row=nrows
        )
        figure.update_yaxes(range=[0, 1], visible=False, row=nrows, col=1)
        figure.update_xaxes(showline=True, linecolor="#bbb", row=nrows, col=1)

    for index in range(1, nrows + 1):
        figure.update_xaxes(
            range=[view_from / divisor, view_to / divisor],
            row=index,
            col=1,
            showgrid=index == nrows,
            gridcolor=GRID,
            color="#666",
            tickfont=dict(size=9),
            # A whole chromosome reads in Mb, where one hover digit hides the
            # difference between two neighbouring bins.
            hoverformat=",.3f" if divisor >= 1_000_000 else ",.1f",
            # A spike through every row on hover is what makes a peak line up
            # with the gene model at the bottom without counting pixels.
            showspikes=True,
            spikemode="across",
            spikesnap="cursor",
            spikethickness=1,
            spikecolor="#999",
            spikedash="solid",
            showline=True,
            linecolor="#bbb",
        )
    figure.update_xaxes(
        title_text=f"{chrom} ({unit})", title_font=dict(size=10), row=nrows, col=1
    )
    figure.update_layout(
        height=plot_height + margin_top + margin_bottom,
        margin=dict(l=margin_left, r=16, t=margin_top, b=margin_bottom),
        hovermode="x unified",
        dragmode="zoom",
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=False,
        font=dict(family='system-ui, "Noto Sans SC", sans-serif', size=10),
        hoverlabel=dict(
            bgcolor="rgba(255,255,255,0.94)", bordercolor="#ccc", font=dict(size=10)
        ),
        **({"uirevision": uid} if uid else {}),
    )
    return figure


def window_stats(frame: pd.DataFrame, order: Sequence[str]) -> pd.DataFrame:
    """Per-track summary of the current window, computed on raw values."""
    records = []
    for name in order:
        values = frame.loc[frame["track"] == name, "signal"]
        if values.empty:
            continue
        records.append(
            {
                "track": name,
                "mean": round(float(values.mean()), 3),
                "max": round(float(values.max()), 3),
                "covered": f"{float((values > 0).mean()) * 100:.0f}%",
            }
        )
    return pd.DataFrame(records, columns=["track", "mean", "max", "covered"])
