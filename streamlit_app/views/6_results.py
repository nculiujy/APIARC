from __future__ import annotations

import os
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

_APP_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _APP_DIR not in sys.path:
    sys.path.insert(0, _APP_DIR)

from utils import datasets, env_paths, i18n, pipeline, runtime, tracks  # noqa: E402
from utils.ui_theme import (  # noqa: E402
    empty_state,
    metrics,
    page_header,
    section_title,
    stage_rows,
)

config = runtime.load_config()
rna_name = datasets.config_dataset(config, "RNAseq")
chip_name = datasets.config_dataset(config, "ChIPseq")
joint = datasets.joint_id(config)
annotation = datasets.annotation_info(config)

SCOPE_KEY = {
    "RNAseq": "res.scope.rna",
    "ChIPseq": "res.scope.chip",
    "Integrated": "res.scope.joint",
}

page_header(
    i18n.t("res.title"),
    chips=[
        (i18n.t("res.chip.rna", name=rna_name), "ok"),
        (i18n.t("res.chip.chip", name=chip_name), "ok"),
        (i18n.t("res.chip.annotation", id=annotation["id"]), "ok"),
        (i18n.t("res.chip.joint", id=joint), "ok"),
    ],
)


def tree_size(path: Path) -> int:
    return _tree_size(str(path))


# The overview asks for the size of the same dozen directories on every rerun,
# and walking them took most of the page's load time.
@st.cache_data(ttl=30, show_spinner=False)
def _tree_size(key: str) -> int:
    path = Path(key)
    if not path.exists():
        return 0
    if path.is_file():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += (Path(root) / name).stat().st_size
            except OSError:
                continue
    return total


PREVIEW_ROWS = 100
PLOT_POINTS = 6000

# ChIPseeker writes the class of a peak in `annotation_category` and a per-peak
# description in `annotation`; the descriptions are nearly all distinct, so a
# count over them gives one bar of height one per peak.
FEATURE_COLUMNS = ("annotation_category", "feature", "region", "class")
DETAIL_COLUMNS = ("annotation",)

FEATURE_KEYS = {
    "intron": "res.feature.intron",
    "exon": "res.feature.exon",
    "intergenic": "res.feature.intergenic",
    "promoter": "res.feature.promoter",
    "5' utr": "res.feature.utr5",
    "3' utr": "res.feature.utr3",
}


def feature_values(frame: pd.DataFrame) -> pd.Series | None:
    """The genomic feature of every peak, or None when the table has no such column."""
    for wanted in FEATURE_COLUMNS:
        match = next((name for name in frame.columns if str(name).lower() == wanted), None)
        if match is not None:
            return frame[match].astype(str).str.strip()
    for wanted in DETAIL_COLUMNS:
        match = next((name for name in frame.columns if str(name).lower() == wanted), None)
        if match is not None:
            text = frame[match].astype(str).str.strip()
            # The description carries the class first, then the transcript in
            # brackets: "Intron (ENSMUST..., intron 3 of 5)".
            return text.str.split(" (", n=1).str[0].str.strip()
    return None


def feature_label(value: str) -> str:
    """The feature name in the reader's language; an unknown one stays as is."""
    key = FEATURE_KEYS.get(str(value).strip().lower())
    if not key:
        return str(value)
    text = i18n.t(key)
    return str(value) if text == key else text


def read_table(path: Path, **kwargs) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    try:
        if path.suffix == ".tsv":
            return pd.read_csv(path, sep="\t", **kwargs)
        if path.suffix == ".csv":
            return pd.read_csv(path, **kwargs)
        if path.suffix in {".bed", ".txt"}:
            return pd.read_csv(path, sep="\t", header=None, **kwargs)
    except Exception:  # noqa: BLE001
        return None
    return None


def download(path: Path, label: str | None = None) -> None:
    if not path.is_file():
        return
    size = i18n.size(path.stat().st_size)
    st.download_button(
        label or path.name,
        # Deferred read: the page reruns on every widget change, and pulling a
        # multi-megabyte result into memory each time stalls the whole tab.
        lambda: path.read_bytes(),
        file_name=path.name,
        key=f"dl-{path}",
        help=f"{env_paths.relative(path)} · {size}",
    )


def preview(frame: pd.DataFrame) -> None:
    """Show only the head of a table; the file keeps every row.

    A differential table runs to tens of thousands of rows and st.dataframe
    ships all of them to the browser, which is what made this page crawl.
    """
    total = len(frame)
    st.dataframe(frame.head(PREVIEW_ROWS), width="stretch", hide_index=True)
    if total > PREVIEW_ROWS:
        st.caption(
            i18n.t(
                "res.preview.more",
                rows=i18n.number(PREVIEW_ROWS),
                total=i18n.number(total),
            )
        )


def pdf_gallery(paths: list[Path], base: Path, key: str, label: str) -> None:
    """Render one PDF at a time.

    An expander still runs its body on every rerun, so rasterising a dozen
    pages stalled the tab before anything below it was drawn. A picker keeps
    the work down to the file on screen.
    """
    if not paths:
        st.caption(i18n.t("res.no_plots"))
        return
    names = [str(path.relative_to(base)) for path in paths]
    picked = st.selectbox(label, names, key=key)
    pdf_preview(paths[names.index(picked)])


def pdf_preview(path: Path) -> None:
    if not path.is_file():
        return
    st.caption(env_paths.relative(path))
    try:
        import fitz  # type: ignore

        with fitz.open(path) as document:
            page = document.load_page(0)
            pixmap = page.get_pixmap(dpi=110)
            st.image(pixmap.tobytes("png"))
    except Exception:  # noqa: BLE001
        st.info(i18n.t("res.pymupdf"))
    download(path)


def stage_overview(scope: str) -> None:
    statuses = pipeline.stage_status(scope, config, markers(scope))
    stages = datasets.STAGE_DIRS[scope] if scope != "Integrated" else datasets.INTEGRATED_STAGES
    rows = [
        (label, f"{stages[key]} · {i18n.size(tree_size(stage_dir(scope, key)))}", status, kind)
        for (label, _detail, status, kind), key in zip(statuses, stages)
    ]
    stage_rows(rows)


section_title(i18n.t("res.section.scope"))


def source_picker(omics: str, name: str, column, label: str) -> tuple[str, str]:
    """Dataset and annotation this page reads for one omics.

    A switched reference or a mock run leaves its results behind under another
    ``<dataset>/<annotation>`` directory; defaulting to the configured pair alone
    would show four empty stages while those files sit on disk.
    """
    options = datasets.result_inventory(omics)
    if not options:
        column.caption(i18n.t("res.pick.none", omics=omics))
        return name, annotation["id"]
    labels = [f"{dataset} · {version}" for dataset, version in options]
    configured = f"{name} · {annotation['id']}"
    default = labels.index(configured) if configured in labels else 0
    picked = column.selectbox(label, labels, index=default, key=f"view-{omics}")
    return options[labels.index(picked)]


left, right = st.columns(2)
rna_dataset, rna_version = source_picker(
    "RNAseq", rna_name, left, i18n.t("res.pick.rna")
)
chip_dataset, chip_version = source_picker(
    "ChIPseq", chip_name, right, i18n.t("res.pick.chip")
)

if (rna_dataset, rna_version) != (rna_name, annotation["id"]) or (
    chip_dataset,
    chip_version,
) != (chip_name, annotation["id"]):
    st.caption(
        i18n.t(
            "res.pick.configured",
            name=rna_name,
            version=annotation["id"],
        )
    )

same_version = rna_version == chip_version
if not same_version:
    st.warning(i18n.t("res.pick.mismatch", rna=rna_version, chip=chip_version))

rna_root = datasets.result_directory("RNAseq", rna_dataset, annotation=rna_version)
chip_root = datasets.result_directory("ChIPseq", chip_dataset, annotation=chip_version)
joint_root = datasets.joint_directory(
    config, rna=rna_dataset, chip=chip_dataset, annotation=rna_version
)


def rna_stage(key: str) -> Path:
    return datasets.stage("RNAseq", rna_dataset, key, annotation=rna_version)


def chip_stage(key: str) -> Path:
    return datasets.stage("ChIPseq", chip_dataset, key, annotation=chip_version)


def joint_stage(key: str) -> Path:
    return datasets.joint_stage(
        config, key, rna=rna_dataset, chip=chip_dataset, annotation=rna_version
    )


def stage_dir(scope: str, key: str) -> Path:
    if scope == "Integrated":
        return joint_stage(key)
    return rna_stage(key) if scope == "RNAseq" else chip_stage(key)


def markers(scope: str) -> dict[str, Path]:
    keys = datasets.INTEGRATED_STAGES if scope == "Integrated" else datasets.STAGE_DIRS[scope]
    return {key: stage_dir(scope, key) / datasets.FINISH for key in keys}


section_title(i18n.t("res.section.overview"))
rna_done, rna_total = pipeline.progress(
    "RNAseq", config, markers("RNAseq"), only=pipeline.active_stages("RNAseq", config)
)
chip_done, chip_total = pipeline.progress(
    "ChIPseq", config, markers("ChIPseq"), only=pipeline.active_stages("ChIPseq", config)
)
if same_version:
    # Only the stages the Snakefile actually compiled in: counting all five would
    # leave the total unreachable on a descriptive design, where 4_common is not
    # included at all.
    joint_done, joint_total = pipeline.progress(
        "Integrated", config, markers("Integrated"),
        only=pipeline.active_stages("Integrated", config),
    )
else:
    joint_done, joint_total = 0, len(datasets.INTEGRATED_STAGES)
metrics(
    [
        (
            i18n.t("res.metric.rna"),
            i18n.ratio(rna_done, rna_total),
            i18n.size(tree_size(rna_root)),
        ),
        (
            i18n.t("res.metric.chip"),
            i18n.ratio(chip_done, chip_total),
            i18n.size(tree_size(chip_root)),
        ),
        (
            i18n.t("res.metric.joint"),
            i18n.ratio(joint_done, joint_total),
            i18n.size(tree_size(joint_root)),
        ),
        (
            i18n.t("res.metric.total"),
            i18n.size(tree_size(rna_root) + tree_size(chip_root) + tree_size(joint_root)),
            f"{rna_version} / {chip_version}",
        ),
    ]
)

with st.expander(i18n.t("res.stage_products"), expanded=False):
    scope = st.radio(
        i18n.t("res.scope.label"),
        ["RNAseq", "ChIPseq", "Integrated"],
        format_func=lambda value: i18n.t(SCOPE_KEY[value]),
        horizontal=True,
        key="overview-scope",
    )
    stage_overview(scope)

tab_rna, tab_chip, tab_joint, tab_track = st.tabs(
    [
        "RNA-seq",
        "ChIP-seq",
        i18n.t("res.tab.joint"),
        i18n.t("res.tab.track"),
    ]
)

with tab_rna:
    section_title(i18n.t("res.rna.title"), i18n.t("res.rna.hint"))
    deg_path = rna_stage("4_DEseq") / "DEG_result.csv"
    frame = read_table(deg_path)
    if frame is None or frame.empty:
        empty_state(i18n.t("res.rna.empty"))
    else:
        frame = frame.dropna(subset=["padj", "log2FoldChange"], how="any")
        columns = st.columns([1, 1, 2])
        alpha = columns[0].slider(
            i18n.t("res.rna.alpha"),
            0.001,
            0.2,
            float(config.get("stats", {}).get("alpha", 0.05)),
            0.005,
        )
        lfc_cut = columns[1].slider(i18n.t("res.rna.lfc"), 0.0, 4.0, 1.0, 0.1)
        significant = frame[(frame["padj"] <= alpha) & (frame["log2FoldChange"].abs() >= lfc_cut)]
        up = int((significant["log2FoldChange"] > 0).sum())
        down = int((significant["log2FoldChange"] < 0).sum())
        metrics(
            [
                (
                    i18n.t("res.rna.metric.genes"),
                    i18n.number(len(frame)),
                    i18n.t("res.rna.metric.genes_hint"),
                ),
                (
                    i18n.t("res.rna.metric.significant"),
                    i18n.number(len(significant)),
                    i18n.t(
                        "res.rna.metric.significant_hint",
                        alpha=i18n.number(alpha, 3),
                        lfc=i18n.number(lfc_cut, 1),
                    ),
                ),
                (
                    i18n.t("res.rna.metric.up"),
                    i18n.number(up),
                    i18n.t("res.rna.metric.up_hint"),
                ),
                (
                    i18n.t("res.rna.metric.down"),
                    i18n.number(down),
                    i18n.t("res.rna.metric.down_hint"),
                ),
            ]
        )

        # Same look as the pipeline's R volcano: small low-alpha points, blue for
        # down and red for up, dashed cutoffs, counts in the title. Zooming a
        # scatter of thousands of genes only ever hid the shape of the cloud.
        plot_frame = frame.copy()
        plot_frame["neglog10"] = -np.log10(plot_frame["padj"].clip(lower=1e-300))
        plot_frame["label"] = "ns"
        plot_frame.loc[(plot_frame["padj"] <= alpha) & (plot_frame["log2FoldChange"] >= lfc_cut), "label"] = "up"
        plot_frame.loc[(plot_frame["padj"] <= alpha) & (plot_frame["log2FoldChange"] <= -lfc_cut), "label"] = "down"

        # Every rerun re-sends the chart payload; thousands of grey points were
        # enough to stall the tab. Signal points all stay, the background is
        # sampled, and the counts above still come from the whole frame.
        if len(plot_frame) > PLOT_POINTS:
            signal = plot_frame[plot_frame["label"] != "ns"]
            background = plot_frame[plot_frame["label"] == "ns"]
            room = max(0, PLOT_POINTS - len(signal))
            if len(background) > room:
                background = background.sample(n=room, random_state=0)
            plot_frame = pd.concat([signal, background])
            st.caption(
                i18n.t(
                    "res.rna.sampled",
                    points=i18n.number(len(plot_frame)),
                    genes=i18n.number(len(frame)),
                )
            )

        color = alt.Color(
            "label:N",
            scale=alt.Scale(domain=["down", "ns", "up"], range=["#2b6cb0", "#b8c2cc", "#c53030"]),
            legend=alt.Legend(title=None, orient="right"),
        )
        # One legend for the pair: the second chart encodes the same three
        # classes and repeated it verbatim beside the first.
        color_plain = alt.Color(
            "label:N",
            scale=alt.Scale(domain=["down", "ns", "up"], range=["#2b6cb0", "#b8c2cc", "#c53030"]),
            legend=None,
        )
        thresholds = alt.Chart(
            pd.DataFrame({"x": [-lfc_cut, lfc_cut], "y": [-np.log10(alpha), -np.log10(alpha)]})
        )
        volcano = (
            alt.layer(
                alt.Chart(plot_frame)
                .mark_circle(size=7, opacity=0.4)
                .encode(
                    x=alt.X("log2FoldChange:Q", title=i18n.t("res.axis.fold")),
                    y=alt.Y("neglog10:Q", title=i18n.t("res.axis.padj")),
                    color=color,
                ),
                thresholds.mark_rule(color="#999", strokeDash=[4, 4]).encode(x="x:Q"),
                thresholds.mark_rule(color="#999", strokeDash=[4, 4]).encode(y="y:Q"),
            )
            .properties(height=340)
        )
        ma = (
            alt.layer(
                alt.Chart(plot_frame)
                .mark_circle(size=7, opacity=0.4)
                .encode(
                    x=alt.X(
                        "baseMean:Q",
                        scale=alt.Scale(type="log", domainMin=0.5),
                        title=i18n.t("res.axis.basemean"),
                    ),
                    y=alt.Y("log2FoldChange:Q", title=i18n.t("res.axis.fold")),
                    color=color_plain,
                ),
                alt.Chart(pd.DataFrame({"y": [0.0]})).mark_rule(color="#999", strokeDash=[4, 4]).encode(y="y:Q"),
            )
            .properties(height=340)
        )
        # Side by side: stacked at full width the second chart pushed the gene
        # table another screen down for the same points.
        plot_left, plot_right = st.columns(2)
        for column, chart in ((plot_left, volcano), (plot_right, ma)):
            column.altair_chart(
                chart.configure_view(stroke=None).configure_axis(
                    grid=False, domainColor="#888", tickColor="#888", labelFontSize=10
                ),
                width="stretch",
            )

        section_title(i18n.t("res.rna.table"), i18n.t("res.rna.table_hint"))
        show = significant.sort_values("padj")
        preview(
            show.assign(
                log2FoldChange=show["log2FoldChange"].round(3),
                pvalue=show["pvalue"].map(lambda value: f"{value:.3e}"),
                padj=show["padj"].map(lambda value: f"{value:.3e}"),
                baseMean=show["baseMean"].round(1),
            )
        )
        download(deg_path)

        tables = rna_stage("4_DEseq") / "tables"
        plots = rna_stage("4_DEseq") / "plots"
        left, right = st.columns(2)
        with left:
            section_title(i18n.t("res.rna.tables"))
            for path in sorted(path for path in tables.rglob("*") if path.is_file()):
                download(path, str(path.relative_to(tables)))
        with right:
            section_title(i18n.t("res.rna.plots"))
            pdf_gallery(sorted(plots.rglob("*.pdf")), plots, "rna-plot", i18n.t("res.plot"))

with tab_chip:
    section_title(i18n.t("res.chip.title"), i18n.t("res.chip.hint"))
    peak_root = chip_stage("4_peak_result")
    summary = read_table(peak_root / "peakcalling" / "peak_calling_summary.tsv")
    if summary is None or summary.empty:
        empty_state(i18n.t("res.chip.empty"))
    else:
        preview(summary)
        download(peak_root / "peakcalling" / "peak_calling_summary.tsv")

        beds = sorted((peak_root / "peakcalling").glob("*/*_peaks.bed"))
        metrics(
            [
                (
                    i18n.t("res.chip.metric.samples"),
                    i18n.number(len(beds)),
                    i18n.t("res.chip.metric.samples_hint"),
                ),
                (
                    i18n.t("res.chip.metric.anno"),
                    i18n.number(len(list((peak_root / "annotation").glob("*/*_peak_anno.csv")))),
                    "peak_anno.csv",
                ),
                (
                    i18n.t("res.chip.metric.bigwig"),
                    i18n.number(len(list((peak_root / "deeptools").glob("*/*.bw")))),
                    i18n.t("res.chip.metric.bigwig_hint"),
                ),
                (
                    i18n.t("res.chip.metric.merged"),
                    i18n.number(
                        len(list((peak_root / "deeptools").glob("*/*_peaks_merged.bed")))
                    ),
                    i18n.t("res.chip.metric.merged_hint"),
                ),
            ]
        )

        anno_files = sorted((peak_root / "annotation").glob("*/*_peak_anno.csv"))
        if anno_files:
            section_title(i18n.t("res.chip.anno.title"), i18n.t("res.chip.anno.hint"))
            names = [path.parent.name for path in anno_files]
            picked = st.selectbox(i18n.t("res.chip.anno.sample"), names, key="anno-sample")
            peak_anno = read_table(anno_files[names.index(picked)])
            if peak_anno is not None and not peak_anno.empty:
                features = feature_values(peak_anno)
                if features is not None:
                    feature_col = i18n.t("res.chip.anno.feature")
                    count_col = i18n.t("res.chip.anno.count")
                    counts = (
                        features.map(feature_label)
                        .value_counts()
                        .rename_axis(feature_col)
                        .reset_index(name=count_col)
                    )
                    st.altair_chart(
                        alt.Chart(counts)
                        .mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4)
                        .encode(
                            x=alt.X(field=feature_col, type="nominal", sort="-y", title=None),
                            y=alt.Y(field=count_col, type="quantitative", title=i18n.t("res.chip.anno.y")),
                            color=alt.value("#1f6feb"),
                            tooltip=[
                                alt.Tooltip(field=feature_col, type="nominal"),
                                alt.Tooltip(field=count_col, type="quantitative"),
                            ],
                        )
                        .properties(height=300),
                        width="stretch",
                    )
                preview(peak_anno)
                download(anno_files[names.index(picked)])

        bws = sorted((peak_root / "deeptools").glob("*/*.bw"))
        if bws:
            section_title(i18n.t("res.chip.signal.title"))
            st.dataframe(
                pd.DataFrame(
                    [
                        {
                            i18n.t("res.column.track"): path.parent.name,
                            i18n.t("res.column.file"): path.name,
                            i18n.t("res.column.size"): i18n.size(path.stat().st_size),
                        }
                        for path in bws
                    ]
                ),
                width="stretch",
                hide_index=True,
            )

        pdf_paths = sorted((peak_root / "deeptools").glob("*/*.pdf"))
        if pdf_paths:
            section_title(i18n.t("res.chip.qc.title"))
            pdf_gallery(pdf_paths, peak_root / "deeptools", "chip-plot", i18n.t("res.plot"))

with tab_joint:
    section_title(i18n.t("res.joint.title"), i18n.t("res.joint.hint"))
    stats_path = joint_stage("0_stats") / "stats_summary.tsv"
    stats = read_table(stats_path)
    if not same_version:
        st.info(i18n.t("res.joint.mismatch", rna=rna_version, chip=chip_version))
    elif stats is None or stats.empty:
        empty_state(i18n.t("res.joint.empty"))
    else:
        preview(stats)
        download(stats_path)

        for label, filename in (
            (i18n.t("res.joint.enrichment"), "enrichment_summary.tsv"),
            (i18n.t("res.joint.motif"), "motif_summary.tsv"),
            (i18n.t("res.joint.enhancer"), "enhancer_summary.tsv"),
        ):
            stage_key = {"enrichment_summary.tsv": "1_enrichment", "motif_summary.tsv": "2_motif", "enhancer_summary.tsv": "3_enhancer"}[filename]
            path = joint_stage(stage_key) / filename
            table = read_table(path)
            if table is None or table.empty:
                continue
            with st.expander(f"{label} · {filename}", expanded=False):
                preview(table)
                download(path)

        motif_edges = read_table(joint_stage("2_motif") / "tf_gene_edges.tsv")
        enhancer_links = read_table(joint_stage("3_enhancer") / "enhancer_gene_links.tsv")
        common_tables = joint_stage("4_common") / "tables"
        genes = read_table(common_tables / "common_genes.tsv")
        metrics(
            [
                (
                    i18n.t("res.joint.metric.regions"),
                    i18n.number(len(stats)),
                    i18n.t("res.joint.metric.regions_hint"),
                ),
                (
                    i18n.t("res.joint.metric.genes"),
                    i18n.number(0 if genes is None else len(genes)),
                    i18n.t("res.joint.metric.genes_hint"),
                ),
                (
                    i18n.t("res.joint.metric.edges"),
                    i18n.number(0 if motif_edges is None else len(motif_edges)),
                    i18n.t("res.joint.metric.edges_hint"),
                ),
                (
                    i18n.t("res.joint.metric.links"),
                    i18n.number(0 if enhancer_links is None else len(enhancer_links)),
                    i18n.t("res.joint.metric.links_hint"),
                ),
                (
                    i18n.t("res.joint.metric.stages"),
                    i18n.ratio(joint_done, joint_total),
                    i18n.t("res.joint.metric.stages_hint"),
                ),
            ]
        )

        section_title(i18n.t("res.joint.common.title"), i18n.t("res.joint.common.hint"))
        if genes is None or genes.empty:
            st.caption(i18n.t("res.joint.common.empty"))
        else:
            filters = st.columns(3, vertical_alignment="bottom")
            marks = filters[0].multiselect(
                i18n.t("res.joint.filter.mark"),
                sorted(genes["mark"].unique()),
                key="common-marks",
            )
            regions = filters[1].multiselect(
                i18n.t("res.joint.filter.region"),
                sorted(genes["region_type"].unique()),
                key="common-regions",
            )
            # Concordance only says whether ChIP and RNA move the same way; it
            # says nothing about effect size, so it gets its own control rather
            # than joining the significance filters, and defaults to "all" so
            # genes moving in opposite directions stay visible.
            directions = [
                i18n.t("res.joint.direction.all"),
                i18n.t("res.joint.direction.same"),
                i18n.t("res.joint.direction.opposite"),
            ]
            concordance = filters[2].selectbox(
                i18n.t("res.joint.filter.direction"), directions, key="common-direction"
            )
            picked = genes
            if marks:
                picked = picked[picked["mark"].isin(marks)]
            if regions:
                picked = picked[picked["region_type"].isin(regions)]
            same = picked["concordant"].astype(str).str.lower() == "true"
            if concordance == directions[1]:
                picked, same = picked[same], same[same]
            elif concordance == directions[2]:
                picked, same = picked[~same], same[~same]
            metrics(
                [
                    (
                        i18n.t("res.joint.metric.rows"),
                        i18n.number(len(picked)),
                        i18n.t("res.joint.metric.rows_hint"),
                    ),
                    (
                        i18n.t("res.joint.metric.distinct"),
                        i18n.number(picked["gene"].nunique()),
                        i18n.t("res.joint.metric.distinct_hint"),
                    ),
                    (
                        i18n.t("res.joint.metric.same"),
                        i18n.number(int(same.sum())),
                        i18n.t("res.joint.metric.same_hint"),
                    ),
                    (
                        i18n.t("res.joint.metric.opposite"),
                        i18n.number(len(picked) - int(same.sum())),
                        i18n.t("res.joint.metric.opposite_hint"),
                    ),
                ]
            )
            wanted = [
                "gene", "mark", "comparison", "region_type", "chip_direction", "rna_direction",
                "n_regions", "chip_lfc_mean", "chip_fdr", "rna_log2FC", "rna_padj",
            ]
            show = picked.sort_values(["mark", "comparison", "region_type", "chip_fdr"])
            show = show[[name for name in wanted if name in show.columns]]
            numeric = {
                "chip_lfc_mean": lambda value: round(float(value), 3),
                "rna_log2FC": lambda value: round(float(value), 3),
                "chip_fdr": lambda value: f"{float(value):.3e}",
                "rna_padj": lambda value: f"{float(value):.3e}",
            }
            preview(
                show.assign(
                    **{name: show[name].map(func) for name, func in numeric.items() if name in show.columns}
                )
            )
            download(common_tables / "common_genes.tsv")
            download(common_tables / "common_stats.tsv")
            with st.expander(
                f"{i18n.t('res.joint.enrich_title')} · go_kegg_common.tsv", expanded=False
            ):
                # Only the columns the panel shows: the whole table carries the
                # gene lists and runs to 6 MB.
                enrichment = read_table(
                    common_tables / "go_kegg_common.tsv",
                    usecols=lambda name: name
                    in {
                        "mark", "comparison", "region_type", "ontology", "term_id",
                        "term", "n_genes", "n_term", "odds_ratio", "p", "FDR",
                    },
                )
                if enrichment is None or enrichment.empty:
                    st.caption(i18n.t("res.joint.enrich_empty"))
                else:
                    st.caption(i18n.t("res.joint.enrich_note"))
                    preview(enrichment.sort_values("FDR"))
                    download(common_tables / "go_kegg_common.tsv")
                    download(common_tables / "go_kegg_all_diff.tsv")
                    download(common_tables / "network_edges.tsv")
                    download(common_tables / "network_nodes.tsv")

        section_title(i18n.t("res.joint.plots"))
        joint_pdfs = sorted(
            path
            for stage_key in pipeline.JOINT_STAGES
            for path in joint_stage(stage_key).glob("**/*.pdf")
        )
        pdf_gallery(joint_pdfs, joint_root, "joint-plot", i18n.t("res.plot"))

with tab_track:
    section_title(i18n.t("res.track.title"), i18n.t("res.track.hint"))

    def reset_jump() -> None:
        """Called before the run starts, so the overview is drawn with no jump."""
        st.session_state["bw-jump"] = 0

    peak_root = chip_stage("4_peak_result")
    bws = sorted((peak_root / "deeptools").glob("*/*.bw"))
    if not bws:
        empty_state(i18n.t("res.track.empty"))
    elif not tracks.have_pybigwig():
        st.warning(i18n.t("res.track.no_pybigwig"))
    elif not tracks.have_plotly():
        st.warning(i18n.t("res.track.no_plotly"))
    else:
        # The label carries the group name: one Input library shows up under
        # every group directory, and the bare file name would not say which
        # comparison it belongs to, or two tables would hold rows with the same
        # name.
        sources = {f"{path.parent.name}/{path.stem}": path for path in bws}
        # The group name is sample_name in metadata.csv; deeptools splits its
        # output into one directory per group.
        groups: dict[str, list[str]] = {}
        for label, path in sources.items():
            groups.setdefault(path.parent.name, []).append(label)
        for labels in groups.values():
            # Input goes last within a group: the default track is always an IP.
            labels.sort(key=lambda label: sources[label].stem.startswith("Input_"))

        peak_files = sorted((peak_root / "peakcalling").glob("*/*_peaks.narrowPeak")) or sorted(
            (peak_root / "deeptools").glob("*/*_peaks_merged.bed")
        )
        peaks = tracks.read_peaks(peak_files)

        sizes: dict[str, int] = {}
        for path in sources.values():
            for chrom_name, length in tracks.chrom_sizes(path).items():
                sizes[chrom_name] = max(sizes.get(chrom_name, 0), length)

        if not sizes:
            empty_state(i18n.t("res.track.no_chrom"))
        else:
            gtf = Path(datasets.annotation_info(config)["gtf"])
            if not gtf.is_absolute():
                gtf = env_paths.BASE_DIR / gtf
            genes = tracks.gene_index(gtf)
            exons = tracks.exon_index(gtf)

            if "bw-region-default" not in st.session_state:
                # Open on the strongest peak this run called, so the tab shows
                # something on the first load instead of an arbitrary window.
                anchor = next(iter(tracks.top_peaks(peaks, 1)), None)
                home = (
                    anchor["chrom"]
                    if anchor and anchor["chrom"] in sizes
                    else max(sizes, key=sizes.get)
                )
                centre = (anchor["start"] + anchor["end"]) // 2 if anchor else sizes[home] // 2
                left_edge = max(0, min(centre - 5_000, sizes[home] - 10_000))
                st.session_state["bw-region-default"] = (
                    f"{home}:{left_edge}-{min(sizes[home], left_edge + 10_000)}"
                )
                st.session_state["bw-region"] = st.session_state["bw-region-default"]
                st.session_state["bw-applied"] = st.session_state["bw-region-default"]

            # What is typed and what is drawn are two states: a half-typed
            # coordinate should neither repaint the figure nor raise, so View
            # applies it. An emptied box falls back to the opening window.
            lookup = st.columns([3.4, 1, 1.2], vertical_alignment="bottom")
            query = lookup[0].text_input(
                i18n.t("res.track.region"),
                key="bw-region",
                placeholder=i18n.t("res.track.region_placeholder"),
            )
            if lookup[1].button(
                i18n.t("res.track.look"), key="bw-look", type="primary"
            ):
                st.session_state["bw-applied"] = (
                    query.strip() or st.session_state["bw-region-default"]
                )
                # Counted so that asking for the same locus twice still brings
                # the picture back to it after the reader has zoomed away.
                st.session_state["bw-jump"] = st.session_state.get("bw-jump", 0) + 1
            if st.session_state.get("bw-jump", 0) > 0:
                lookup[2].button(
                    i18n.t("res.track.whole"), key="bw-whole", on_click=reset_jump
                )

            options = st.columns([1.1, 1, 1, 1.3], vertical_alignment="bottom")
            smooth_bp = options[0].select_slider(
                i18n.t("res.track.smooth"),
                [0, 25, 50, 100, 200, 500, 1000],
                100,
                key="bw-smooth",
                help=i18n.t("res.track.smooth_help"),
            )
            unified = options[1].checkbox(
                i18n.t("res.track.unified"),
                True,
                key="bw-shared",
                help=i18n.t("res.track.unified_help"),
            )
            show_peaks = options[2].checkbox(
                i18n.t("res.track.peaks"),
                False,
                key="bw-peaks",
                help=i18n.t("res.track.peaks_help"),
            )
            bins = options[3].slider(
                i18n.t("res.track.bins"), 100, 2000, tracks.DEFAULT_BINS, 50
            )

            try:
                locus = tracks.parse_locus(st.session_state.get("bw-applied", ""), sizes, genes)
            except ValueError as error:
                st.error(str(error))
            else:
                chrom = locus["chrom"]
                start, end = locus["start"], locus["end"]
                gene = locus["gene"]
                chrom_len = sizes.get(chrom, end)
                # The chart opens on the whole chromosome, the way IGV does;
                # View only moves the window onto the locus.
                overview = st.session_state.get("bw-jump", 0) == 0
                if gene:
                    # Whether the gene body itself is called is the point of
                    # looking at a gene name, so it belongs in the line, not
                    # only in the picture.
                    inside = tracks.peaks_in(
                        peaks, gene["chrom"], gene["start"], gene["end"], limit=999
                    )
                    parts = [
                        f"GENCODE {gene['name']}",
                        f"{gene['chrom']}:{i18n.number(gene['start'])}-{i18n.number(gene['end'])}",
                        i18n.t("res.track.strand", strand=gene["strand"]),
                        str(gene["biotype"]),
                        (
                            i18n.t("res.track.gene_inside", n=i18n.number(len(inside)))
                            if inside
                            else i18n.t("res.track.gene_outside")
                        ),
                    ]
                    st.caption(" · ".join(parts))

                table = st.session_state.get("bw-track-rows")
                if (
                    table is None
                    or "track" not in table.columns
                    or set(table["track"]) != set(sources)
                ):
                    defaults = tracks.comparison_defaults(
                        groups,
                        [
                            (row.get("reference", ""), row.get("treatment", ""))
                            for row in datasets.comparison_rows("ChIPseq", chip_name)
                        ],
                    )
                    rows_in_order = defaults + [label for label in sources if label not in defaults]
                    table = pd.DataFrame(
                        [
                            {
                                "show": label in defaults,
                                "order": index + 1 if label in defaults else len(defaults) + index + 1,
                                "track": label,
                            }
                            for index, label in enumerate(rows_in_order)
                        ]
                    )
                    st.session_state["bw-track-rows"] = table
                checked = int(table["show"].fillna(False).sum())
                with st.expander(
                    i18n.t(
                        "res.track.selected",
                        n=i18n.number(checked),
                        total=i18n.number(len(sources)),
                    ),
                    expanded=False,
                ):
                    edited = st.data_editor(
                        table,
                        key="bw-track-editor",
                        hide_index=True,
                        width="stretch",
                        column_config={
                            "show": st.column_config.CheckboxColumn(
                                i18n.t("res.track.column.show")
                            ),
                            "order": st.column_config.NumberColumn(
                                i18n.t("res.track.column.order"), min_value=1, step=1
                            ),
                            "track": st.column_config.TextColumn(
                                i18n.t("res.track.column.track"), disabled=True
                            ),
                        },
                    )
                st.session_state["bw-track-rows"] = edited
                picked = [
                    str(label)
                    for label in edited.loc[edited["show"].fillna(False)].sort_values("order")["track"]
                    if str(label) in sources
                ]

                if not picked:
                    st.warning(i18n.t("res.track.pick_one"))
                else:
                    trackset = {name: sources[name] for name in picked}
                    # The chromosome is read end to end every time, so zooming out
                    # of the opened locus still shows the rest of the curve. One
                    # bin of that overview is 227 kb wide, which is why the
                    # window is read a second time at bin resolution and spliced
                    # in: without it the jump would land between two bins.
                    context, notes = tracks.read_tracks(
                        trackset,
                        chrom=chrom,
                        start=0,
                        end=chrom_len,
                        bins=bins,
                        smooth_bp=smooth_bp,
                        # Hundreds of kb averaged together looks the same
                        # everywhere; the max keeps the peaks in the overview.
                        aggregate="max",
                    )
                    frame = context
                    stats_frame = context
                    if not overview:
                        detail, detail_notes = tracks.read_tracks(
                            trackset,
                            chrom=chrom,
                            start=start,
                            end=end,
                            bins=bins,
                            smooth_bp=smooth_bp,
                        )
                        notes.extend(detail_notes)
                        if not detail.empty:
                            frame = tracks.merge_resolution(
                                context, detail, start=start, end=end
                            )
                            stats_frame = detail
                    for note in notes:
                        st.warning(note)
                    if frame.empty:
                        empty_state(i18n.t("res.track.no_bins"))
                    else:
                        order = [name for name in picked if name in set(frame["track"])]
                        view_start, view_end = (0, chrom_len) if overview else (start, end)
                        bands = (
                            tracks.peaks_in(
                                peaks,
                                chrom,
                                0 if overview else max(0, start - tracks.WIDE_FLANK),
                                chrom_len if overview else end + tracks.WIDE_FLANK,
                                limit=4000,
                            )
                            if show_peaks
                            else []
                        )
                        visible = [
                            band
                            for band in bands
                            if band["start"] < view_end and band["end"] > view_start
                        ]
                        # On the whole chromosome a score-sorted cut would shade
                        # one hotspot and nothing else, so the bands are spread
                        # over the span instead.
                        drawn = tracks.spread_peaks(
                            visible, tracks.peak_draw_cap(len(order))
                        )
                        gene_rows, exon_rows = tracks.locus_genes(
                            genes, exons, chrom, start, end
                        )
                        # A gene body is a fraction of a pixel across a whole
                        # chromosome, so the annotation row only appears once
                        # the locus is opened up.
                        show_model = not overview
                        if gene:
                            focus = (gene["start"], gene["end"])
                        elif overview:
                            # A region query has no gene body to point at, so the
                            # window itself is what the marker has to locate.
                            focus = (start, end)
                        else:
                            focus = None
                        st.plotly_chart(
                            tracks.build_chart(
                                frame,
                                chrom=chrom,
                                start=0,
                                end=chrom_len,
                                view=(view_start, view_end),
                                order=order,
                                peaks=drawn,
                                focus=focus,
                                genes=gene_rows if show_model else None,
                                exons=exon_rows if show_model else None,
                                unified=unified,
                                # Only View or Whole chromosome moves the picture;
                                # any other rerun keeps the zoom the reader dragged.
                                uid=f"{st.session_state.get('bw-applied', '')}"
                                f"|{st.session_state.get('bw-jump', 0)}",
                            ),
                            width="stretch",
                            config={"displaylogo": False},
                        )
                        span = view_end - view_start
                        caption = [
                            f"{chrom}:{i18n.number(view_start)}-{i18n.number(view_end)}",
                            i18n.t(
                                "res.track.span",
                                span=i18n.number(span),
                                bins=i18n.number(bins),
                                per_bin=i18n.number(span / bins, 1),
                            ),
                            (
                                i18n.t("res.track.bin_max")
                                if overview
                                else i18n.t(
                                    "res.track.context",
                                    per_bin=i18n.number(chrom_len / bins),
                                )
                            ),
                        ]
                        if smooth_bp:
                            caption.append(i18n.t("res.track.smoothed", bp=i18n.number(smooth_bp)))
                        if show_peaks:
                            caption.append(
                                i18n.t(
                                    "res.track.peak_count",
                                    seen=i18n.number(len(visible)),
                                    drawn=i18n.number(len(drawn)),
                                )
                            )
                        if overview:
                            caption.append(
                                i18n.t(
                                    "res.track.window",
                                    start=i18n.number(start),
                                    end=i18n.number(end),
                                )
                            )
                        st.caption(" · ".join(caption))
                        with st.expander(
                            i18n.t(
                                "res.track.stats_chrom"
                                if overview
                                else "res.track.stats_window"
                            ),
                            expanded=False,
                        ):
                            stats = tracks.window_stats(stats_frame, order)
                            stats.columns = [
                                i18n.t(f"res.track.stat.{name}") for name in stats.columns
                            ]
                            st.dataframe(stats, hide_index=True, width="stretch")
