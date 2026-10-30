
COMMON_ROOT = f"{JOINT_ROOT}/4_common"
COMMON_TABLES = f"{COMMON_ROOT}/tables"
COMMON_GENES_TABLE = f"{COMMON_TABLES}/common_genes.tsv"
COMMON_EXPRESSION = f"{COMMON_TABLES}/expression_matrix.tsv"
COMMON_SIGNAL = f"{COMMON_TABLES}/signal_matrix.tsv"
COMMON_SIGNAL_BY_REGION = f"{COMMON_TABLES}/signal_by_region.tsv"
COMMON_SAMPLE_MAP = f"{COMMON_TABLES}/sample_map.tsv"
COMMON_STATS = f"{COMMON_TABLES}/common_stats.tsv"
COMMON_PLOTS = f"{COMMON_ROOT}/plots"


def normalized_count_inputs():
    """normalized_counts.tsv per comparison, in the order 5_diffbind writes them.

    The directory is named after the treatment group ("H3K27ac_TGFb12h"), not
    after the mark: the comparison tuple carries the group in its first slot and
    using it as a mark silently produced a path that never exists.
    """
    paths = []
    for _namespace, reference, treatment in CHIP_COMPARISONS:
        label = f"{treatment}_vs_{reference}"
        paths.append(f"{DIFF_ROOT}/{treatment}/{label}/normalized_counts.tsv")
    return paths


def normalized_count_args():
    """--normalized-counts comparison=path, the comparison spelled as in the tables."""
    args = []
    for path in normalized_count_inputs():
        comparison = path.split("/")[-2]
        args.append(f"--normalized-counts {comparison}={path}")
    return " ".join(args)


COMPARISON_TABLES = normalized_count_inputs()
# 3_enhancer may be switched off; the distance column is then absent rather than
# the rule unbuildable.
ENHANCER_LINKS_FOR_COMMON = [ENHANCER_LINKS] if "ENHANCER_LINKS" in globals() else []
COMMON_SETTINGS = config.get("integrated", {}).get("common", {})


rule integrated_common_genes:
    input:
        gene_table=STATS_GENE_TABLE,
        regions=STATS_REGIONS,
        rna=comparison_rna_tables(),
        normalized=COMPARISON_TABLES,
        chip_metadata=CHIP_METADATA,
        rna_metadata=RNA_METADATA,
        enhancer=ENHANCER_LINKS_FOR_COMMON
    output:
        genes=COMMON_GENES_TABLE,
        expression=COMMON_EXPRESSION,
        signal=COMMON_SIGNAL,
        signal_by_region=COMMON_SIGNAL_BY_REGION,
        sample_map=COMMON_SAMPLE_MAP,
        stats=COMMON_STATS
    params:
        script="workflow/scripts/Integrated/common/common_genes.py",
        rna_args=rna_args_labelled(),
        counts_args=normalized_count_args(),
        enhancer_args="--enhancer-links " + ENHANCER_LINKS if ENHANCER_LINKS_FOR_COMMON else "",
        rna_alpha=COMMON_SETTINGS.get("rna_alpha", 0.05),
        min_base_mean=config.get("integrated", {}).get("min_expression_baseMean", 10),
        outdir=COMMON_TABLES
    threads: 2
    log:
        "logs/Integrated/4_common_genes.log"
    conda:
        "../../envs/Integrated_stats.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \\
            --gene-table {input.gene_table} \\
            --regions {input.regions} \\
            {params.rna_args} \\
            {params.counts_args} \\
            --chip-metadata {input.chip_metadata} \\
            --rna-metadata {input.rna_metadata} \\
            {params.enhancer_args} \\
            --rna-alpha {params.rna_alpha} \\
            --min-base-mean {params.min_base_mean} \\
            --outdir {params.outdir} > {log} 2>&1
        """


rule integrated_common_figures:
    input:
        genes=COMMON_GENES_TABLE,
        expression=COMMON_EXPRESSION,
        signal=COMMON_SIGNAL,
        stats=COMMON_STATS,
        sample_map=COMMON_SAMPLE_MAP,
        rna=comparison_rna_tables()
    output:
        flag=f"{COMMON_PLOTS}/figures_finished.txt"
    params:
        script="workflow/scripts/Integrated/common/common_figures.R",
        outdir=COMMON_PLOTS,
        heatmap_rows=COMMON_SETTINGS.get("heatmap_max_rows", 100),
        rna_args=rna_args_plain()
    threads: 1
    log:
        "logs/Integrated/4_common_figures.log"
    conda:
        "../../envs/Integrated_stats.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.script} \\
            --common {input.genes} \\
            --expression {input.expression} \\
            --signal {input.signal} \\
            --stats {input.stats} \\
            --sample-map {input.sample_map} \\
            {params.rna_args} \\
            --heatmap-rows {params.heatmap_rows} \\
            --outdir {params.outdir} > {log} 2>&1

        touch {output.flag}
        """


COMMON_ENRICH = {
    "common": f"{COMMON_TABLES}/go_kegg_common.tsv",
    "all_diff": f"{COMMON_TABLES}/go_kegg_all_diff.tsv",
    "compare": f"{COMMON_TABLES}/term_compare.tsv",
    "edges": f"{COMMON_TABLES}/gene_term_edges.tsv",
    "summary": f"{COMMON_TABLES}/gene_term_summary.tsv",
    "network_edges": f"{COMMON_TABLES}/network_edges.tsv",
    "network_nodes": f"{COMMON_TABLES}/network_nodes.tsv",
}


rule integrated_common_enrichment:
    input:
        genes=COMMON_GENES_TABLE,
        gene_table=STATS_GENE_TABLE,
        rna=comparison_rna_tables()
    output:
        common=COMMON_ENRICH["common"],
        all_diff=COMMON_ENRICH["all_diff"],
        compare=COMMON_ENRICH["compare"],
        edges=COMMON_ENRICH["edges"],
        summary=COMMON_ENRICH["summary"],
        network_edges=COMMON_ENRICH["network_edges"],
        network_nodes=COMMON_ENRICH["network_nodes"]
    params:
        script="workflow/scripts/Integrated/common/common_enrichment.R",
        species=SPECIES,
        alpha=COMMON_SETTINGS.get("enrich_alpha", 0.05),
        min_gene_set=config.get("integrated", {}).get("enrichment_min_gene_set", 10),
        max_terms=COMMON_SETTINGS.get("network_max_terms", 12),
        min_base_mean=config.get("integrated", {}).get("min_expression_baseMean", 10),
        outdir=COMMON_TABLES,
        rna_args=rna_args_plain()
    threads: 4
    log:
        "logs/Integrated/4_common_enrichment.log"
    conda:
        "../../envs/Integrated_enrichment.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.script} \\
            --common {input.genes} \\
            --gene-table {input.gene_table} \\
            {params.rna_args} \\
            --species {params.species} \\
            --alpha {params.alpha} \\
            --min-gene-set {params.min_gene_set} \\
            --max-terms {params.max_terms} \\
            --min-base-mean {params.min_base_mean} \\
            --outdir {params.outdir} > {log} 2>&1
        """


rule integrated_common_network:
    input:
        edges=COMMON_ENRICH["network_edges"],
        nodes=COMMON_ENRICH["network_nodes"],
        enrich=COMMON_ENRICH["common"]
    output:
        figure=f"{COMMON_PLOTS}/network.pdf",
        bar=f"{COMMON_PLOTS}/enrichment_bar.pdf",
        flag=f"{COMMON_PLOTS}/network_finished.txt"
    params:
        script="workflow/scripts/Integrated/common/network_figures.R",
        outdir=COMMON_PLOTS,
        max_terms=COMMON_SETTINGS.get("network_max_terms", 12),
        max_genes=COMMON_SETTINGS.get("network_max_genes", 25),
        alpha=COMMON_SETTINGS.get("enrich_alpha", 0.05),
        max_bar_terms=COMMON_SETTINGS.get("enrichment_bar_top_terms", 20)
    threads: 1
    log:
        "logs/Integrated/4_common_network.log"
    conda:
        "../../envs/Integrated_stats.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.script} \\
            --edges {input.edges} \\
            --nodes {input.nodes} \\
            --enrich {input.enrich} \\
            --alpha {params.alpha} \\
            --max_terms {params.max_terms} \\
            --max_genes {params.max_genes} \\
            --max-bar-terms {params.max_bar_terms} \\
            --outdir {params.outdir} > {log} 2>&1

        touch {output.flag}
        """


rule integrated_common_finished:
    input:
        genes=COMMON_GENES_TABLE,
        figures=f"{COMMON_PLOTS}/figures_finished.txt",
        network=f"{COMMON_PLOTS}/network_finished.txt",
        enrich_common=COMMON_ENRICH["common"],
        enrich_edges=COMMON_ENRICH["edges"],
        network_edges=COMMON_ENRICH["network_edges"]
    output:
        flag=f"{COMMON_ROOT}/finished.txt"
    run:
        def count_rows(path):
            with open(path) as handle:
                return max(0, sum(1 for _ in handle) - 1)

        with open(output.flag, "w") as handle:
            handle.write(f"common_gene_rows\t{count_rows(input.genes)}\n")
            handle.write(f"enrichment_rows\t{count_rows(input.enrich_common)}\n")
            handle.write(f"gene_term_edges\t{count_rows(input.enrich_edges)}\n")
            handle.write(f"network_edges\t{count_rows(input.network_edges)}\n")
            handle.write(f"tables\t{COMMON_TABLES}\n")
            handle.write(f"plots\t{COMMON_PLOTS}\n")
            handle.write("direction_filter\tboth (concordant column distinguishes)\n")
