
STATS_REGIONS = f"{STATS_ROOT}/regions.tsv"
STATS_REGION_SETS = f"{STATS_ROOT}/region_sets.tsv"
STATS_BEDDIR = f"{STATS_ROOT}/bed"
STATS_GENE_TABLE = f"{STATS_ROOT}/chip_gene_table.tsv"
STATS_SUMMARY = f"{STATS_ROOT}/stats_summary.tsv"
STATS_DOSE = f"{STATS_ROOT}/dose_response.tsv"
STATS_PLOTS = f"{STATS_ROOT}/plots"
RNA_DEG_TABLE = RNA_DEG_FALLBACK


def gene_model_bed():
    return gene_bed(ASSEMBLY)


def differential_region_tables():
    return [
        f"{DIFF_ROOT}/{mark}/{treatment}_vs_{reference}/differential_regions.tsv"
        for mark, reference, treatment in CHIP_COMPARISONS
    ]


def merged_peak_beds():
    """Merged consensus peaks per mark/condition group (descriptive mode)."""
    return [f"{DEEPT_DIR}/{group}/{group}_peaks_merged.bed" for group in CHIP_GROUPS]


def negative_control_tables():
    return [
        f"{DIFF_ROOT}/{mark}/{treatment}_vs_{reference}/negative_control.tsv"
        for mark, reference, treatment in CHIP_COMPARISONS
    ]


rule integrated_regions:
    input:
        genes=gene_model_bed(),
        differential=differential_region_tables(),
        merged=merged_peak_beds()
    output:
        regions=STATS_REGIONS,
        region_sets=STATS_REGION_SETS,
        beds=directory(STATS_BEDDIR)
    params:
        script="workflow/scripts/Integrated/stats/0_1_annotate_regions.py",
        mode=CHIP_MODE,
        differential_args=lambda wildcards, input: " ".join(
            f"--differential {path}" for path in input.differential
        ),
        bound_args=lambda wildcards, input: " ".join(
            f"--bound-bed {group}={path}"
            for group, path in zip(CHIP_GROUPS, input.merged)
        ),
        promoter_window=config.get("integrated", {}).get("promoter_window", 2000),
        max_distance=config.get("integrated", {}).get("enhancer_max_distance", 1000000)
    threads: 2
    log:
        "logs/Integrated/0_1_annotate_regions.log"
    conda:
        "../../envs/Integrated_stats.yaml"
    shell:
        """
        mkdir -p {STATS_ROOT}
        python {params.script} \\
            --mode {params.mode} \\
            --genes {input.genes} \\
            --promoter-window {params.promoter_window} \\
            --max-distance {params.max_distance} \\
            {params.differential_args} \\
            {params.bound_args} \\
            --out {output.regions} \\
            --beddir {output.beds} > {log} 2>&1
        """


rule integrated_gene_table:
    input:
        regions=STATS_REGIONS,
        rna=comparison_rna_tables()
    output:
        table=STATS_GENE_TABLE
    params:
        script="workflow/scripts/Integrated/stats/0_2_gene_table.py",
        min_base_mean=config.get("integrated", {}).get("min_expression_baseMean", 10),
        rna_args=rna_args_labelled()
    threads: 2
    log:
        "logs/Integrated/0_2_gene_table.log"
    conda:
        "../../envs/Integrated_stats.yaml"
    shell:
        """
        python {params.script} \\
            --regions {input.regions} \\
            {params.rna_args} \\
            --min-base-mean {params.min_base_mean} \\
            --out {output.table} > {log} 2>&1
        """


rule integrated_tests:
    input:
        gene_table=STATS_GENE_TABLE,
        rna=comparison_rna_tables(),
        negative=negative_control_tables()
    output:
        summary=STATS_SUMMARY,
        dose=STATS_DOSE,
        families=f"{STATS_ROOT}/analysis_families.tsv"
    params:
        script="workflow/scripts/Integrated/stats/0_3_tests.R",
        mode=CHIP_MODE,
        permutations=config.get("stats", {}).get("permutations", 10000),
        bootstrap=config.get("stats", {}).get("bootstrap", 2000),
        seed=config.get("stats", {}).get("seed", 20261023),
        alpha=config.get("stats", {}).get("alpha", 0.05),
        min_base_mean=config.get("integrated", {}).get("min_expression_baseMean", 10),
        outdir=STATS_ROOT,
        rna_args=rna_args_labelled(),
        negative_args=lambda wildcards, input: " ".join(
            f"--neg-control {path}" for path in input.negative
        )
    threads: 2
    log:
        "logs/Integrated/0_3_tests.log"
    conda:
        "../../envs/Integrated_stats.yaml"
    shell:
        """
        Rscript {params.script} \\
            --gene-table {input.gene_table} \\
            {params.rna_args} \\
            --mode {params.mode} \\
            {params.negative_args} \\
            --permutations {params.permutations} \\
            --bootstrap {params.bootstrap} \\
            --seed {params.seed} \\
            --alpha {params.alpha} \\
            --min-base-mean {params.min_base_mean} \\
            --outdir {params.outdir} > {log} 2>&1
        """


rule integrated_figures:
    input:
        stats=STATS_SUMMARY,
        dose=STATS_DOSE,
        gene_table=STATS_GENE_TABLE
    output:
        flag=f"{STATS_PLOTS}/figures_finished.txt"
    params:
        script="workflow/scripts/Integrated/stats/0_4_plots.R",
        outdir=STATS_PLOTS
    threads: 2
    log:
        "logs/Integrated/0_4_plots.log"
    conda:
        "../../envs/Integrated_stats.yaml"
    shell:
        """
        Rscript {params.script} \\
            --stats {input.stats} \\
            --dose {input.dose} \\
            --gene-table {input.gene_table} \\
            --outdir {params.outdir} > {log} 2>&1

        touch {output.flag}
        """


rule integrated_stats_finished:
    input:
        summary=STATS_SUMMARY,
        figures=f"{STATS_PLOTS}/figures_finished.txt"
    output:
        flag=f"{STATS_ROOT}/finished.txt",
        mode=f"{STATS_ROOT}/mode.txt"
    run:
        with open(output.mode, "w") as handle:
            handle.write(
                f"mode: {CHIP_MODE}\n"
                f"comparisons: {len(CHIP_COMPARISONS)}\n"
                "differential outputs are only produced in differential mode\n"
            )
        with open(output.flag, "w") as handle:
            handle.write(f"mode\t{CHIP_MODE}\n")
