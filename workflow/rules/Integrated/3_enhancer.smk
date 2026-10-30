
ENHANCER_ROOT = f"{JOINT_ROOT}/3_enhancer"
ENHANCER_LINKS = f"{ENHANCER_ROOT}/enhancer_gene_links.tsv"
ENHANCER_BINS = f"{ENHANCER_ROOT}/distance_bins.tsv"
ENHANCER_SUMMARY = f"{ENHANCER_ROOT}/enhancer_summary.tsv"


def enhancer_bins():
    bins = config.get("integrated", {}).get("enhancer_distance_bins", [0, 10000, 50000, 100000, 500000, 1000000])
    return ",".join(str(int(value)) for value in bins)


rule integrated_enhancer:
    input:
        regions=STATS_REGIONS,
        genes=gene_model_bed(),
        rna=comparison_rna_tables()
    output:
        links=ENHANCER_LINKS,
        bins=ENHANCER_BINS,
        summary=ENHANCER_SUMMARY
    params:
        script="workflow/scripts/Integrated/enhancer/link_enhancers.py",
        # The script pairs each ChIP comparison with the RNA contrast of the
        # same time point, so it needs the comparison names, not just paths.
        rna_args=rna_args_labelled(),
        bins=enhancer_bins(),
        promoter_window=config.get("integrated", {}).get("promoter_window", 2000),
        max_distance=config.get("integrated", {}).get("enhancer_max_distance", 1000000),
        random_per_region=config.get("integrated", {}).get("enhancer_random_per_region", 5),
        min_base_mean=config.get("integrated", {}).get("min_expression_baseMean", 10),
        seed=config.get("stats", {}).get("seed", 20261023),
        outdir=ENHANCER_ROOT
    threads: 2
    log:
        "logs/Integrated/3_enhancer.log"
    conda:
        "../../envs/Integrated_enhancer.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \\
            --regions {input.regions} \\
            --genes {input.genes} \\
            {params.rna_args} \\
            --bins {params.bins} \\
            --promoter-window {params.promoter_window} \\
            --max-distance {params.max_distance} \\
            --random-per-region {params.random_per_region} \\
            --min-base-mean {params.min_base_mean} \\
            --seed {params.seed} \\
            --outdir {params.outdir} > {log} 2>&1
        """


rule integrated_enhancer_finished:
    input:
        links=ENHANCER_LINKS,
        bins=ENHANCER_BINS,
        summary=ENHANCER_SUMMARY,
        genes=gene_model_bed()
    output:
        flag=f"{ENHANCER_ROOT}/finished.txt"
    params:
        distance_bins=enhancer_bins()
    run:
        with open(output.flag, "w") as handle:
            handle.write(f"gene_model\t{input.genes}\n")
            handle.write(f"distance_bins\t{params.distance_bins}\n")
            handle.write("null\tdistance-matched random positions\n")
