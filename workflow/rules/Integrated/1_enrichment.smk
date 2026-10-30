
ENRICH_ROOT = f"{JOINT_ROOT}/1_enrichment"
ENRICH_SUMMARY = f"{ENRICH_ROOT}/enrichment_summary.tsv"
ENRICH_UNIVERSE = f"{ENRICH_ROOT}/universe.txt"
ENRICH_CONTRAST = f"{ENRICH_ROOT}/differential_vs_static.tsv"


rule integrated_enrichment:
    input:
        region_sets=STATS_REGION_SETS,
        rna=comparison_rna_tables()
    output:
        summary=ENRICH_SUMMARY,
        universe=ENRICH_UNIVERSE
    params:
        script="workflow/scripts/Integrated/enrichment/go_kegg.R",
        universe_name=config.get("integrated", {}).get("universe", "expressed"),
        species=SPECIES,
        alpha=config.get("stats", {}).get("alpha", 0.05),
        min_gene_set=config.get("integrated", {}).get("enrichment_min_gene_set", 10),
        outdir=ENRICH_ROOT,
        rna_args=rna_args_plain()
    threads: 4
    log:
        "logs/Integrated/1_enrichment.log"
    conda:
        "../../envs/Integrated_enrichment.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.script} \\
            --region-sets {input.region_sets} \\
            {params.rna_args} \\
            --universe {params.universe_name} \\
            --species {params.species} \\
            --alpha {params.alpha} \\
            --min-gene-set {params.min_gene_set} \\
            --outdir {params.outdir} > {log} 2>&1
        """


rule integrated_enrichment_finished:
    input:
        summary=ENRICH_SUMMARY,
        universe=ENRICH_UNIVERSE
    output:
        flag=f"{ENRICH_ROOT}/finished.txt"
    params:
        universe_name=config.get("integrated", {}).get("universe", "expressed")
    run:
        import os

        contrast = "present" if os.path.exists(ENRICH_CONTRAST) else "absent"
        with open(output.flag, "w") as handle:
            handle.write(f"universe\t{params.universe_name}\n")
            handle.write(f"differential_vs_static\t{contrast}\n")
