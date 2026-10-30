
MOTIF_ROOT = f"{JOINT_ROOT}/2_motif"
MOTIF_PLOTS = f"{MOTIF_ROOT}/plots"
MOTIF_RESOURCE = "workflow/resources/Motif_tf_anno/JASPAR_CORE"
MOTIF_MAP = "workflow/resources/Motif_tf_anno/JASPAR_CORE_ID_to_NAME.txt"


def tf_list_path():
    name = "Homo_sapiens_TF.txt" if SPECIES == "homo" else "Mus_musculus_TF.txt"
    return f"workflow/resources/Motif_tf_anno/{name}"


MOTIF_HITS = f"{MOTIF_ROOT}/motif_hits.tsv"
MOTIF_SUMMARY = f"{MOTIF_ROOT}/motif_summary.tsv"
MOTIF_COMPARE = f"{MOTIF_ROOT}/motif_background_compare.tsv"
MOTIF_VERSION = f"{MOTIF_ROOT}/motif_database.meme.version"
FIMO_HITS = f"{MOTIF_ROOT}/fimo_hits.tsv"
TF_GENE_EDGES = f"{MOTIF_ROOT}/tf_gene_edges.tsv"
TF_GENE_SUMMARY = f"{MOTIF_ROOT}/tf_gene_summary.tsv"
TF_GENE_VERSION = f"{MOTIF_ROOT}/tf_gene_network.version"
TF_GENE_FIGURE = f"{MOTIF_PLOTS}/tf_gene_network.pdf"

# The per-region edges only exist for the common genes; the path is spelled here
# because 4_common.smk is included after this file.
COMMON_GENES_FOR_MOTIF = f"{JOINT_ROOT}/4_common/tables/common_genes.tsv"

MOTIF_SETTINGS = config.get("integrated", {})


def genome_fasta_input():
    return genome_fasta(ASSEMBLY)


rule integrated_motif:
    input:
        region_sets=STATS_REGION_SETS,
        rna=comparison_rna_tables(),
        genome=genome_fasta_input(),
        motif_map=MOTIF_MAP
    output:
        hits=MOTIF_HITS,
        summary=MOTIF_SUMMARY,
        compare=MOTIF_COMPARE,
        version=MOTIF_VERSION
    params:
        script="workflow/scripts/Integrated/motif/ame_enrichment.py",
        motif_dir=MOTIF_RESOURCE,
        tf_list=tf_list_path(),
        # both: a motif has to beat the static peak set and shuffled DNA, so a
        # matrix that only reflects the GC content of the set drops out.
        background=MOTIF_SETTINGS.get("motif_background", "both"),
        kmer=MOTIF_SETTINGS.get("motif_kmer", 2),
        alpha=config.get("stats", {}).get("alpha", 0.05),
        min_regions=MOTIF_SETTINGS.get("motif_min_regions", 50),
        min_base_mean=config.get("integrated", {}).get("min_expression_baseMean", 10),
        max_background=MOTIF_SETTINGS.get("motif_max_background", 20000),
        seed=config.get("stats", {}).get("seed", 20261023),
        outdir=MOTIF_ROOT,
        rna_args=rna_args_plain()
    # AME is single-threaded, but the sets are independent: one run per worker
    # turns the 34 runs from hours into minutes. Memory is the limit, not cores.
    threads: 32
    log:
        "logs/Integrated/2_motif.log"
    conda:
        "../../envs/Integrated_motif.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} \\
            --region-sets {input.region_sets} \\
            {params.rna_args} \\
            --motif-dir {params.motif_dir} \\
            --tf-list {params.tf_list} \\
            --genome {input.genome} \\
            --background {params.background} \\
            --kmer {params.kmer} \\
            --alpha {params.alpha} \\
            --min-regions {params.min_regions} \\
            --min-base-mean {params.min_base_mean} \\
            --max-background {params.max_background} \\
            --seed {params.seed} \\
            --threads {threads} \\
            --outdir {params.outdir} > {log} 2>&1
        """


if COMMON_ENABLED:
    rule integrated_motif_network:
        input:
            hits=MOTIF_HITS,
            compare=MOTIF_COMPARE,
            common=COMMON_GENES_FOR_MOTIF,
            region_sets=STATS_REGION_SETS,
            genome=genome_fasta_input()
        output:
            fimo=FIMO_HITS,
            edges=TF_GENE_EDGES,
            summary=TF_GENE_SUMMARY,
            version=TF_GENE_VERSION
        params:
            script="workflow/scripts/Integrated/motif/tf_gene_network.py",
            motif_dir=MOTIF_RESOURCE,
            alpha=config.get("stats", {}).get("alpha", 0.05),
            max_tp_pct=MOTIF_SETTINGS.get("motif_max_tp_pct", 0.9),
            max_motifs=MOTIF_SETTINGS.get("motif_network_max_motifs", 20),
            max_regions=MOTIF_SETTINGS.get("motif_fimo_max_regions", 4000),
            window=MOTIF_SETTINGS.get("fimo_window", 200),
            thresh=MOTIF_SETTINGS.get("fimo_thresh", 1e-4),
            seed=config.get("stats", {}).get("seed", 20261023),
            outdir=MOTIF_ROOT
        threads: 1
        log:
            "logs/Integrated/2_motif_network.log"
        conda:
            "../../envs/Integrated_motif.yaml"
        shell:
            """
            mkdir -p {params.outdir}
            python {params.script} \\
                --motif-hits {input.hits} \\
                --background-compare {input.compare} \\
                --common {input.common} \\
                --region-sets {input.region_sets} \\
                --motif-dir {params.motif_dir} \\
                --genome {input.genome} \\
                --alpha {params.alpha} \\
                --max-tp-pct {params.max_tp_pct} \\
                --max-motifs {params.max_motifs} \\
                --max-regions {params.max_regions} \\
                --window {params.window} \\
                --thresh {params.thresh} \\
                --seed {params.seed} \\
                --outdir {params.outdir} > {log} 2>&1
            """

    rule integrated_motif_figures:
        input:
            edges=TF_GENE_EDGES,
            summary=TF_GENE_SUMMARY
        output:
            figure=TF_GENE_FIGURE,
            flag=f"{MOTIF_PLOTS}/figures_finished.txt"
        params:
            script="workflow/scripts/Integrated/motif/tf_gene_figures.R",
            outdir=MOTIF_PLOTS,
            min_regions=MOTIF_SETTINGS.get("network_min_hit_regions", 2),
            max_tfs=MOTIF_SETTINGS.get("network_max_tfs", 12),
            max_genes=MOTIF_SETTINGS.get("network_max_genes", 24)
        threads: 1
        log:
            "logs/Integrated/2_motif_figures.log"
        # Base graphics: the motif environment has no R at all and the statistics
        # environment carries only r-base and optparse.
        conda:
            "../../envs/Integrated_stats.yaml"
        shell:
            """
            mkdir -p {params.outdir}
            Rscript {params.script} \\
                --edges {input.edges} \\
                --summary {input.summary} \\
                --min-hit-regions {params.min_regions} \\
                --max-tfs {params.max_tfs} \\
                --max-genes {params.max_genes} \\
                --outdir {params.outdir} > {log} 2>&1

            touch {output.flag}
            """

    MOTIF_NETWORK_OUTPUTS = [FIMO_HITS, TF_GENE_EDGES, TF_GENE_SUMMARY, TF_GENE_VERSION,
                             TF_GENE_FIGURE]
else:
    MOTIF_NETWORK_OUTPUTS = []


rule integrated_motif_finished:
    input:
        hits=MOTIF_HITS,
        summary=MOTIF_SUMMARY,
        compare=MOTIF_COMPARE,
        version=MOTIF_VERSION,
        network=MOTIF_NETWORK_OUTPUTS
    output:
        flag=f"{MOTIF_ROOT}/finished.txt"
    run:
        with open(input.hits) as handle:
            cells = max(0, sum(1 for _ in handle) - 1)
        # The panel and the acceptance check both need the per-region edge count;
        # reading it here keeps it from being re-derived from a 1 MB table.
        try:
            with open(TF_GENE_EDGES) as handle:
                edges = max(0, sum(1 for _ in handle) - 1)
        except OSError:
            edges = 0
        robust = "per-region (fimo)" if input.network else "absent"
        with open(output.flag, "w") as handle:
            handle.write(f"motif_database\t{MOTIF_RESOURCE}\n")
            handle.write(f"background\t{config.get('integrated', {}).get('motif_background', 'both')}\n")
            handle.write(f"matrix_set_rows\t{cells}\n")
            handle.write(f"tf_gene_edges\t{edges}\n")
            handle.write(f"evidence\t{robust}\n")
            handle.write(f"database_version\t{input.version}\n")
