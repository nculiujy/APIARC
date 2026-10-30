
DEG_RESULT = f"{RNA_ROOT}/4_deseq/DEG_result.csv"

COUNT_MATRIX = config.get("RNAseq", {}).get(
    "count_matrix", f"{RNA_ROOT}/3_align_quant/gene_count_matrix.csv"
)


def read_rna_groups() -> dict[str, str]:
    """sample -> group, read at parse time so the shell knows the contrasts.

    Keyed by the run ID: the count matrix columns carry the same ID, and a
    sample_name is shared by the replicates of one group.
    """
    import csv

    try:
        with open(RNA_METADATA, newline="") as handle:
            rows = list(csv.DictReader(handle))
    except OSError:
        return {}
    picked = {}
    for row in rows:
        name = (row.get("sample") or "").strip() or (row.get("sample_name") or "").strip()
        if name:
            picked[name] = (row.get("group") or "").strip()
    return picked


def rna_contrasts() -> list[tuple[str, str]]:
    """(treatment, reference) pairs, mirroring the rule in 4_1_RNAseq_DEseq.R."""
    groups = sorted({group for group in read_rna_groups().values() if group})
    if not groups:
        return []
    if "control" in groups:
        reference = "control"
    elif len(groups) == 2:
        reference = groups[0]
    else:
        return []
    return [(treat, reference) for treat in groups if treat != reference]


RNA_GROUPS_ARG = ",".join(f"{name}={group}" for name, group in read_rna_groups().items())
RNA_CONTRASTS_ARG = ",".join(f"{treat}={ref}" for treat, ref in rna_contrasts())


rule rnaseq_4_deseq:
    input:
        gene_matrix=COUNT_MATRIX,
        metadata=RNA_METADATA,
        gene_bed=gene_bed(ASSEMBLY)
    output:
        deg_result=DEG_RESULT,
        deseq_flag=f"{RNA_ROOT}/4_deseq/finished.txt"
    params:
        out_dir=f"{RNA_ROOT}/4_deseq",
        table_dir=f"{RNA_ROOT}/4_deseq/tables",
        plot_dir=f"{RNA_ROOT}/4_deseq/plots",
        deseq_script="workflow/scripts/RNAseq_pipline/4_1_RNAseq_DEseq.R",
        plot_script="workflow/scripts/RNAseq_pipline/4_2_DEseq_result_plot_pipline.R",
        groups=RNA_GROUPS_ARG,
        contrasts=RNA_CONTRASTS_ARG
    threads: 2
    log:
        "logs/RNAseq/4_DEseq.log"
    conda:
        "../../envs/RNAseq_4_DEseq.yaml"
    shell:
        """
        mkdir -p {params.out_dir}

        Rscript {params.deseq_script} \\
            --input {input.gene_matrix} \\
            --metadata {input.metadata} \\
            --output {output.deg_result} > {log} 2>&1

        # One plot set per contrast: the tables and figures of two treatments
        # describe different biological questions and must not overwrite.
        for pair in $(echo "{params.contrasts}" | tr ',' ' '); do
            treat=${{pair%%=*}}
            ref=${{pair##*=}}
            Rscript {params.plot_script} \\
                --deg {params.out_dir}/DEG_result_${{treat}}_vs_${{ref}}.csv \\
                --groups "{params.groups}" \\
                --outdir_table {params.table_dir}/${{treat}}_vs_${{ref}} \\
                --outdir_plot {params.plot_dir}/${{treat}}_vs_${{ref}} \\
                --bed {input.gene_bed} >> {log} 2>&1
        done

        touch {output.deseq_flag}
        """
