
VALIDATION_REPORT = f"{PROVENANCE_ROOT}/validation_report.tsv"

COUNT_MATRIX = config.get("RNAseq", {}).get(
    "count_matrix", f"{RNA_ROOT}/3_align_quant/gene_count_matrix.csv"
)


rule validate_design:
    input:
        chip_metadata=CHIP_METADATA,
        rna_metadata=RNA_METADATA
    output:
        report=VALIDATION_REPORT
    params:
        script="workflow/scripts/common/validate_config.py",
        mode=config.get("integrated", {}).get("mode", "auto"),
        count_matrix=COUNT_MATRIX,
        root=PROVENANCE_ROOT
    threads: 1
    log:
        "logs/validate_design.log"
    shell:
        """
        mkdir -p {params.root}
        python {params.script} \\
            --chip-metadata {input.chip_metadata} \\
            --rna-metadata {input.rna_metadata} \\
            --count-matrix {params.count_matrix} \\
            --mode {params.mode} \\
            --out {output.report} > {log} 2>&1
        """
