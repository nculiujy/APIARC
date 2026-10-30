
QC_SUMMARY = f"{CHIP_ROOT}/3_align/qc_summary.tsv"
CHIP_SAMPLE_IDS = [s.sample for s in CHIP_SAMPLES]
CHIP_PEAK_BEDS = [f"{PEAK_DIR}/{sid}/{sid}_peaks.bed" for sid in CHIP_SAMPLE_IDS]
MERGED_PEAK_BEDS = [
    f"{DEEPT_DIR}/{group}/{group}_peaks_merged.bed" for group in CHIP_GROUPS
]

rule chipseq_4_peakcalling:
    input:
        align_flag=f"{CHIP_ROOT}/3_align/finished.txt"
    output:
        manifest=f"{PEAK_DIR}/bam_manifest.tsv",
        summary=f"{PEAK_DIR}/peak_calling_summary.tsv",
        beds=CHIP_PEAK_BEDS
    params:
        script="workflow/scripts/ChIPseq_pipline/4_1_peakcalling.py",
        metadata=CHIP_METADATA,
        bamdir=f"{CHIP_ROOT}/3_align",
        outdir=PEAK_DIR,
        config_file=primary_config_file(),
        assembly=ASSEMBLY
    threads: config.get("threads", 8)
    log:
        "logs/ChIPseq/4_1_peakcalling.log"
    conda:
        "../../envs/ChIPseq_4_peak_result.yaml"
    shell:
        """
        python {params.script} \\
            --metadata {params.metadata} \\
            --bamdir {params.bamdir} \\
            --outdir {params.outdir} \\
            --config {params.config_file} \\
            --assembly {params.assembly} \\
            --threads {threads} > {log} 2>&1
        """

rule chipseq_4_merged_peaks:
    input:
        beds=CHIP_PEAK_BEDS
    output:
        merged=MERGED_PEAK_BEDS,
        summary=f"{DEEPT_DIR}/merged_peaks_summary.tsv"
    params:
        script="workflow/scripts/ChIPseq_pipline/4_5_merge_peaks.py",
        metadata=CHIP_METADATA,
        peakdir=PEAK_DIR,
        outdir=DEEPT_DIR,
        min_support=config.get("chip", {}).get("min_peak_support", 1)
    threads: 2
    log:
        "logs/ChIPseq/4_5_merge_peaks.log"
    conda:
        "../../envs/ChIPseq_4_peak_result.yaml"
    shell:
        """
        python {params.script} \\
            --metadata {params.metadata} \\
            --peakdir {params.peakdir} \\
            --outdir {params.outdir} \\
            --min-support {params.min_support} > {log} 2>&1
        """

rule chipseq_4_annotation:
    input:
        peak=f"{PEAK_DIR}/{{sample}}/{{sample}}_peaks.bed"
    output:
        table=f"{ANNO_DIR}/{{sample}}/{{sample}}_peak_anno.csv",
        pie=f"{ANNO_DIR}/{{sample}}/{{sample}}_peak_pie.pdf"
    params:
        script="workflow/scripts/ChIPseq_pipline/4_2_annoChIPPeaks.R",
        config_file=primary_config_file(),
        assembly=ASSEMBLY,
        txdb=ANNO["txdb"],
        orgdb=ANNO["orgdb"]
    log:
        "logs/ChIPseq/4_2_annotation_{sample}.log"
    conda:
        "../../envs/ChIPseq_4_peak_result.yaml"
    shell:
        """
        Rscript {params.script} \\
            --peakfile {input.peak} \\
            --output {ANNO_DIR}/{wildcards.sample} \\
            --name {wildcards.sample} \\
            --config {params.config_file} \\
            --assembly {params.assembly} \\
            --txdb {params.txdb} \\
            --orgdb {params.orgdb} > {log} 2>&1
        """

rule chipseq_4_coverage:
    input:
        manifest=f"{PEAK_DIR}/bam_manifest.tsv",
        merged=MERGED_PEAK_BEDS
    output:
        flag=f"{DEEPT_DIR}/deeptools_finished.txt"
    params:
        script="workflow/scripts/ChIPseq_pipline/4_3_deeptools.py",
        metadata=CHIP_METADATA,
        merged_dir=DEEPT_DIR,
        outdir=DEEPT_DIR,
        norm=config.get("chip", {}).get("coverage_norm", "BPM")
    threads: config.get("threads", 8)
    log:
        "logs/ChIPseq/4_3_deeptools.log"
    conda:
        "../../envs/ChIPseq_4_peak_result.yaml"
    shell:
        """
        python {params.script} \\
            --metadata {params.metadata} \\
            --manifest {input.manifest} \\
            --merged-dir {params.merged_dir} \\
            --outdir {params.outdir} \\
            --norm {params.norm} \\
            --threads {threads} > {log} 2>&1

        touch {output.flag}
        """

rule chipseq_4_qc:
    input:
        manifest=f"{PEAK_DIR}/bam_manifest.tsv",
        summary=f"{PEAK_DIR}/peak_calling_summary.tsv",
        beds=CHIP_PEAK_BEDS
    output:
        qc=QC_SUMMARY
    params:
        script="workflow/scripts/ChIPseq_pipline/4_4_qc_summary.py",
        metadata=CHIP_METADATA,
        peakdir=PEAK_DIR,
        bamdir=f"{CHIP_ROOT}/3_align",
        outdir=f"{CHIP_ROOT}/3_align",
        frip=config.get("qc", {}).get("thresholds", {}).get("frip", 0.05),
        # Thresholds belong to the project, not to the script defaults: a project
        # that relaxes the read gate must actually get the relaxed gate.
        min_mapped=config.get("qc", {}).get("thresholds", {}).get("min_mapped_reads", 10_000_000),
        fail_on_threshold_arg="--fail-on-threshold" if config.get("qc", {}).get("fail_on_threshold", False) else ""
    threads: 2
    log:
        "logs/ChIPseq/4_4_qc_summary.log"
    conda:
        "../../envs/ChIPseq_4_peak_result.yaml"
    shell:
        """
        python {params.script} \\
            --metadata {params.metadata} \\
            --manifest {input.manifest} \\
            --peakdir {params.peakdir} \\
            --bamdir {params.bamdir} \\
            --outdir {params.outdir} \\
            --frip-threshold {params.frip} \\
            --min-mapped {params.min_mapped} \\
            {params.fail_on_threshold_arg} > {log} 2>&1
        """

rule chipseq_4_peak_result_finished:
    input:
        qc=QC_SUMMARY,
        annotation=expand(f"{ANNO_DIR}/{{sample}}/{{sample}}_peak_anno.csv",
                          sample=CHIP_SAMPLE_IDS),
        coverage=f"{DEEPT_DIR}/deeptools_finished.txt",
        merged=MERGED_PEAK_BEDS
    output:
        flag=f"{PEAK_ROOT}/finished.txt"
    run:
        with open(output.flag, "w") as handle:
            handle.write("input\t" + "\n".join([str(path) for path in input]) + "\n")
