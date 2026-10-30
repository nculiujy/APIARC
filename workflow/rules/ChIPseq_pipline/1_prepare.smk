rule chipseq_1_prepare:
    input:
        # Same reason as rnaseq_1_prepare: the config was listed as an input,
        # so a joint-stage setting made this rule and the whole peak chain
        # stale. The metadata file is the only content prepare reads.
        metadata=CHIP_METADATA
    output:
        manifest=f"{CHIP_ROOT}/1_rawdata/rawdata_manifest.csv",
        finished_flag=f"{CHIP_ROOT}/1_rawdata/finished.txt"
    params:
        omics="ChIPseq",
        dataset=CHIP_DATASET,
        outdir=f"{CHIP_ROOT}/1_rawdata",
        script="workflow/scripts/common/prepare_rawdata.py"
    log:
        "logs/ChIPseq/1_prepare.log"
    shell:
        """
        python {params.script} \\
            --omics {params.omics} \\
            --dataset {params.dataset} \\
            --outdir {params.outdir} > {log} 2>&1
        """
