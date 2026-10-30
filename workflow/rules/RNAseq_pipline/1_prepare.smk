rule rnaseq_1_prepare:
    input:
        # The config file used to be an input here, so editing any knob - even
        # one that only the joint stage reads - invalidated prepare and pulled
        # the whole clean/align/quant chain into the DAG. Everything prepare
        # actually reads is either the metadata file (an input already) or is
        # baked into the output path below (omics, dataset), so the config is
        # not needed here.
        metadata=RNA_METADATA
    output:
        manifest=f"{RNA_ROOT}/1_rawdata/rawdata_manifest.csv",
        finished_flag=f"{RNA_ROOT}/1_rawdata/finished.txt"
    params:
        omics="RNAseq",
        dataset=RNA_DATASET,
        outdir=f"{RNA_ROOT}/1_rawdata",
        script="workflow/scripts/common/prepare_rawdata.py"
    log:
        "logs/RNAseq/1_prepare.log"
    shell:
        """
        python {params.script} \\
            --omics {params.omics} \\
            --dataset {params.dataset} \\
            --outdir {params.outdir} > {log} 2>&1
        """
