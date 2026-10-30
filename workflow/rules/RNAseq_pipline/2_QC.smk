rule rnaseq_2_qc:
    input:
        raw_flag=f"{RNA_ROOT}/1_rawdata/finished.txt"
    output:
        qc_flag=f"{RNA_ROOT}/2_cleandata/finished.txt"
    params:
        input_dir=f"{RNA_ROOT}/1_rawdata",
        out_dir=f"{RNA_ROOT}/2_cleandata",
        script="workflow/scripts/RNAseq_pipline/2_QC_fqfile.pl",
        max_parallel=config.get("max_parallel", 4)
    threads: config.get("threads", 8)
    conda:
        "../../envs/RNAseq_2_QC.yaml"
    log:
        "logs/RNAseq/2_QC.log"
    shell:
        """
        mkdir -p {params.out_dir}
        perl {params.script} -i {params.input_dir} -o {params.out_dir} --threads {threads} --maxparallel {params.max_parallel} > {log} 2>&1
        if [ -z "$(find {params.out_dir} -name '*.clean.fastq*' -o -name '*.clean.fq*' | head -n 1)" ]; then
            echo "[ERROR] QC produced no clean reads" | tee -a {log}
            exit 1
        fi
        touch {output.qc_flag}
        """
