rule chipseq_2_qc:
    input:
        raw_flag=f"{CHIP_ROOT}/1_rawdata/finished.txt"
    output:
        qc_flag=f"{CHIP_ROOT}/2_cleandata/finished.txt"
    params:
        input_dir=f"{CHIP_ROOT}/1_rawdata",
        out_dir=f"{CHIP_ROOT}/2_cleandata",
        script="workflow/scripts/ChIPseq_pipline/2_QC.pl"
    threads: config.get("threads", 8)
    conda:
        "../../envs/ChIPseq_2_QC.yaml"
    log:
        "logs/ChIPseq/2_QC.log"
    shell:
        """
        mkdir -p {params.out_dir}
        perl {params.script} -i {params.input_dir} -o {params.out_dir} --threads {threads} > {log} 2>&1
        if [ -z "$(find {params.out_dir} -name '*.fq.gz' -o -name '*.fastq.gz' | head -n 1)" ]; then
            echo "[ERROR] QC produced no clean reads" | tee -a {log}
            exit 1
        fi
        touch {output.qc_flag}
        """
