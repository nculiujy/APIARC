
def get_bowtie2_index():
    return chip_index(ASSEMBLY)


rule chipseq_3_cc_pipline:
    input:
        qc_flag=f"{CHIP_ROOT}/2_cleandata/finished.txt"
    output:
        align_flag=f"{CHIP_ROOT}/3_align/finished.txt"
    params:
        align_script="workflow/scripts/ChIPseq_pipline/3_1_ChIPseq.pl",
        bowtie2_index=get_bowtie2_index(),
        picard_dir="workflow/resources/picard-2.18.2",
        input_dir=f"{CHIP_ROOT}/2_cleandata",
        out_dir=f"{CHIP_ROOT}/3_align"
    threads: config.get("threads", 8)
    log:
        "logs/ChIPseq/3_CC_pipline.log"
    conda:
        "../../envs/ChIPseq_3_CC_pipline.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        perl {params.align_script} \\
            --inputdir {params.input_dir} \\
            --outputdir {params.out_dir} \\
            --indexdir {params.bowtie2_index} \\
            --picarddir {params.picard_dir} \\
            --threads {threads} > {log} 2>&1

        touch {output.align_flag}
        """
