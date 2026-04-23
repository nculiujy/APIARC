import yaml

with open("workflow/envs/ChIPseq_2_QC.yaml") as f:
    env_cfg_chipseq_2 = yaml.safe_load(f)

rule chipseq_2_qc:
    input:
        input_dir="result/ChIPseq_pipline/1_Rawdata"
    output:
        out_dir=directory("result/ChIPseq_pipline/2_Cleandata")
    params:
        script="workflow/scripts/ChIPseq_pipline/2_QC.pl"
    threads: env_cfg_chipseq_2.get("threads", 8)
    conda:
        "../../envs/ChIPseq_2_QC.yaml"
    log:
        "logs/ChIPseq_pipline/2_QC.log"
    shell:
        """
        perl {params.script} -i {input.input_dir} -o {output.out_dir} --threads {threads} > {log} 2>&1
        """
