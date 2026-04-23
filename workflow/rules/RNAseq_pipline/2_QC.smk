import yaml

with open("workflow/envs/RNAseq_2_QC.yaml") as f:
    env_cfg_rnaseq_2 = yaml.safe_load(f)

rule rnaseq_2_qc:
    input:
        finished_flag="result/RNAseq_pipline/1_Rawdata/finished.txt"
    output:
        out_dir=directory("result/RNAseq_pipline/2_Cleandata")
    params:
        input_dir="result/RNAseq_pipline/1_Rawdata",
        script="workflow/scripts/RNAseq_pipline/2_QC_fqfile.pl",
        max_parallel=env_cfg_rnaseq_2.get("max_parallel", 4)
    threads: env_cfg_rnaseq_2.get("threads", 8)
    conda:
        "../../envs/RNAseq_2_QC.yaml"
    log:
        "logs/RNAseq_pipline/2_QC.log"
    shell:
        """
        perl {params.script} -i {params.input_dir} -o {output.out_dir} --threads {threads} --maxparallel {params.max_parallel} > {log} 2>&1
        """
