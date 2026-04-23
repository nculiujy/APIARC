import yaml

with open("workflow/envs/RNAseq_1_download.yaml") as f:
    env_cfg_rnaseq_1 = yaml.safe_load(f)

rule rnaseq_1_download:
    input:


        "config/config.yaml"
    output:
        result_dir=directory("result/RNAseq_pipline/1_Rawdata"),
        finished_flag="result/RNAseq_pipline/1_Rawdata/finished.txt"
    params:
        species=config["species"],
        experiment=config.get("experiment", "GSEXXXXXX"),
        script="workflow/scripts/RNAseq_pipline/1_download.py",

        rawdata_dir=f"config/{config.get('experiment', 'GSEXXXXXX')}/RNAseq"
    threads: env_cfg_rnaseq_1.get("threads", 8)
    conda:
        "../../envs/RNAseq_1_download.yaml"
    log:
        "logs/RNAseq_pipline/1_download.log"
    shell:
        """
        export DOWNLOAD_THREADS={threads}
        python {params.script} {params.rawdata_dir} {output.result_dir} > {log} 2>&1
        touch {output.finished_flag}
        """
