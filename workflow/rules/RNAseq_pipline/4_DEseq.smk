import yaml
import os
import json


import pandas as pd
def get_deseq_config_json():

    meta = pd.read_csv("config/RNAseq_metadata.csv")
    deseq_config = dict(zip(meta['sample_name'], meta['group']))
    return json.dumps(deseq_config).replace('"', '\\"')


def get_bed_file():
    species = config["species"]
    if species == "homo":

        return "workflow/resources/anno/RNAseq_anno/GRCh38/gencode.v44.annotation.bed"
    elif species == "mm":
        return "workflow/resources/anno/RNAseq_anno/GRCm38/gencode.vM25.annotation.bed"
    return ""

rule rnaseq_4_deseq:
    input:
        gene_matrix="result/RNAseq_pipline/3_RC_pipline/gene_count_matrix.csv",
        align_log="logs/RNAseq_pipline/3_RC_pipline.log"
    output:
        out_dir=directory("result/RNAseq_pipline/4_DEseq"),
        deg_result="result/RNAseq_pipline/4_DEseq/DEG_result.csv",
        table_dir=directory("result/RNAseq_pipline/4_DEseq/tables"),
        plot_dir=directory("result/RNAseq_pipline/4_DEseq/plots"),
        deseq_flag="result/RNAseq_pipline/4_DEseq/DEseq_finished.txt"
    params:
        deseq_script="workflow/scripts/RNAseq_pipline/4_1_RNAseq_DEseq.R",
        plot_script="workflow/scripts/RNAseq_pipline/4_2_DEseq_result_plot_pipline.R",
        config_json=get_deseq_config_json(),
        gtf_file=get_gtf_file(), # using the same get_gtf_file from 3_RC_pipline.smk context
        bed_file=get_bed_file()
    log:
        "logs/RNAseq_pipline/4_DEseq.log"
    shell:
        """
        mkdir -p {output.out_dir}
        

        if [ ! -f {params.bed_file} ]; then

            awk '$3 == "gene" {{
                gene_id=""; gene_name="";
                for(i=9; i<=NF; i++) {{
                    if($i=="gene_id") gene_id=$(i+1);
                    if($i=="gene_name") gene_name=$(i+1);
                }}
                gsub(/"|;/, "", gene_id);
                gsub(/"|;/, "", gene_name);
                print $1"\\t"$4"\\t"$5"\\t"gene_id"\\t0\\t"$7"\\t"gene_name
            }}' {params.gtf_file} > {params.bed_file}
        fi
        
        # Step 1: Run DESeq2
        Rscript {params.deseq_script} \\
            --input {input.gene_matrix} \\
            --config "{params.config_json}" \\
            --output {output.deg_result} > {log} 2>&1
            
        # Step 2: Plot results
        Rscript {params.plot_script} \\
            --log {input.align_log} \\
            --deg {output.deg_result} \\
            --outdir_table {output.table_dir} \\
            --outdir_plot {output.plot_dir} \\
            --bed {params.bed_file} >> {log} 2>&1
            
        touch {output.deseq_flag}
        """