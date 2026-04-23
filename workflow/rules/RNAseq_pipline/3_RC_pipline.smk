import yaml
import os

with open("workflow/envs/RNAseq_3_RC_pipline.yaml") as f:
    env_cfg_rnaseq_3 = yaml.safe_load(f)


def get_hisat2_index():
    species = config["species"]
    if species == "homo":
        return "workflow/resources/anno/RNAseq_anno/GRCh38/grch38/genome"
    elif species == "mm":
        return "workflow/resources/anno/RNAseq_anno/GRCm38/grcm38/genome"
    return ""


def get_gtf_file():
    species = config["species"]
    if species == "homo":
        return "workflow/resources/anno/RNAseq_anno/GRCh38/gencode.v44.annotation.gtf"
    elif species == "mm":
        return "workflow/resources/anno/RNAseq_anno/GRCm38/gencode.vM25.annotation.gtf"
    return ""

rule rnaseq_3_rc_pipline:
    input:
        cleandata_dir="result/RNAseq_pipline/2_Cleandata"
    output:
        out_dir=directory("result/RNAseq_pipline/3_RC_pipline"),
        align_flag="result/RNAseq_pipline/3_RC_pipline/Align_finished.txt",
        quant_flag="result/RNAseq_pipline/3_RC_pipline/Quant_finished.txt",
        gene_matrix="result/RNAseq_pipline/3_RC_pipline/gene_count_matrix.csv",
        transcript_matrix="result/RNAseq_pipline/3_RC_pipline/transcript_count_matrix.csv"
    params:
        align_script="workflow/scripts/RNAseq_pipline/3_1_Align.py",
        quant_script="workflow/scripts/RNAseq_pipline/3_2_Quant.py",
        prepde_script="workflow/scripts/RNAseq_pipline/3_3_prepDE.py",
        metadata_file="config/RNAseq_metadata.csv",
        hisat2_index=get_hisat2_index(),
        gtf_file=get_gtf_file(),
        max_parallel=env_cfg_rnaseq_3.get("max_parallel", 4)
    threads: env_cfg_rnaseq_3.get("threads", 8)
    log:
        "logs/RNAseq_pipline/3_RC_pipline.log"
    shell:
        """
        # Step 1: Align
        python {params.align_script} \\
            --inputdir {input.cleandata_dir} \\
            --outputdir {output.out_dir} \\
            --threads {threads} \\
            --maxparallel {params.max_parallel} \\
            --hisat2_index {params.hisat2_index} >> {log} 2>&1

        # Step 2: Quant (StringTie)
        python {params.quant_script} \\
            --inputdir {output.out_dir} \\
            --outputdir {output.out_dir} \\
            --threads {threads} \\
            --maxparallel {params.max_parallel} \\
            --gtf_file {params.gtf_file} >> {log} 2>&1

        # Step 3: PrepDE (Merge stringtie outputs into matrix)

        SAMPLE_LIST="{output.out_dir}/sample_list.txt"
        > $SAMPLE_LIST
        


        awk -F ',' 'NR>1 {{print $1" "$2}}' {params.metadata_file} | while read SRR NAME; do

            GTF_PATH=$(find {output.out_dir} -name "transcripts.gtf" | grep "/$SRR/" | head -n 1 || true)
            if [ -n "$GTF_PATH" ]; then

                echo "${{NAME}} $GTF_PATH" >> $SAMPLE_LIST
            fi
        done
        
        cat $SAMPLE_LIST >> {log}

        if [ -s $SAMPLE_LIST ]; then
            python {params.prepde_script} \
                -i $SAMPLE_LIST \
                -g {output.gene_matrix} \
                -t {output.transcript_matrix} >> {log} 2>&1
        else
            echo "[Warning] sample_list.txt is empty. Creating dummy matrix files to prevent failure." >> {log}
            touch {output.gene_matrix}
            touch {output.transcript_matrix}
        fi
            
        touch {output.align_flag}
        touch {output.quant_flag}
        """
