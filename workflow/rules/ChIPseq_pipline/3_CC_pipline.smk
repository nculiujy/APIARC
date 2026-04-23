import yaml
import os

with open("workflow/envs/ChIPseq_3_CC_pipline.yaml") as f:
    env_cfg_chipseq_3 = yaml.safe_load(f)


def get_bowtie2_index():
    species = config["species"]
    if species == "homo":
        return "workflow/resources/anno/ChIPseq_anno/GRCh38/grch38"
    elif species == "mm":
        return "workflow/resources/anno/ChIPseq_anno/GRCm38/grcm38"
    return ""

rule chipseq_3_cc_pipline:
    input:
        cleandata_dir="result/ChIPseq_pipline/2_Cleandata"
    output:
        out_dir=directory("result/ChIPseq_pipline/3_CC_pipline"),
        align_flag="result/ChIPseq_pipline/3_CC_pipline/Align_finished.txt",
        bw_flag="result/ChIPseq_pipline/3_CC_pipline/BW_finished.txt"
    params:
        align_script="workflow/scripts/ChIPseq_pipline/3_1_ChIPseq.pl",
        bw_script="workflow/scripts/ChIPseq_pipline/3_2_bamtobwfile.pl",
        bowtie2_index=get_bowtie2_index(),
        picard_dir="workflow/resources/picard-2.18.2"
    threads: env_cfg_chipseq_3.get("threads", 8)
    log:
        "logs/ChIPseq_pipline/3_CC_pipline.log"
    shell:
        """
        # Step 1: Align and MarkDuplicates
        perl {params.align_script} \\
            --inputdir {input.cleandata_dir} \\
            --outputdir {output.out_dir} \\
            --indexdir {params.bowtie2_index} \\
            --picarddir {params.picard_dir} \\
            --threads {threads} > {log} 2>&1
            
        touch {output.align_flag}

        # Step 2: Convert BAM to BigWig
        perl {params.bw_script} \\
            --inputdir {output.out_dir} \\
            --threads {threads} >> {log} 2>&1
            
        touch {output.bw_flag}
        """