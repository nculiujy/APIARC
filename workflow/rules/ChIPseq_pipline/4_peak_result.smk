import yaml
import os

with open("workflow/envs/ChIPseq_4_peak_result.yaml", "w") as f:
    yaml.dump({"threads": 8, "max_parallel": 4}, f)

rule chipseq_4_peak_result:
    input:
        bw_flag="result/ChIPseq_pipline/3_CC_pipline/BW_finished.txt",
        bam_dir="result/ChIPseq_pipline/3_CC_pipline"
    output:
        out_dir=directory("result/ChIPseq_pipline/4_peak_result"),
        peak_flag="result/ChIPseq_pipline/4_peak_result/Peak_finished.txt"
    params:
        peakcalling_script="workflow/scripts/ChIPseq_pipline/4_1_peakcalling.pl",
        anno_script="workflow/scripts/ChIPseq_pipline/4_2_annoChIPPeaks.R",
        deeptools_script="workflow/scripts/ChIPseq_pipline/4_3_deeptools.py",
        metadata="config/ChIPseq_metadata.csv",
        genome=config.get("species", "mm"),
        config_file="config/config.yaml"
    threads: 8
    log:
        "logs/ChIPseq_pipline/4_peak_result.log"
    shell:
        """
        mkdir -p {output.out_dir}
        
        # Step 1: MACS2 Peak Calling
        perl {params.peakcalling_script} \\
            --metadata {params.metadata} \\
            --bamdir {input.bam_dir} \\
            --outdir {output.out_dir}/peakcalling \\
            --picdir {output.out_dir}/peakcalling/plots \\
            --genome {params.genome} \\
            --threads {threads} > {log} 2>&1
            
        # Step 2: Peak Annotation

        find {output.out_dir}/peakcalling -name "*_peaks_tab.bed" | while read bed_file; do
            group_name=$(basename $(dirname $bed_file))
            anno_outdir="{output.out_dir}/annotation/$group_name"
            
            Rscript {params.anno_script} \\
                --peakfile $bed_file \\
                --output $anno_outdir \\
                --name $group_name \\
                --config {params.config_file} >> {log} 2>&1
        done
        
        # Step 3: Deeptools (Heatmap & Profile)
        python {params.deeptools_script} \\
            --metadata {params.metadata} \\
            --bamdir {input.bam_dir} \\
            --peakdir {output.out_dir}/peakcalling \\
            --outdir {output.out_dir}/deeptools \\
            --threads {threads} >> {log} 2>&1
            
        touch {output.peak_flag}
        """