rule integrated_1_genes:
    input:
        rnaseq_deg="result/RNAseq_pipline/4_DEseq/DEG_result.csv",
        chipseq_flag="result/ChIPseq_pipline/4_peak_result/Peak_finished.txt"
    output:
        finished_flag="result/Integrated/Promoter/1_Integrated_genes/finished.txt"
    params:
        out_dir="result/Integrated/Promoter/1_Integrated_genes",
        plots_dir="result/Integrated/Promoter/1_Integrated_genes/plots",
        chipseq_anno_dir="result/ChIPseq_pipline/4_peak_result/annotation",
        script="workflow/scripts/Integrated/Promoter/1_Promoter_plot.R",
        yaml_file="config/config.yaml",
        matrix_list="result/Integrated/Promoter/1_Integrated_genes/matrix_list.txt"
    log:
        "logs/Integrated/Promoter/1_Integrated_genes.log"
    shell:
        """
        mkdir -p {params.plots_dir}
        

        > {params.matrix_list}
        find result/ChIPseq_pipline/4_peak_result/deeptools -name "matrix_*.gz" | while read matrix_file; do

            base=$(basename $matrix_file)
            group_name=$(echo $base | sed 's/matrix_//g' | sed 's/_vs_.*\\.gz//g')
            echo -e "$group_name\t$matrix_file" >> {params.matrix_list}
        done
        

        if [ ! -s {params.matrix_list} ]; then
            echo "No computeMatrix results found. Please check ChIPseq module 4." > {log}
            touch {output.finished_flag}
            exit 0
        fi
        
        grep "finished" {input.rnaseq_deg} > /dev/null || true
        
        Rscript {params.script} \
            --RNAseq {input.rnaseq_deg} \
            --ChIPseq {params.chipseq_anno_dir} \
            --output {params.out_dir} \
            --picture {params.plots_dir} \
            --yaml {params.yaml_file} \
            --matrix {params.matrix_list} > {log} 2>&1
            
        touch {output.finished_flag}
        """