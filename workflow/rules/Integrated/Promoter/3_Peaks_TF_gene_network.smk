rule integrated_3_motif_enrichment:
    input:
        yaml = "config/config.yaml",
        bed_files_flag = "result/Integrated/Promoter/1_Integrated_genes/finished.txt"
    output:
        finished_flag = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/motif_enrichment.finished"
    params:
        bed_dir = "result/Integrated/Promoter/1_Integrated_genes",
        script = "workflow/scripts/Integrated/Promoter/3_1_TF_motif_enrichment_pipeline.pl",
        motifdir = "workflow/resources/Motif_tf_anno/JASPAR_CORE",
        motifmapfile = "workflow/resources/Motif_tf_anno/JASPAR_CORE_ID_to_NAME.txt",
        bw_dir = "result/ChIPseq_pipline/3_CC_pipline",
        out_dir = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}"
    log:
        "logs/Integrated/Promoter/3_Peaks_TF_gene_network/{group}_3_1_TF_motif_enrichment.log"
    shell:
        """
        mkdir -p {params.out_dir}
        perl {params.script} \
            --motifdir {params.motifdir} \
            --motifmapfile {params.motifmapfile} \
            --groupfile {input.yaml} \
            --bwfile {params.bw_dir} \
            --beddir {params.bed_dir} \
            --groupname {wildcards.group} \
            --outputdir {params.out_dir} \
            --outpic {params.out_dir} > {log} 2>&1
            
        touch {output.finished_flag}
        """

rule integrated_3_extract_motifs:
    input:
        yaml = "config/config.yaml",
        motif_finished = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/motif_enrichment.finished"
    output:
        tf_gene_list = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/Motifs_list.txt",
        finished_flag = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/extract_motifs.finished"
    params:
        script = "workflow/scripts/Integrated/Promoter/3_2_extr_enriched_motifs.pl",
        motif_indir = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}"
    log:
        "logs/Integrated/Promoter/3_Peaks_TF_gene_network/{group}_3_2_extr_enriched_motifs.log"
    shell:
        """
        perl {params.script} \
            --motifdir {params.motif_indir} \
            --yamlfile {input.yaml} \
            --outfile {output.tf_gene_list} > {log} 2>&1
            
        touch {output.finished_flag}
        """

rule integrated_3_plot_motifs:
    input:
        extracted_finished = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/extract_motifs.finished",
        rnaseq_deg = "result/RNAseq_pipline/4_DEseq/DEG_result.csv",
        tf_gene_list = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/Motifs_list.txt"
    output:
        finished_flag = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/plot_motifs.finished"
    params:
        script = "workflow/scripts/Integrated/Promoter/3_3_TF_motif_plot.R",
        out_dir = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}"
    log:
        "logs/Integrated/Promoter/3_Peaks_TF_gene_network/{group}_3_3_TF_motif_plot.log"
    shell:
        """
        Rscript {params.script} \
            --RNAseq {input.rnaseq_deg} \
            --TF_gene {input.tf_gene_list} \
            --output {params.out_dir}/{wildcards.group}_TF_gene_merged.txt \
            --picture {params.out_dir}/{wildcards.group}_TF_motif_heatmap.pdf > {log} 2>&1
            
        touch {output.finished_flag}
        """

rule integrated_3_format_motifs:
    input:
        plots_finished = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/plot_motifs.finished",
        motif_file = "workflow/resources/Motif_tf_anno/JASPAR_CORE_ID_to_NAME.txt"
    output:
        finished_flag = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/finished.txt"
    params:
        script = "workflow/scripts/Integrated/Promoter/3_4_format_TF_motifs.pl",
        out_dir = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}"
    log:
        "logs/Integrated/Promoter/3_Peaks_TF_gene_network/{group}_3_4_format_TF_motifs.log"
    shell:
        """
        perl {params.script} \
            --motiffile {input.motif_file} \
            --outdir {params.out_dir} > {log} 2>&1
            
        touch {output.finished_flag}
        """

rule integrated_3_all:
    input:
        expand("result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/finished.txt", group=INTEGRATED_GROUPS)

