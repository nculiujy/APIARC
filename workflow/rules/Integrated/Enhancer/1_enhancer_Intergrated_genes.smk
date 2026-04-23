rule enhancer_integrated_1_distance:
    input:
        yaml = "config/config.yaml",
        promoter_flag = "result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/finished.txt"
    output:
        enh_gene = "result/Integrated/Enhancer/1_Integrated_genes/{group}/enh_gene_within_1mb.bed",
        rand_enh = "result/Integrated/Enhancer/1_Integrated_genes/{group}/random_enh_gene_within_1mb.bed"
    params:
        genefile = "result/Integrated/Promoter/1_Integrated_genes/{group}_all_peaks.bed",
        script = "workflow/scripts/Integrated/Enhancer/1_1_calculate_eRNA_gene_distance.pl",
        enhfile = "workflow/resources/Enhancer_anno/refGene.ncbiRefSeq.genecode.filtered.bed"
    log:
        "logs/Integrated/Enhancer/1_Integrated_genes/{group}_distance.log"
    shell:
        """
        mkdir -p $(dirname {output.enh_gene})
        mkdir -p $(dirname {log})
        
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting calculate_eRNA_gene_distance.pl for {wildcards.group}..." > {log}
        
        perl {params.script} \
            --enhfile {params.enhfile} \
            --genefile {params.genefile} \
            --outfile {output.enh_gene} \
            --randbed {output.rand_enh} >> {log} 2>&1
            
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished calculate_eRNA_gene_distance.pl successfully." >> {log}
        """

rule enhancer_integrated_1_signal:
    input:
        enh = "result/Integrated/Enhancer/1_Integrated_genes/{group}/enh_gene_within_1mb.bed",
        peak_flag = "result/ChIPseq_pipline/4_peak_result/Peak_finished.txt"
    output:
        matrix_output = "result/Integrated/Enhancer/1_Integrated_genes/{group}/enh_gene_within_1mb_matrix.gz",
        tab_output = "result/Integrated/Enhancer/1_Integrated_genes/{group}/enh_gene_within_1mb_signal.tab"
    params:
        bw_dir = "result/ChIPseq_pipline/4_peak_result/deeptools/{group}_vs_Input"
    log:
        "logs/Integrated/Enhancer/1_Integrated_genes/{group}_signal.log"
    threads: 10
    shell:
        """
        export PATH="/home/jyliu/miniconda3/envs/ChIPseq_Pipline/bin:$PATH"
        export MPLCONFIGDIR=/tmp/matplotlib_cache
        mkdir -p $MPLCONFIGDIR
        mkdir -p $(dirname {log})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting computeMatrix for TRUE Enhancer-Gene signal ({wildcards.group})..." > {log}
        bwfiles=$(ls {params.bw_dir}/*.bw | tr '\n' ' ')
        echo "Found bw files: $bwfiles" >> {log}
        
        computeMatrix reference-point -R {input.enh} -S $bwfiles \
          -b 1000 -a 1000 --binSize 50 -p {threads} -o {output.matrix_output} --outFileNameMatrix {output.tab_output} >> {log} 2>&1
          
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished computeMatrix successfully." >> {log}
        """

rule enhancer_integrated_1_rand_signal:
    input:
        randbed = "result/Integrated/Enhancer/1_Integrated_genes/{group}/random_enh_gene_within_1mb.bed",
        peak_flag = "result/ChIPseq_pipline/4_peak_result/Peak_finished.txt"
    output:
        matrix_output = "result/Integrated/Enhancer/1_Integrated_genes/{group}/random_enh_gene_within_1mb_matrix.gz",
        tab_output = "result/Integrated/Enhancer/1_Integrated_genes/{group}/random_enh_gene_within_1mb_signal.tab"
    params:
        bw_dir = "result/ChIPseq_pipline/4_peak_result/deeptools/{group}_vs_Input"
    log:
        "logs/Integrated/Enhancer/1_Integrated_genes/{group}_rand_signal.log"
    threads: 10
    shell:
        """
        export PATH="/home/jyliu/miniconda3/envs/ChIPseq_Pipline/bin:$PATH"
        export MPLCONFIGDIR=/tmp/matplotlib_cache
        mkdir -p $MPLCONFIGDIR
        mkdir -p $(dirname {log})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting computeMatrix for RANDOM Enhancer-Gene signal ({wildcards.group})..." > {log}
        bwfiles=$(ls {params.bw_dir}/*.bw | tr '\n' ' ')
        echo "Found bw files: $bwfiles" >> {log}
        
        computeMatrix reference-point -R {input.randbed} -S $bwfiles \
          -b 1000 -a 1000 --binSize 50 -p {threads} -o {output.matrix_output} --outFileNameMatrix {output.tab_output} >> {log} 2>&1
          
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished computeMatrix for random signal successfully." >> {log}
        """

rule enhancer_integrated_1_true_enhancer:
    input:
        matrix1 = "result/Integrated/Enhancer/1_Integrated_genes/{group}/enh_gene_within_1mb_signal.tab",
        matrix2 = "result/Integrated/Enhancer/1_Integrated_genes/{group}/random_enh_gene_within_1mb_signal.tab",
        enh = "result/Integrated/Enhancer/1_Integrated_genes/{group}/enh_gene_within_1mb.bed"
    output:
        true_enh = "result/Integrated/Enhancer/1_Integrated_genes/{group}/Ture_enh_gene_within_1mb.bed"
    params:
        script = "workflow/scripts/Integrated/Enhancer/1_2_Calculate_signal.R"
    log:
        "logs/Integrated/Enhancer/1_Integrated_genes/{group}_true_enhancer.log"
    shell:
        """
        # The original R script expects comma separated sample names to process
        samples=$(ls result/ChIPseq_pipline/4_peak_result/deeptools/{wildcards.group}_vs_Input/*.bw | xargs -n 1 basename | sed 's/.bw//' | tr '\n' ',' | sed 's/,$//')
        
        Rscript {params.script} \
            -m {input.matrix1} \
            -n {input.matrix2} \
            -s $samples \
            -a {input.enh} \
            -o {output.true_enh} > {log} 2>&1
        """

rule enhancer_integrated_1_true_signal:
    input:
        true_enh = "result/Integrated/Enhancer/1_Integrated_genes/{group}/Ture_enh_gene_within_1mb.bed",
        peak_flag = "result/ChIPseq_pipline/4_peak_result/Peak_finished.txt"
    output:
        matrix_output = "result/Integrated/Enhancer/1_Integrated_genes/{group}/Ture_enh_gene_within_1mb_matrix.gz",
        tab_output = "result/Integrated/Enhancer/1_Integrated_genes/{group}/Ture_enh_gene_within_1mb_signal.tab"
    params:
        bw_dir = "result/ChIPseq_pipline/4_peak_result/deeptools/{group}_vs_Input"
    log:
        "logs/Integrated/Enhancer/1_Integrated_genes/{group}_true_signal.log"
    threads: 10
    shell:
        """
        export PATH="/home/jyliu/miniconda3/envs/ChIPseq_Pipline/bin:$PATH"
        export MPLCONFIGDIR=/tmp/matplotlib_cache
        mkdir -p $MPLCONFIGDIR
        bwfiles=$(ls {params.bw_dir}/*.bw | tr '\n' ' ')
        computeMatrix reference-point -R {input.true_enh} -S $bwfiles \
          -b 1000 -a 1000 --binSize 50 -p {threads} -o {output.matrix_output} --outFileNameMatrix {output.tab_output} >> {log} 2>&1
        """

rule enhancer_integrated_1_plot_profile:
    input:
        matrix = "result/Integrated/Enhancer/1_Integrated_genes/{group}/Ture_enh_gene_within_1mb_matrix.gz"
    output:
        plot = "result/Integrated/Enhancer/1_Integrated_genes/{group}/plots/Ture_enh_gene_within_1mb_profile.pdf"
    log:
        "logs/Integrated/Enhancer/1_Integrated_genes/{group}_plot_profile.log"
    shell:
        """
        export PATH="/home/jyliu/miniconda3/envs/ChIPseq_Pipline/bin:$PATH"
        export MPLCONFIGDIR=/tmp/matplotlib_cache
        mkdir -p $MPLCONFIGDIR
        samples=$(ls result/ChIPseq_pipline/4_peak_result/deeptools/{wildcards.group}_vs_Input/*.bw | xargs -n 1 basename | sed 's/.bw//' | tr '\n' ' ')
        colors=$(echo "$samples" | awk '{{for(i=1;i<=NF;i++) print "#66c2a5"}}')
        
        mkdir -p $(dirname {output.plot})
        plotProfile -m {input.matrix} -out {output.plot} \
            --plotType fill --colors $colors \
            --plotTitle "Ture_enhancer_from_peaks" \
            --perGroup \
            --samplesLabel $samples \
            --plotHeight 6 --plotWidth 15 > {log} 2>&1
        """

rule enhancer_integrated_1_cor_plot:
    input:
        rnaseq_deg = "result/RNAseq_pipline/4_DEseq/DEG_result.csv",
        chipseq_flag = "result/ChIPseq_pipline/4_peak_result/Peak_finished.txt",
        matrix = "result/Integrated/Enhancer/1_Integrated_genes/{group}/Ture_enh_gene_within_1mb_signal.tab",
        enh = "result/Integrated/Enhancer/1_Integrated_genes/{group}/Ture_enh_gene_within_1mb.bed"
    output:
        out_dir = directory("result/Integrated/Enhancer/1_Integrated_genes/{group}/correlation"),
        plots_dir = directory("result/Integrated/Enhancer/1_Integrated_genes/{group}/plots/correlation"),
        finished_flag = "result/Integrated/Enhancer/1_Integrated_genes/{group}/finished.txt"
    params:
        chipseq_anno_dir = "result/ChIPseq_pipline/4_peak_result/annotation",
        script = "workflow/scripts/Integrated/Enhancer/1_3_Enhancer_plot.R",
        yaml_file = "config/config.yaml"
    log:
        "logs/Integrated/Enhancer/1_Integrated_genes/{group}_cor_plot.log"
    shell:
        """
        mkdir -p {output.out_dir} {output.plots_dir}
        
        Rscript {params.script} \
            --RNAseq {input.rnaseq_deg} \
            --ChIPseq {params.chipseq_anno_dir} \
            --matrix {input.matrix} \
            --yaml {params.yaml_file} \
            --enh {input.enh} \
            --output {output.out_dir} \
            --picture {output.plots_dir} >> {log} 2>&1
            
        touch {output.finished_flag}
        """

rule enhancer_integrated_1_all:
    input:
        expand("result/Integrated/Enhancer/1_Integrated_genes/{group}/finished.txt", group=INTEGRATED_GROUPS)
