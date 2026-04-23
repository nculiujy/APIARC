
rule enhancer_integrated_2_kegg_go:
    input:
        enhancer_dir="result/Integrated/Enhancer/1_Integrated_genes/{group}/correlation",
        finished_flag="result/Integrated/Enhancer/1_Integrated_genes/{group}/finished.txt"
    output:
        finished_flag="result/Integrated/Enhancer/2_KEGG_GO_network/{group}/finished.txt"
    params:
        out_dir="result/Integrated/Enhancer/2_KEGG_GO_network/{group}",
        script="workflow/scripts/Integrated/Enhancer/2_KEGG_GO_network.R",
        yaml_file="config/config.yaml"
    log:
        "logs/Integrated/Enhancer/2_KEGG_GO_network/{group}_KEGG_GO.log"
    shell:
        """
        mkdir -p {params.out_dir}
        
        Rscript {params.script} \
            --indir {input.enhancer_dir} \
            --outdir {params.out_dir} \
            --yaml {params.yaml_file} > {log} 2>&1
            
        touch {output.finished_flag}
        """

rule enhancer_integrated_2_all:
    input:
        expand("result/Integrated/Enhancer/2_KEGG_GO_network/{group}/finished.txt", group=INTEGRATED_GROUPS)