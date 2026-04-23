rule integrated_2_kegg_go:
    input:
        finished_flag="result/Integrated/Promoter/1_Integrated_genes/finished.txt"
    output:
        finished_flag="result/Integrated/Promoter/2_KEGG_GO_network/finished.txt"
    params:
        promoter_dir="result/Integrated/Promoter/1_Integrated_genes",
        out_dir="result/Integrated/Promoter/2_KEGG_GO_network",
        script="workflow/scripts/Integrated/Promoter/2_KEGG_GO_network.R",
        yaml_file="config/config.yaml"
    log:
        "logs/Integrated/Promoter/2_KEGG_GO_network.log"
    shell:
        """
        Rscript {params.script} \
            --indir {params.promoter_dir} \
            --outdir {params.out_dir} \
            --yaml {params.yaml_file} > {log} 2>&1
            
        touch {output.finished_flag}
        """
