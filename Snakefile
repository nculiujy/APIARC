configfile: "config/config.yaml"

import pandas as pd
import re

def get_integrated_groups():
    try:
        df = pd.read_csv("config/ChIPseq_metadata.csv")
        # Extract group names by removing _repX from IP_name
        groups = set()
        for ip_name in df['IP_name'].dropna():
            group = re.sub(r'_rep\d+$', '', str(ip_name))
            groups.add(group)
        return list(groups)
    except Exception as e:
        return []

INTEGRATED_GROUPS = get_integrated_groups()

if config.get("RNAseq_modules", {}).get("1_download", True):
    include: "workflow/rules/RNAseq_pipline/1_download.smk"

if config.get("RNAseq_modules", {}).get("2_QC", True):
    include: "workflow/rules/RNAseq_pipline/2_QC.smk"

if config.get("RNAseq_modules", {}).get("3_RC", True):
    include: "workflow/rules/RNAseq_pipline/3_RC_pipline.smk"

if config.get("RNAseq_modules", {}).get("4_DEseq", True):
    include: "workflow/rules/RNAseq_pipline/4_DEseq.smk"

if config.get("ChIPseq_modules", {}).get("1_download", True):
    include: "workflow/rules/ChIPseq_pipline/1_download.smk"

if config.get("ChIPseq_modules", {}).get("2_QC", True):
    include: "workflow/rules/ChIPseq_pipline/2_QC.smk"

if config.get("ChIPseq_modules", {}).get("3_CC", True):
    include: "workflow/rules/ChIPseq_pipline/3_CC_pipline.smk"

if config.get("ChIPseq_modules", {}).get("4_peak_result", True):
    include: "workflow/rules/ChIPseq_pipline/4_peak_result.smk"

if config.get("Integrated_modules", {}).get("1_Integrated_genes", True):
    include: "workflow/rules/Integrated/Promoter/1_Intergrated_genes.smk"
    include: "workflow/rules/Integrated/Promoter/2_KEGG_GO_network.smk"
    include: "workflow/rules/Integrated/Promoter/3_Peaks_TF_gene_network.smk"

if config.get("Integrated_modules", {}).get("2_Enhancer_genes", True):
    include: "workflow/rules/Integrated/Enhancer/1_enhancer_Intergrated_genes.smk"
    include: "workflow/rules/Integrated/Enhancer/2_KEGG_GO_network.smk"
    include: "workflow/rules/Integrated/Enhancer/3_Peaks_TF_gene_network.smk"

rule all:
    input:
        "result/RNAseq_pipline/1_Rawdata/finished.txt" if config.get("RNAseq_modules", {}).get("1_download", True) and not config.get("RNAseq_modules", {}).get("2_QC", True) and not config.get("RNAseq_modules", {}).get("3_RC", True) and not config.get("RNAseq_modules", {}).get("4_DEseq", True) else [],
        "result/RNAseq_pipline/2_Cleandata" if config.get("RNAseq_modules", {}).get("2_QC", True) and not config.get("RNAseq_modules", {}).get("3_RC", True) and not config.get("RNAseq_modules", {}).get("4_DEseq", True) else [],
        "result/RNAseq_pipline/3_RC_pipline/Quant_finished.txt" if config.get("RNAseq_modules", {}).get("3_RC", True) and not config.get("RNAseq_modules", {}).get("4_DEseq", True) else [],
        "result/RNAseq_pipline/4_DEseq/DEseq_finished.txt" if config.get("RNAseq_modules", {}).get("4_DEseq", True) else [],
        "result/ChIPseq_pipline/1_Rawdata/finished.txt" if config.get("ChIPseq_modules", {}).get("1_download", True) and not config.get("ChIPseq_modules", {}).get("2_QC", True) and not config.get("ChIPseq_modules", {}).get("3_CC", True) and not config.get("ChIPseq_modules", {}).get("4_peak_result", True) else [],
        "result/ChIPseq_pipline/2_Cleandata" if config.get("ChIPseq_modules", {}).get("2_QC", True) and not config.get("ChIPseq_modules", {}).get("3_CC", True) and not config.get("ChIPseq_modules", {}).get("4_peak_result", True) else [],
        "result/ChIPseq_pipline/3_CC_pipline/BW_finished.txt" if config.get("ChIPseq_modules", {}).get("3_CC", True) and not config.get("ChIPseq_modules", {}).get("4_peak_result", True) and not config.get("Integrated_modules", {}).get("1_Integrated_genes", True) else [],
        "result/ChIPseq_pipline/4_peak_result/Peak_finished.txt" if config.get("ChIPseq_modules", {}).get("4_peak_result", True) and not config.get("Integrated_modules", {}).get("1_Integrated_genes", True) else [],
        "result/Integrated/Promoter/1_Integrated_genes/finished.txt",
        "result/Integrated/Promoter/2_KEGG_GO_network/finished.txt",
        expand("result/Integrated/Promoter/3_Peaks_TF_gene_network/{group}/finished.txt", group=INTEGRATED_GROUPS) if config.get("Integrated_modules", {}).get("1_Integrated_genes", True) else [],
        expand("result/Integrated/Enhancer/1_Integrated_genes/{group}/finished.txt", group=INTEGRATED_GROUPS) if config.get("Integrated_modules", {}).get("2_Enhancer_genes", True) else [],
        expand("result/Integrated/Enhancer/2_KEGG_GO_network/{group}/finished.txt", group=INTEGRATED_GROUPS) if config.get("Integrated_modules", {}).get("2_Enhancer_genes", True) else [],
        expand("result/Integrated/Enhancer/3_Peaks_TF_gene_network/{group}/finished.txt", group=INTEGRATED_GROUPS) if config.get("Integrated_modules", {}).get("2_Enhancer_genes", True) else []
