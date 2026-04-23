#!/usr/bin/env Rscript
rm(list = ls())

install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% c("org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "GOSemSim")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = repo)
      }
      BiocManager::install(package, update = FALSE)
    } else {
      install.packages(package, repos = repo, dependencies = TRUE)
    }
  }
  library(package, character.only = TRUE)
}

required_packages <- c(
  "dplyr", "org.Mm.eg.db", "org.Hs.eg.db", "ggplot2",
  "clusterProfiler", "enrichplot", "yaml", "optparse", "futile.logger"
)
lapply(required_packages, install_and_load)

option_list <- list(
  make_option(c("-i", "--indir"), type = "character", help="Input directory containing bed files and csv"),
  make_option(c("-o", "--outdir"), type = "character", help="Output directory"),
  make_option(c("-y", "--yaml"), type = "character", help="Config yaml file")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

indir <- opt$indir
outdir <- opt$outdir
yaml_file <- opt$yaml

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

if (!file.exists(yaml_file)) {
  stop("YAML config file not found.")
}

yaml_info <- yaml.load_file(yaml_file)
genome_type <- "mm" # default
if (!is.null(yaml_info$RNAseq) && !is.null(yaml_info$RNAseq$species)) {
  genome_type <- yaml_info$RNAseq$species
}

if (genome_type == "homo" || genome_type == "hg") {
  OrgDb <- org.Hs.eg.db
  species_kegg <- "hsa"
} else {
  OrgDb <- org.Mm.eg.db
  species_kegg <- "mmu"
}

csv_file <- file.path(indir, "Promoter_gene_all.csv")
if (!file.exists(csv_file)) {
  cat("Promoter_gene_all.csv not found in", indir, "\n")
  q()
}

gene_all <- read.csv(csv_file, stringsAsFactors = FALSE)

bed_files <- list.files(indir, pattern = "_common_up_peaks\\.bed$|_common_down_peaks\\.bed$", full.names = TRUE)

if (length(bed_files) == 0) {
  cat("No common peak bed files found.\n")
  q()
}

run_enrichment <- function(bed_file) {
  base_name <- sub("\\.bed$", "", basename(bed_file))
  cat("Running analysis for:", base_name, "\n")
  
  export_cytoscape <- function(enrich_res, out_file) {
    df <- as.data.frame(enrich_res)
    if (nrow(df) == 0) return()
    edges <- data.frame(Source=character(), Target=character(), stringsAsFactors=FALSE)
    for (i in 1:nrow(df)) {
      term <- df$Description[i]
      genes <- unlist(strsplit(as.character(df$geneID[i]), "/"))
      if (length(genes) > 0) {
        edges <- rbind(edges, data.frame(Source=term, Target=genes, stringsAsFactors=FALSE))
      }
    }
    write.csv(edges, out_file, row.names=FALSE, quote=FALSE)
  }
  
  bed_data <- tryCatch(read.table(bed_file, sep="\t", stringsAsFactors=FALSE, header=FALSE), error=function(e) NULL)
  if (is.null(bed_data) || nrow(bed_data) == 0) return(NULL)
  
  genes <- unique(bed_data$V4)
  if (length(genes) == 0) return(NULL)
  
  # Convert SYMBOL to ENTREZID
  gene_df <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb), error=function(e) NULL)
  if (is.null(gene_df) || nrow(gene_df) == 0) return(NULL)
  
  # For coloring cnetplot, we use the RNAseq log2FoldChange for these genes
  sub_gene_all <- gene_all[gene_all$gene %in% genes, ]
  geneList <- sub_gene_all$log2FoldChange
  names(geneList) <- sub_gene_all$gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  # GO Enrichment
  ego <- tryCatch({
    enrichGO(
      gene          = gene_df$ENTREZID,
      OrgDb         = OrgDb,
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )
  }, error=function(e) NULL)
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(outdir, paste0(base_name, "_GO.csv")), row.names = FALSE)
    
    p_bar <- barplot(ego, showCategory=10, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    ggsave(file.path(outdir, paste0(base_name, "_GO_barplot.pdf")), p_bar, width=10, height=8)
    
    p_dot <- dotplot(ego, showCategory=10, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    ggsave(file.path(outdir, paste0(base_name, "_GO_dotplot.pdf")), p_dot, width=10, height=8)
    
    # cnetplot: needs foldChange vector with SYMBOL names
    p_cnet <- cnetplot(ego, foldChange=geneList, showCategory=5)
    ggsave(file.path(outdir, paste0(base_name, "_GO_cnetplot.pdf")), p_cnet, width=12, height=10, bg="white")
    
    export_cytoscape(ego, file.path(outdir, paste0(base_name, "_GO_cytoscape.csv")))
  } else {
    cat("No significant GO terms for", base_name, "\n")
  }
  
  # KEGG Enrichment
  ekegg <- tryCatch({
    enrichKEGG(
      gene          = gene_df$ENTREZID,
      organism      = species_kegg,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH"
    )
  }, error=function(e) NULL)
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_read <- setReadable(ekegg, OrgDb = OrgDb, keyType="ENTREZID")
    write.csv(as.data.frame(ekegg_read), file.path(outdir, paste0(base_name, "_KEGG.csv")), row.names = FALSE)
    
    p_bar_kegg <- barplot(ekegg_read, showCategory=15)
    ggsave(file.path(outdir, paste0(base_name, "_KEGG_barplot.pdf")), p_bar_kegg, width=8, height=6)
    
    p_dot_kegg <- dotplot(ekegg_read, showCategory=15)
    ggsave(file.path(outdir, paste0(base_name, "_KEGG_dotplot.pdf")), p_dot_kegg, width=8, height=6)
    
    p_cnet_kegg <- cnetplot(ekegg_read, foldChange=geneList, showCategory=5)
    ggsave(file.path(outdir, paste0(base_name, "_KEGG_cnetplot.pdf")), p_cnet_kegg, width=10, height=8, bg="white")
    
    export_cytoscape(ekegg_read, file.path(outdir, paste0(base_name, "_KEGG_cytoscape.csv")))
  } else {
    cat("No significant KEGG pathways for", base_name, "\n")
  }
}

for (bed in bed_files) {
  run_enrichment(bed)
}

cat("KEGG and GO Enrichment Analysis Completed.\n")
