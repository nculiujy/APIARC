#!/usr/bin/env Rscript

rm(list = ls())

install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% installed.packages()) {
      suppressWarnings(library(package, character.only = TRUE))
    } else {
      if (package %in% c("org.Mm.eg.db")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = repo)
        }
        BiocManager::install(package, ask = FALSE)
      } else {
        install.packages(package, repos = repo)
      }
    }
  }
  suppressWarnings(suppressPackageStartupMessages(library(package, character.only = TRUE)))
}

required_packages <- c(
  "dplyr",
  "ComplexHeatmap",
  "circlize",
  "optparse"
)
lapply(required_packages, install_and_load)

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(optparse)

option_list <- list(
  make_option(c("-r", "--RNAseq"), type = "character"),
  make_option(c("-c", "--TF_gene"), type = "character"),
  make_option(c("-o", "--output"), type = "character"),
  make_option(c("-p", "--picture"), type = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$RNAseq) | is.null(opt$TF_gene) | is.null(opt$output) | is.null(opt$picture)) {
  cat("Missing arguments\n")
  quit(status=1)
}

RNAseq_gene_file <- opt$RNAseq
TF_gene_file <- opt$TF_gene
output_folder <- opt$output
heatmap_pdf <- opt$picture

cyto_file <- sub("\\.txt$", "_cytoscape.csv", output_folder)

if (file.info(TF_gene_file)$size == 0) {
  cat("TF gene list file is empty:", TF_gene_file, "\n")

  file.create(output_folder)
  file.create(cyto_file)
  pdf(heatmap_pdf, width = 4, height = 4)
  plot.new()
  text(0.5, 0.5, "No enriched motifs found", cex = 1.5)
  dev.off()
  quit(status = 0)
}

TF_motif_data <- read.table(TF_gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(TF_motif_data) <- c("Ensembl", "TF", "pvalue")
TF_unique <- TF_motif_data[!duplicated(TF_motif_data$TF), ]
TF_filtered <- TF_unique[TF_unique$pvalue < 0.05, ]
TF_all_filtered <- TF_motif_data[TF_motif_data$pvalue < 0.05, ]

RNAseq_gene <- read.csv(RNAseq_gene_file, stringsAsFactors = FALSE)
RNAseq_gene <- RNAseq_gene[!is.na(RNAseq_gene$padj), ]
RNAseq_gene$gene_id_noversion <- sub("\\..*", "", RNAseq_gene$gene_id)

non_sample_columns <- c("SYMBOL", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_id", "gene_id_noversion")
samples <- setdiff(colnames(RNAseq_gene), non_sample_columns)

TF_gene_data_merged <- merge(
  TF_filtered,
  RNAseq_gene[, c("gene_id_noversion", "log2FoldChange", "padj", samples, "SYMBOL")],
  by.x = "Ensembl",
  by.y = "gene_id_noversion",
  all.x = TRUE
)

# Export full TF-gene Cytoscape network
cyto_file <- sub("\\.txt$", "_cytoscape.csv", output_folder)
TF_network_all <- merge(
  TF_all_filtered,
  RNAseq_gene[, c("gene_id_noversion", "SYMBOL", "log2FoldChange", "padj")],
  by.x = "Ensembl",
  by.y = "gene_id_noversion",
  all.x = TRUE
)
TF_network_all <- na.omit(TF_network_all)
if(nrow(TF_network_all) > 0) {
  cyto_df <- data.frame(
    Source = TF_network_all$TF,
    Target = TF_network_all$SYMBOL,
    TF_pvalue = TF_network_all$pvalue,
    Gene_log2FC = TF_network_all$log2FoldChange,
    Gene_padj = TF_network_all$padj,
    stringsAsFactors = FALSE
  )
  write.csv(cyto_df, cyto_file, row.names = FALSE, quote = FALSE)
} else {
  file.create(cyto_file)
}

if (nrow(TF_gene_data_merged) == 0) {
  cat("No merged TF data available, exiting gracefully\n")
  file.create(heatmap_pdf)
  quit(status = 0)
}
TF_gene_data_merged <- na.omit(TF_gene_data_merged)
if (nrow(TF_gene_data_merged) == 0) {
  cat("No merged TF data after NA omit, exiting gracefully\n")
  file.create(heatmap_pdf)
  quit(status = 0)
}

expData <- TF_gene_data_merged[, samples, drop = FALSE]
rownames(expData) <- TF_gene_data_merged$TF
expData <- na.omit(expData)
if (nrow(expData) < 2 || ncol(expData) == 0) {
  cat("Not enough expData, exiting gracefully\n")
  file.create(heatmap_pdf)
  quit(status = 0)
}

valscaling <- function(x) { x / sum(x) }
zscore <- function(x) { (x - mean(x)) / sd(x) }

expData_scaled <- apply(expData, 2, valscaling)
if (is.null(dim(expData_scaled))) expData_scaled <- t(as.matrix(expData_scaled))
expData_scaled <- t(apply(expData_scaled, 1, zscore))
if (is.null(dim(expData_scaled))) expData_scaled <- t(as.matrix(expData_scaled))
expData_scaled <- as.matrix(expData_scaled)
rownames(expData_scaled) <- rownames(expData)

if (nrow(expData_scaled) < 2) {
  cat("Not enough expData_scaled, exiting gracefully\n")
  file.create(heatmap_pdf)
  quit(status = 0)
}

TF_gene_data_merged <- TF_gene_data_merged[TF_gene_data_merged$TF %in% rownames(expData_scaled), ]
UpTFs <- rownames(expData_scaled)[which(TF_gene_data_merged$log2FoldChange > 1)]
DownTFs <- rownames(expData_scaled)[which(TF_gene_data_merged$log2FoldChange < -1)]
matchIndexes <- match(c(UpTFs, DownTFs), rownames(expData_scaled))

f1 <- colorRamp2(seq(min(expData_scaled, na.rm=TRUE), max(expData_scaled, na.rm=TRUE), length = 4),
                 c("#CCFFFF", "#99CCFF", "#FF9966", "#CC3333"))

pdf(heatmap_pdf, width = 4, height = 10)
# Make sure matchIndexes has no NA and is within bounds
valid_matches <- matchIndexes[!is.na(matchIndexes)]
if (length(valid_matches) > 0) {
  ha <- rowAnnotation(foo = anno_mark(at = valid_matches, labels = rownames(expData_scaled)[valid_matches]))
} else {
  ha <- NULL
}

if (nrow(expData_scaled) >= 2) {
  km_num <- min(2, nrow(expData_scaled))
  if (km_num >= 2 && nrow(expData_scaled) > 2) {
    if (!is.null(ha)) {
      hm <- Heatmap(expData_scaled, name = "z-score",
              cluster_columns = FALSE,
              show_row_names = FALSE,
              right_annotation = ha,
              col = f1,
              row_names_gp = gpar(fontsize = 8),
              row_km = km_num)
    } else {
      hm <- Heatmap(expData_scaled, name = "z-score",
              cluster_columns = FALSE,
              show_row_names = FALSE,
              col = f1,
              row_names_gp = gpar(fontsize = 8),
              row_km = km_num)
    }
  } else {
    if (!is.null(ha)) {
      hm <- Heatmap(expData_scaled, name = "z-score",
              cluster_columns = FALSE,
              show_row_names = FALSE,
              right_annotation = ha,
              col = f1,
              row_names_gp = gpar(fontsize = 8))
    } else {
      hm <- Heatmap(expData_scaled, name = "z-score",
              cluster_columns = FALSE,
              show_row_names = FALSE,
              col = f1,
              row_names_gp = gpar(fontsize = 8))
    }
  }
  draw(hm)
}
dev.off()

