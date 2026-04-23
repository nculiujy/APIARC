#!/usr/bin/env Rscript


if (!require("optparse")) install.packages("optparse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!require("yaml")) install.packages("yaml", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!require("RColorBrewer")) install.packages("RColorBrewer", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!require("ChIPseeker")) BiocManager::install("ChIPseeker")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")

library(optparse)
library(yaml)
library(ChIPseeker)
library(ggplot2)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(graphics)
library(RColorBrewer)

rm(list = ls())


option_list <- list(
  make_option(c("-p", "--peakfile"), type = "character", help = "Path to the input Peak file (.narrowPeak or .bed)"),
  make_option(c("-o", "--output"), type = "character", help = "Path to the output directory", default = NULL),
  make_option(c("-n", "--name"), type = "character", help = "Prefix name for output files", default = "PeakAnno"),
  make_option(c("-c", "--config"), type = "character", help = "Path to the config.yaml file", default = "config/config.yaml")
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (is.null(opt$peakfile)) {
  stop("--peakfile must be provided.")
}


if (is.null(opt$output)) {
  opt$output <- getwd()
  cat("No output path provided, using current working directory:", opt$output, "\n")
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (!file.exists(opt$peakfile)) {
  cat("Warning: The input Peak file does not exist. MACS2 may have been skipped or failed for this sample.\n")
  

  peak_anno_file <- file.path(opt$output, paste0(opt$name, "_peak_anno.csv"))
  write.csv(data.frame(Message="Peak file not found"), peak_anno_file, row.names = FALSE)
  
  pdf_file <- file.path(opt$output, paste0(opt$name, "_pie_bp±2000.pdf"))
  pdf(pdf_file, width = 7.5, height = 7.5)
  plot.new()
  text(0.5, 0.5, "Peak file not found")
  dev.off()
  
  cat("Annotation skipped due to missing peak file. Empty outputs generated.\n")
  quit(save = "no", status = 0)
}

if (!file.exists(opt$config)) {
  stop("The config file does not exist.")
}


config <- yaml.load_file(opt$config)


species <- "mm"  # default
if (!is.null(config$species)) {
  species <- config$species
}

cat("Detected species from config:", species, "\n")


txdb_name <- NULL
orgdb_name <- NULL

if (species == "TAIR") {
  txdb_name <- "TxDb.Athaliana.BioMart.plantsmart28"
  orgdb_name <- "org.At.tair.db"
} else if (species == "mm") {
  txdb_name <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
  orgdb_name <- "org.Mm.eg.db"
} else if (species == "homo") {
  txdb_name <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
  orgdb_name <- "org.Hs.eg.db"
} else {
  stop("Unsupported species:", species)
}

cat("Using TxDb:", txdb_name, "\n")
cat("Using OrgDb:", orgdb_name, "\n")


if (!require(txdb_name, character.only = TRUE)) {
  BiocManager::install(txdb_name, ask = FALSE)
  library(txdb_name, character.only = TRUE)
}
if (!require(orgdb_name, character.only = TRUE)) {
  BiocManager::install(orgdb_name, ask = FALSE)
  library(orgdb_name, character.only = TRUE)
}

txdb <- get(txdb_name)
annoDb_name <- orgdb_name


options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)


if (file.info(opt$peakfile)$size == 0) {
  cat("Warning: The input Peak file is empty (size 0). MACS2 may not have found any peaks under the current cutoff.\n")
  

  peak_anno_file <- file.path(opt$output, paste0(opt$name, "_peak_anno.csv"))
  write.csv(data.frame(Message="No peaks found"), peak_anno_file, row.names = FALSE)
  
  pdf_file <- file.path(opt$output, paste0(opt$name, "_pie_bp±2000.pdf"))
  pdf(pdf_file, width = 7.5, height = 7.5)
  plot.new()
  text(0.5, 0.5, "No peaks found by MACS2")
  dev.off()
  
  cat("Annotation skipped due to empty peak file. Empty outputs generated.\n")
  quit(save = "no", status = 0)
}


peak <- tryCatch({
  readPeakFile(opt$peakfile)
}, error = function(e) {
  cat("Warning: readPeakFile failed. The peak file might be malformed or lack necessary columns.\n")
  cat("Error message:", conditionMessage(e), "\n")
  
  peak_anno_file <- file.path(opt$output, paste0(opt$name, "_peak_anno.csv"))
  write.csv(data.frame(Message="Peak file parsing failed"), peak_anno_file, row.names = FALSE)
  
  pdf_file <- file.path(opt$output, paste0(opt$name, "_pie_bp±2000.pdf"))
  pdf(pdf_file, width = 7.5, height = 7.5)
  plot.new()
  text(0.5, 0.5, "Peak file parsing failed")
  dev.off()
  
  quit(save = "no", status = 0)
})


peakAnno <- annotatePeak(
  peak,
  tssRegion = c(-2000, 2000),
  TxDb = txdb,
  annoDb = annoDb_name
)


annotationData <- as.data.frame(peakAnno)


total_peaks <- nrow(annotationData)
cat("Total peaks:", total_peaks, "\n")


annotationData$annotation_category <- ifelse(grepl("Promoter", annotationData$annotation), "Promoter",
                                    ifelse(grepl("Intron", annotationData$annotation), "Intron",
                                           ifelse(grepl("Exon", annotationData$annotation), "Exon",
                                                  ifelse(grepl("3' UTR", annotationData$annotation), "3' UTR",
                                                         ifelse(grepl("Downstream", annotationData$annotation) | grepl("Distal Intergenic", annotationData$annotation), "Intergenic", annotationData$annotation)))))


annotationCounts <- table(annotationData$annotation_category)
annotationCounts <- as.data.frame(annotationCounts)
colnames(annotationCounts) <- c("Category", "Count")


annotationCounts <- annotationCounts[annotationCounts$Category %in% c("Promoter", "Intergenic", "Exon", "Intron", "3' UTR"), ]


if(sum(annotationCounts$Count) > 0) {
  annotationCounts$Percentage <- annotationCounts$Count / sum(annotationCounts$Count) * 100
} else {
  annotationCounts$Percentage <- 0
}


pdf_file <- file.path(opt$output, paste0(opt$name, "_pie_bp±2000.pdf"))
pdf(pdf_file, width = 7.5, height = 7.5)

if (nrow(annotationCounts) > 0 && sum(annotationCounts$Count) > 0) {

  labels <- annotationCounts$Category
  values <- annotationCounts$Count
  percentages <- paste(format(round(annotationCounts$Percentage, 2), nsmall = 2), "%", sep = "")
  label_text <- paste(labels, "\n(", values, ", ", percentages, ")", sep = "")
  colors <- brewer.pal(max(3, length(annotationCounts$Category)), "Set3")[1:length(annotationCounts$Category)]

  pie(values, labels = label_text, col = colors, border = "white", main = paste(opt$name, "- Total peaks:", total_peaks))
} else {
  plot.new()
  text(0.5, 0.5, "No peaks matched required categories")
}
dev.off()


peak_anno_file <- file.path(opt$output, paste0(opt$name, "_peak_anno.csv"))
write.csv(annotationData, peak_anno_file, row.names = FALSE)


promoter_peaks <- annotationData[grepl("Promoter", annotationData$annotation_category), ]
promoter_file <- file.path(opt$output, paste0(opt$name, "_promoter_peaks.csv"))
write.csv(promoter_peaks, promoter_file, row.names = FALSE)

cat("Annotation finished successfully. Output saved to:", opt$output, "\n")