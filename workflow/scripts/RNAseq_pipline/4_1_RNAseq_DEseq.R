#!/usr/bin/env Rscript

if(!require("yaml")) install.packages("yaml", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("optparse")) install.packages("optparse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
library(DESeq2)
library(yaml)
library(optparse)

.libPaths(c("/usr/local/R/lib64/R/library", .libPaths()))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to the gene count matrix CSV file"),
  make_option(c("-c", "--config"), type = "character", help = "Path to the config JSON string"),
  make_option(c("-o", "--output"), type = "character", help = "Path to the output result CSV file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (is.null(opt$input) || is.null(opt$config) || is.null(opt$output)) {
  stop("Error: Input, config, and output file paths must be provided. Use --input, --config, and --output options.")
}

cat("[INFO] file: ", opt$input, "\n")
cat("[INFO] file: ", opt$config, "\n")
cat("[INFO] file: ", opt$output, "\n")


cat("[INFO] Reading gene count matrix...\n")
if (file.info(opt$input)$size == 0) {
  cat("[WARNING] Gene count matrix is empty! Creating dummy output and exiting gracefully.\n")
  file.create(opt$output)
  quit(status = 0)
}
database_all <- read.table(file = opt$input, sep = ",", header = TRUE)
cat("[INFO] Gene count matrix content:\n")
print(head(database_all))


cat("[INFO] file...\n")
config <- yaml.load(opt$config)
cat("[INFO] Config content:\n")
print(config)


cat("[INFO] Processing expression matrix...\n")
expdata <- database_all[, -1]
rownames(expdata) <- database_all[, 1]
cat("[INFO] Expression matrix dimension: ", dim(expdata), "\n")


cat("[INFO] Extracting sample names and conditions...\n")
sample_names <- names(config)
conditions <- factor(unlist(config), levels = c("T", "P"))
cat("[INFO] Sample:\n")
print(sample_names)
cat("[INFO] Sample:\n")
print(conditions)


if (!all(sample_names %in% colnames(expdata))) {
  cat("[ERROR] Sample names do not match expression matrix columns！\n")
  cat("Expression matrix columns:\n")
  print(colnames(expdata))
  cat("Sample names in config:\n")
  print(sample_names)
  stop("Sample，file！")
}


cat("[INFO] Building colData dataframe...\n")
coldata <- data.frame(conditions = conditions, row.names = sample_names)
cat("[INFO] colData content:\n")
print(coldata)
expdata <- expdata[, sample_names]

# Creating DESeq dataset
cat("[INFO] Removing rows with NA values...\n")
expdata <- expdata[complete.cases(expdata), ]
cat("[INFO] Creating DESeq dataset...\n")
dds <- DESeqDataSetFromMatrix(countData = expdata, colData = coldata, design = ~ conditions)


cat("[INFO] Running DESeq2 analysis...\n")
dds <- tryCatch({
  DESeq(dds)
}, error = function(e) {
  cat("[WARNING] DESeq2 failed (likely due to zero counts across all genes). Exiting gracefully.\n")
  file.create(opt$output)
  quit(status = 0)
})


cat("[INFO] analysis...\n")
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
row.names(resdata) <- resdata[, 1]

names(resdata)[names(resdata) == "Row.names"] <- "SYMBOL"
resdata$gene_id <- sub("\\|[^|]*$", "", resdata$SYMBOL)
resdata$SYMBOL <- sub("^.*\\|", "", resdata$SYMBOL)
resdata <- resdata[!is.na(resdata$padj), ]

cat("[INFO] analysiscontent:\n")
print(head(resdata))


cat("[INFO] file...\n")
write.csv(resdata, opt$output, row.names = FALSE)

cat("[INFO] analysis，: ", opt$output, "\n")

# Rscript RNAseq_DEseq.R --input /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632/mRNA/gene_count_matrix.csv --config "{'mESC_0hDOX_re1': 'T', 'mESC_0hDOX_re2': 'T', 'mESC_24hDOX_re1': 'P', 'mESC_24hDOX_re2': 'P'}" --output /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632/mRNA/DEG_result.csv

