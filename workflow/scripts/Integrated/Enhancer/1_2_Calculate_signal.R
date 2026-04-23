#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-m", "--matrix1"), type="character"),
  make_option(c("-n", "--matrix2"), type="character"),
  make_option(c("-o", "--output"), type="character"),
  make_option(c("-s", "--names"), type="character"),
  make_option(c("-a", "--anno"), type="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$matrix1) | is.null(opt$matrix2) | is.null(opt$output) | is.null(opt$names) | is.null(opt$anno)) quit(status=1)

col_names <- strsplit(opt$names, ",")[[1]]
n_cols_per_group <- 40 # 2000bp / 50bp binSize = 40 bins

process_chipseq_data_by_names <- function(file_path, col_names, n_per_group=40, skip_lines=3) {
  if (!file.exists(file_path)) quit(status=1)
  ChIP <- read.table(file_path, header=FALSE, stringsAsFactors=FALSE, skip=skip_lines)
  ChIP[is.na(ChIP)] <- 0
  total_expected <- n_per_group * length(col_names)
  if (ncol(ChIP) < total_expected) quit(status=1)
  group_means <- lapply(seq_along(col_names), function(i) {
    start_col <- (i-1)*n_per_group + 1
    end_col <- i*n_per_group
    rowMeans(ChIP[,start_col:end_col, drop=FALSE])
  })
  ChIP_mean <- do.call(cbind, group_means)
  colnames(ChIP_mean) <- col_names
  as.data.frame(ChIP_mean)
}

enh_anno <- read.table(opt$anno, header=FALSE, stringsAsFactors=FALSE)
colnames(enh_anno) <- c("Chr", "start", "end","enhancer",  "SYMBOL", "distance")
ChIPseqData_mean <- process_chipseq_data_by_names(opt$matrix1, col_names, n_cols_per_group)
ChIPseqData_random_mean <- process_chipseq_data_by_names(opt$matrix2, col_names, n_cols_per_group)

quantile_threshold <- 0.95
thresholds <- sapply(col_names, function(name) quantile(ChIPseqData_random_mean[[name]], quantile_threshold))

selected_idx <- which(apply(ChIPseqData_mean, 1, function(row) {
  all(row > thresholds)
}))

if (nrow(enh_anno) != nrow(ChIPseqData_mean)) quit(status=1)
enh_allinfo <- cbind(enh_anno, ChIPseqData_mean)
filtered_enh <- enh_allinfo[selected_idx, ]
Ture_enh_bedfile <- filtered_enh[, c("Chr", "start", "end", "enhancer","SYMBOL")]

write.table(
  Ture_enh_bedfile,
  file = opt$output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)