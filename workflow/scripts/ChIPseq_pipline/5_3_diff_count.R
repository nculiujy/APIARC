#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

option_list <- list(
  make_option("--counts", type = "character"),
  make_option("--metadata", type = "character"),
  make_option("--support", type = "character"),
  make_option("--mark", type = "character"),
  make_option("--reference", type = "character"),
  make_option("--treatment", type = "character"),
  # optparse keeps the dash in the name, so opt$min_support would be NULL
  # and the support filter would drop every region.
  make_option("--min-support", type = "integer", default = 2L, dest = "min_support"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--lfc", type = "double", default = 1.0),
  make_option("--out", type = "character"),
  make_option("--norm-out", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (required in c("counts", "metadata", "support", "mark", "reference", "treatment", "out")) {
  if (is.null(opt[[required]])) stop("missing --", required)
}

meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)
if ("role" %in% names(meta)) {
  # Input rows are controls, not contrast members.
  meta <- meta[tolower(trimws(meta$role)) == "ip", , drop = FALSE]
}
if ("input_sample" %in% names(meta)) {
  # Role-aware schema: sample_name names the group and condition is only the
  # replicate tag, so the contrast is between sample_names.
  meta$replicate <- as.integer(sub("^rep", "", meta$condition))
  meta$condition <- meta$sample_name
} else {
  meta <- meta[meta$mark == opt$mark, ]
}
meta <- meta[meta$condition %in% c(opt$reference, opt$treatment), ]
meta <- meta[order(meta$condition, meta$replicate), ]
if (nrow(meta) == 0) stop("no samples for mark ", opt$mark)

counts <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
support <- read.delim(opt$support, check.names = FALSE, stringsAsFactors = FALSE)
counts <- counts[counts$region_id %in% support$region_id[support$n_support >= opt$min_support], ]
if (nrow(counts) == 0) stop("no regions pass the replicate support filter")

rownames(counts) <- counts$region_id
mat <- as.matrix(counts[, meta$sample, drop = FALSE])
mat <- mat[rowSums(mat) > 0, , drop = FALSE]
if (ncol(mat) < 2 || nrow(mat) < 10) stop("not enough data for DESeq2")

condition <- factor(meta$condition, levels = c(opt$reference, opt$treatment))
coldata <- data.frame(condition = condition, row.names = meta$sample)

dds <- DESeqDataSetFromMatrix(countData = round(mat), colData = coldata, design = ~condition)
dds <- DESeq(dds, quiet = TRUE)

res <- results(dds, contrast = c("condition", opt$treatment, opt$reference), alpha = opt$fdr)
shrunk <- tryCatch(
  lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm", quiet = TRUE),
  error = function(e) NULL
)
if (!is.null(shrunk)) {
  res$log2FoldChange <- shrunk$log2FoldChange
  res$lfcSE <- shrunk$lfcSE
  shrinkage <- "apeglm"
} else {
  shrinkage <- "none"
}

out <- data.frame(
  region_id = rownames(res),
  baseMean = res$baseMean,
  log2FC = res$log2FoldChange,
  lfcSE = res$lfcSE,
  stat = res$stat,
  pvalue = res$pvalue,
  FDR = res$padj,
  stringsAsFactors = FALSE
)
out <- out[!is.na(out$FDR), ]
out$direction <- ifelse(out$log2FC > 0, "up", "down")
out$significant <- out$FDR < opt$fdr & abs(out$log2FC) >= opt$lfc
out$method <- "count"
out$shrinkage <- shrinkage

write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(opt[["norm-out"]])) {
  norm <- counts(dds, normalized = TRUE)
  write.table(
    data.frame(region_id = rownames(norm), norm, check.names = FALSE),
    opt[["norm-out"]], sep = "\t", quote = FALSE, row.names = FALSE
  )
}

cat(sprintf(
  "[diff-count] %s %s_vs_%s: %d regions tested, %d significant (FDR<%.3g, |log2FC|>=%.2f), shrinkage=%s\n",
  opt$mark, opt$treatment, opt$reference, nrow(out), sum(out$significant, na.rm = TRUE),
  opt$fdr, opt$lfc, shrinkage
))
