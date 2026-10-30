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
  # Dashed names need an explicit dest: optparse keeps the dash, so
  # opt$min_support and opt$max_false_positive would both read NULL.
  make_option("--min-support", type = "integer", default = 2L, dest = "min_support"),
  make_option("--permutations", type = "integer", default = 20L),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--lfc", type = "double", default = 1.0),
  make_option("--max-false-positive", type = "double", default = 5,
              dest = "max_false_positive"),
  make_option("--seed", type = "integer", default = 20261023L),
  make_option("--out", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (required in c("counts", "metadata", "support", "mark", "reference", "treatment", "out")) {
  if (is.null(opt[[required]])) stop("missing --", required)
}
set.seed(opt$seed)

counts_all <- read.delim(opt$counts, check.names = FALSE, stringsAsFactors = FALSE)
support <- read.delim(opt$support, check.names = FALSE, stringsAsFactors = FALSE)
counts_all <- counts_all[counts_all$region_id %in% support$region_id[support$n_support >= opt$min_support], ]
rownames(counts_all) <- counts_all$region_id

meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)
if ("role" %in% names(meta)) {
  # Input rows are controls, not contrast members.
  meta <- meta[tolower(trimws(meta$role)) == "ip", , drop = FALSE]
}
if ("input_sample" %in% names(meta)) {
  # Role-aware schema: sample_name names the group and condition is only the
  # replicate tag.
  meta$replicate <- as.integer(sub("^rep", "", meta$condition))
  meta$condition <- meta$sample_name
} else {
  meta <- meta[meta$mark == opt$mark, ]
}
meta <- meta[meta$condition %in% c(opt$reference, opt$treatment), ]

run_model <- function(mat, labels) {
  labels <- factor(labels, levels = unique(labels))
  if (nlevels(labels) < 2 || min(table(labels)) < 2) return(NULL)
  coldata <- data.frame(condition = labels, row.names = colnames(mat))
  dds <- tryCatch(
    DESeqDataSetFromMatrix(round(mat), coldata, design = ~condition),
    error = function(e) NULL
  )
  if (is.null(dds)) return(NULL)
  dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) NULL)
  if (is.null(dds)) return(NULL)
  res <- tryCatch(
    results(dds, contrast = c("condition", levels(labels)[2], levels(labels)[1]), alpha = opt$fdr),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)
  sum(res$padj < opt$fdr & abs(res$log2FoldChange) >= opt$lfc, na.rm = TRUE)
}

rows <- list()
ip_samples <- meta$sample[meta$sample %in% colnames(counts_all)]
ip_mat <- as.matrix(counts_all[, ip_samples, drop = FALSE])
ip_mat <- ip_mat[rowSums(ip_mat) > 0, , drop = FALSE]

n_ref <- sum(meta$condition == opt$reference)
n_treat <- sum(meta$condition == opt$treatment)
labels <- c(rep(opt$reference, n_ref), rep(opt$treatment, n_treat))
n_permutations <- min(opt$permutations, choose(length(labels), n_treat))
perm_counts <- c()
if (n_permutations >= 3 && nrow(ip_mat) > 0) {
  for (i in seq_len(n_permutations)) {
    shuffled <- sample(labels)
    observed <- run_model(ip_mat, shuffled)
    if (!is.null(observed)) perm_counts <- c(perm_counts, observed)
  }
}
if (length(perm_counts) > 0) {
  rows[[length(rows) + 1]] <- data.frame(
    control = "condition_label_permutation",
    n_replicates = length(perm_counts),
    median_false_positive = median(perm_counts),
    max_false_positive = max(perm_counts),
    threshold = opt$max_false_positive,
    status = ifelse(max(perm_counts) <= opt$max_false_positive, "PASS", "FAIL"),
    note = "empirical null for the number of significant regions"
  )
} else {
  rows[[length(rows) + 1]] <- data.frame(
    control = "condition_label_permutation",
    n_replicates = 0L,
    median_false_positive = NA_real_,
    max_false_positive = NA_real_,
    threshold = opt$max_false_positive,
    status = "SKIPPED",
    note = "not enough replicates to permute labels (need >= 3 distinct permutations)"
  )
}

input_cols <- grep("^Input_", colnames(counts_all), value = TRUE)
if (length(input_cols) >= 4) {
  half <- floor(length(input_cols) / 2)
  pseudo <- c(rep("A", half), rep("B", length(input_cols) - half))
  input_mat <- as.matrix(counts_all[, input_cols, drop = FALSE])
  input_mat <- input_mat[rowSums(input_mat) > 0, , drop = FALSE]
  observed <- run_model(input_mat, pseudo)
  rows[[length(rows) + 1]] <- data.frame(
    control = "input_vs_input",
    n_replicates = 1L,
    median_false_positive = observed,
    max_false_positive = observed,
    threshold = opt$max_false_positive,
    status = ifelse(!is.null(observed) && observed <= opt$max_false_positive, "PASS", "FAIL"),
    note = "unexpected significant regions indicate a calibration problem"
  )
} else {
  rows[[length(rows) + 1]] <- data.frame(
    control = "input_vs_input",
    n_replicates = 0L,
    median_false_positive = NA_real_,
    max_false_positive = NA_real_,
    threshold = opt$max_false_positive,
    status = "SKIPPED",
    note = sprintf("only %d Input libraries available (need >= 4)", length(input_cols))
  )
}

out <- do.call(rbind, rows)
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[negative-control] %s: %s\n", opt$mark, paste(out$status, out$control, sep = ":", collapse = ", ")))
