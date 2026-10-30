#!/usr/bin/env Rscript

.required <- c("DESeq2", "yaml", "optparse")
.missing <- .required[!vapply(.required, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing) > 0) {
  stop(
    "missing R package(s): ", paste(.missing, collapse = ", "),
    "\nThey must come from the pinned conda environment; run snakemake with --use-conda.",
    call. = FALSE
  )
}
library(DESeq2)
library(yaml)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to the gene count matrix CSV file"),
  make_option(c("-c", "--config"), type = "character",
              help = "Optional YAML/JSON mapping of sample_name -> group; overrides --metadata"),
  make_option(c("-m", "--metadata"), type = "character",
              help = "RNA-seq metadata CSV with sample_name and group columns"),
  make_option(c("-o", "--output"), type = "character", help = "Path to the output result CSV file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (is.null(opt$input) || is.null(opt$output) ||
    (is.null(opt$config) && is.null(opt$metadata))) {
  stop("Error: --input, --output and one of --config/--metadata are required.")
}

cat("[INFO] file: ", opt$input, "\n")
cat("[INFO] file: ", opt$output, "\n")


cat("[INFO] Reading gene count matrix...\n")
if (file.info(opt$input)$size == 0) {
  stop("gene count matrix is empty: ", opt$input,
       "\nThe quantification step did not produce counts; fix it instead of continuing.",
       call. = FALSE)
}
database_all <- read.table(file = opt$input, sep = ",", header = TRUE)
cat("[INFO] Gene count matrix content:\n")
print(head(database_all))


if (!is.null(opt$metadata)) {
  meta <- read.csv(opt$metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("sample_name", "group") %in% colnames(meta))) {
    stop("metadata must have sample_name and group columns: ", opt$metadata, call. = FALSE)
  }
  # Replicates of one group share a sample_name, so the run ID is the only
  # column that can key a DESeq2 colData; fall back to sample_name for the
  # older two-column metadata that names every library once.
  key_col <- if ("sample" %in% colnames(meta)) "sample" else "sample_name"
  keys <- as.character(meta[[key_col]])
  if (anyDuplicated(keys)) {
    stop("metadata column ", key_col, " is not unique: ", opt$metadata, call. = FALSE)
  }
  config <- setNames(as.list(meta$group), keys)
} else {
  config <- yaml.load(opt$config)
}
cat("[INFO] Config content:\n")
print(config)


cat("[INFO] Processing expression matrix...\n")
expdata <- database_all[, -1]
rownames(expdata) <- database_all[, 1]
cat("[INFO] Expression matrix dimension: ", dim(expdata), "\n")


cat("[INFO] Extracting sample names and conditions...\n")
sample_names <- names(config)
groups <- sort(unique(unlist(config)))
known_groups <- c("control", "treatA", "treatB", "treatC", "treatD")
unknown <- setdiff(groups, known_groups)
if (length(unknown) > 0) {
  stop("unknown RNA-seq group(s): ", paste(unknown, collapse = ", "),
       "\nexpected control or treatA/treatB/treatC/treatD", call. = FALSE)
}
if ("control" %in% groups) {
  reference <- "control"
} else if (length(groups) == 2) {
  # Two-group designs written before control/treatA-D keep the alphabetical
  # reference they always had.
  reference <- groups[1]
} else {
  stop("RNA-seq groups must include 'control'; found: ",
       paste(groups, collapse = ", "), call. = FALSE)
}
treats <- setdiff(groups, reference)
if (length(treats) == 0) {
  stop("no treat group to compare against ", reference, call. = FALSE)
}
sizes <- table(unlist(config))
small <- names(sizes)[sizes < 2]
if (length(small) > 0) {
  stop("group(s) with fewer than 2 replicates: ", paste(small, collapse = ", "),
       call. = FALSE)
}
conditions <- factor(unlist(config), levels = c(reference, treats))
cat("[INFO] Reference: ", reference, "\n")
cat("[INFO] Contrasts: ", paste(treats, "vs", reference, collapse = ", "), "\n")
cat("[INFO] Sample:\n")
print(sample_names)
cat("[INFO] Conditions:\n")
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

cat("[INFO] Removing rows with NA values...\n")
expdata <- expdata[complete.cases(expdata), ]
cat("[INFO] Creating DESeq dataset...\n")
dds <- DESeqDataSetFromMatrix(countData = expdata, colData = coldata, design = ~ conditions)


cat("[INFO] Running DESeq2 analysis...\n")
dds <- DESeq(dds)


cat("[INFO] Extracting contrasts...\n")
normalised <- counts(dds, normalized = TRUE)
dir_out <- dirname(opt$output)
primary <- NULL
for (treat in treats) {
  contrast_name <- paste0(treat, "_vs_", reference)
  cat("[INFO] Contrast: ", contrast_name, "\n")
  res <- results(dds, contrast = c("conditions", treat, reference))
  members <- sample_names[as.character(conditions) %in% c(reference, treat)]
  # Merge on an explicit gene_id column: merge(by = "Row.names") is not the
  # spelling R accepts (it wants lowercase), so every contrast died here.
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  norm_df <- as.data.frame(normalised[, members, drop = FALSE])
  norm_df$gene_id <- rownames(norm_df)
  resdata <- merge(res_df, norm_df, by = "gene_id", sort = FALSE)

  resdata$SYMBOL <- resdata$gene_id
  resdata$gene_id <- sub("\\|[^|]*$", "", resdata$SYMBOL)
  resdata$SYMBOL <- sub("^.*\\|", "", resdata$SYMBOL)
  resdata <- resdata[!is.na(resdata$padj), ]

  contrast_path <- file.path(dir_out, paste0("DEG_result_", contrast_name, ".csv"))
  write.csv(resdata, contrast_path, row.names = FALSE)
  cat("[INFO] ", contrast_name, ": ", nrow(resdata), " genes -> ", contrast_path, "\n")
  if (is.null(primary)) primary <- resdata
}

# DEG_result.csv stays the first contrast so the joint analysis and the results
# page keep one well-known table; the other contrasts sit beside it.
write.csv(primary, opt$output, row.names = FALSE)

cat("[INFO] analysis: ", opt$output, "\n")


