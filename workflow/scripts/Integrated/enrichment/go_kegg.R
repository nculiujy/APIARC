#!/usr/bin/env Rscript

# GO/KEGG over-representation for every region set of 0_stats.
#
# Two things are fixed relative to the previous version, both in enrich_lib.R:
# the gene-list size is the number of mapped symbols clusterProfiler used (not
# the number submitted), and the background is the same for every set (it used
# to be unioned with the set's own genes, which made each set test against a
# different universe). Rows whose overload still cannot be represented are
# dropped rather than written.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--region-sets", type = "character", help = "region_sets.tsv from 0_stats"),
  # Repeated --rna values are read with flag_values() in enrich_lib.R: this
  # optparse version keeps only the last one through `action = "list"`.
  make_option("--rna", type = "character", default = "",
              help = "DESeq2 tables, one per contrast; their union is the universe"),
  make_option("--universe", type = "character", default = "expressed",
              help = "expressed | genome | peaks"),
  make_option("--species", type = "character", default = "mm"),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--min-gene-set", type = "integer", default = 10L),
  make_option("--min-base-mean", type = "double", default = 10),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$`region-sets`) || is.null(opt$outdir)) {
  stop("missing --region-sets or --outdir")
}

required <- c("clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db")
installed <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!installed[["clusterProfiler"]]) {
  stop("clusterProfiler is missing; it must come from the pinned environment")
}

orgdb_name <- if (opt$species == "homo") "org.Hs.eg.db" else "org.Mm.eg.db"
kegg_organism <- if (opt$species == "homo") "hsa" else "mmu"
if (!requireNamespace(orgdb_name, quietly = TRUE)) {
  stop("missing annotation package ", orgdb_name)
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(orgdb_name, character.only = TRUE)
})
orgdb <- get(orgdb_name)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(sub("^--file=", "", script_arg[1]))
} else {
  "."
}
source(file.path(script_dir, "enrich_lib.R"))

rna_paths <- flag_values("--rna")
if (length(rna_paths) == 0) stop("missing --rna")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
}

symbols <- universe_from_rna(rna_paths, opt$`min-base-mean`)
universe <- switch(
  opt$universe,
  expressed = symbols$expressed,
  genome = symbols$all,
  peaks = NULL,
  stop("unknown --universe: ", opt$universe)
)

sets <- read_tsv(opt$`region-sets`)
if (nrow(sets) == 0) stop("region_sets.tsv is empty")

genes_of <- function(bed_path) {
  if (!file.exists(bed_path)) return(character())
  lines <- readLines(bed_path, warn = FALSE)
  values <- vapply(strsplit(lines, "\t", fixed = TRUE), function(f) {
    if (length(f) >= 4) f[4] else ""
  }, character(1))
  unique(values[nzchar(values)])
}

all_rows <- list()
for (i in seq_len(nrow(sets))) {
  set_row <- sets[i, ]
  genes <- genes_of(set_row$bed)
  if (length(genes) < opt$`min-gene-set`) next
  set_universe <- if (is.null(universe)) genes else universe
  frame <- enrich_one(genes, set_universe, set_row$set_id, orgdb, kegg_organism)
  if (is.null(frame)) next
  set_dir <- file.path(opt$outdir, set_row$set_id)
  dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(frame, file.path(set_dir, "enrichment.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
  all_rows[[length(all_rows) + 1]] <- frame
}

if (length(all_rows) == 0) {
  stop("no enrichment result was produced; check the region sets and the universe")
}
summary_table <- do.call(rbind, all_rows)
write.table(summary_table, file.path(opt$outdir, "enrichment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

contrast_rows <- list()
keys <- unique(sets[, c("mark", "comparison", "region_type")])
for (i in seq_len(nrow(keys))) {
  key <- keys[i, ]
  selected <- summary_table[
    summary_table$set %in% sets$set_id[
      sets$mark == key$mark & sets$comparison == key$comparison &
        sets$region_type == key$region_type
    ],
  ]
  if (nrow(selected) == 0) next
  selected$group <- ifelse(
    grepl("diff|bound", selected$set), "differential", "static"
  )
  differential <- selected[selected$group == "differential", ]
  static <- selected[selected$group == "static", ]
  if (nrow(static) == 0) next
  joined <- merge(
    differential[, c("ontology", "term_id", "term", "odds_ratio", "FDR", "n_genes", "overlap")],
    static[, c("term_id", "odds_ratio", "FDR", "overlap")],
    by = "term_id", suffixes = c("_diff", "_static")
  )
  if (nrow(joined) == 0) next
  joined$mark <- key$mark
  joined$comparison <- key$comparison
  joined$region_type <- key$region_type
  contrast_rows[[length(contrast_rows) + 1]] <- joined
}
if (length(contrast_rows) > 0) {
  contrast <- do.call(rbind, contrast_rows)
  write.table(
    contrast, file.path(opt$outdir, "differential_vs_static.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

# Impossible rows would mean the odds ratio is still built on the wrong K.
n_impossible <- sum(
  !is.na(summary_table$overlap) & summary_table$overlap > summary_table$n_term,
  na.rm = TRUE
)
n_negative_or <- sum(!is.na(summary_table$odds_ratio) & summary_table$odds_ratio < 0, na.rm = TRUE)
n_corrected <- sum(summary_table$odds_corrected, na.rm = TRUE)
n_no_overlap <- sum(!is.na(summary_table$overlap) & summary_table$overlap == 0, na.rm = TRUE)

writeLines(
  c(
    paste("universe:", opt$universe),
    paste("universe_size:", length(universe)),
    paste("universe_min_base_mean:", opt$`min-base-mean`),
    paste("rna_tables:", paste(basename(rna_paths), collapse = ",")),
    paste("region_sets:", nrow(sets)),
    paste("sets_enriched:", length(all_rows)),
    paste("alpha:", opt$alpha),
    paste("rows_overlap_gt_term:", n_impossible),
    paste("rows_negative_odds_ratio:", n_negative_or),
    paste("rows_haldane_corrected:", n_corrected),
    paste("rows_no_overlap:", n_no_overlap)
  ),
  file.path(opt$outdir, "universe.txt")
)

cat(sprintf(
  "[enrichment] %d sets, %d rows, universe=%s (%d genes), impossible rows=%d\n",
  length(all_rows), nrow(summary_table), opt$universe, length(universe), n_impossible
))
