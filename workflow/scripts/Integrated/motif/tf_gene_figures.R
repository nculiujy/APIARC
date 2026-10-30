#!/usr/bin/env Rscript

# TF-gene network drawn from the per-region fimo edges.
#
# One page per mark x comparison x region type x category, because the two peak
# directions make different claims and would otherwise be averaged on one page.
# The gene cap is 24: past that the right-hand labels collide at legible size,
# and the table beside the figure holds the full edge list.
# Edge width follows the best fimo p-value of the TF on that gene, not the motif
# q-value: q belongs to the matrix and is identical for every edge a motif
# contributes, so it cannot separate them.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--edges", type = "character"),
  make_option("--summary", type = "character", default = ""),
  make_option("--min-hit-regions", type = "integer", default = 2L, dest = "min_hit_regions"),
  make_option("--max-tfs", type = "integer", default = 12L, dest = "max_tfs"),
  make_option("--max-genes", type = "integer", default = 24L, dest = "max_genes"),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$edges) || is.null(opt$outdir)) {
  stop("missing --edges or --outdir")
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(sub("^--file=", "", script_arg[1]))
} else {
  "."
}
source(file.path(script_dir, "..", "common", "network_lib.R"))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
figure <- file.path(opt$outdir, "tf_gene_network.pdf")

edges <- if (file.exists(opt$edges)) read_tsv_plain(opt$edges) else NULL

placeholder <- function(message) {
  grDevices::pdf(figure, width = 7.2, height = 7.6)
  graphics::plot.new()
  graphics::text(0.5, 0.5, message, cex = 0.9)
  grDevices::dev.off()
  cat("[tf-gene] ", message, "; wrote a placeholder page\n", sep = "")
}

if (is.null(edges) || nrow(edges) == 0) {
  placeholder("no TF-gene edge passed the per-region hit filter")
  quit(save = "no")
}

edges$n_hit_regions <- suppressWarnings(as.integer(edges$n_hit_regions))
edges$best_p <- suppressWarnings(as.numeric(edges$best_p))
edges <- edges[!is.na(edges$n_hit_regions) & edges$n_hit_regions >= opt$min_hit_regions, ,
               drop = FALSE]
if (nrow(edges) == 0) {
  placeholder(sprintf("no TF-gene edge reached %d hit region(s)",
                      opt$min_hit_regions))
  quit(save = "no")
}

# Category is the peak direction, spelled so the page title reads as a claim.
category_label <- function(value) {
  switch(value,
         diff_up = "peaks up",
         diff_down = "peaks down",
         value)
}

short <- function(value, mark) {
  parts <- strsplit(value, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(value)
  strip <- function(x) sub(paste0("^", mark, "_"), "", x)
  sprintf("%s vs %s", strip(parts[1]), strip(parts[2]))
}

keys <- unique(edges[, c("mark", "comparison", "region_type", "category")])
keys <- keys[order(keys$mark, keys$comparison, keys$region_type, keys$category), ,
             drop = FALSE]

drawn <- 0L
grDevices::pdf(figure, width = 7.2, height = 7.6, onefile = TRUE)
for (index in seq_len(nrow(keys))) {
  key <- keys[index, ]
  page <- edges[
    edges$mark == key$mark & edges$comparison == key$comparison &
      edges$region_type == key$region_type & edges$category == key$category, ,
    drop = FALSE
  ]
  if (nrow(page) == 0) next

  # Degree on each side comes from the filtered page, so a TF kept on the page
  # has at least one edge to a gene that is also kept.
  tf_degree <- table(page$tf)
  gene_degree <- table(page$gene)
  tfs <- sort(tf_degree, decreasing = TRUE)
  genes <- sort(gene_degree, decreasing = TRUE)
  tfs <- tfs[seq_len(min(length(tfs), opt$max_tfs))]
  genes <- genes[seq_len(min(length(genes), opt$max_genes))]

  page <- page[page$tf %in% names(tfs) & page$gene %in% names(genes), , drop = FALSE]
  if (nrow(page) == 0) next

  left <- node_positions(tfs, 0.97, 0.03)
  right <- node_positions(genes, 0.97, 0.03)
  edge_list <- data.frame(
    left = page$tf, right = page$gene,
    weight = -log10(pmax(page$best_p, 1e-320)),
    stringsAsFactors = FALSE
  )
  draw_bipartite(
    left, right, edge_list,
    title = sprintf("%s  %s  %s  (%s)", key$mark, short(key$comparison, key$mark),
                    key$region_type, category_label(key$category)),
    subtitle = sprintf("%d TFs x %d genes drawn, >= %d hit region(s), edge width = -log10(best p)",
                       length(tfs), length(genes), opt$min_hit_regions),
    max_left_labels = opt$max_tfs,
    max_right_labels = opt$max_genes,
    cex_left = 0.62,
    cex_right = 0.6,
    # Edges must stay readable as individual claims: at the previous strength a
    # 12 x 24 page still merged into a grey band, and the width mapping carried
    # almost none of the contrast.
    # Room for the labels on both sides: gene symbols such as C530043K16Rik run
    # past the default right margin and are cut by the page edge.
    xlim = c(-0.30, 1.30),
    edge_alpha = c(0.12, 0.5),
    edge_lwd = c(0.3, 1.3)
  )
  drawn <- drawn + 1L
}
grDevices::dev.off()

if (drawn == 0L) {
  placeholder("every TF-gene page was empty after the filters")
  quit(save = "no")
}

if (nzchar(opt$summary) && file.exists(opt$summary)) {
  summary <- read_tsv_plain(opt$summary)
  cat(sprintf("[tf-gene] %d page(s); summary covers %d set(s), %d edges\n",
              drawn, nrow(summary),
              if ("n_edges" %in% names(summary)) sum(summary$n_edges) else 0L))
} else {
  cat(sprintf("[tf-gene] %d page(s) -> %s\n", drawn, figure))
}
