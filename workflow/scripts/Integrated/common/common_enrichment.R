#!/usr/bin/env Rscript

# GO/KEGG for the common genes, with the differential-only list beside it.
#
# The enrichment stage tests region sets, i.e. every gene that lost or gained
# binding. The question here is narrower: what does the intersection with
# differential expression add? Both gene lists are tested against the same
# universe, so the two tables differ only in the gene list, and a term that
# gains strength in the common list is a term whose genes are also regulated.
#
# The gene-to-term edges for the network come from the common list only; the
# all-differential list is a control, not a second network.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--common", type = "character", help = "common_genes.tsv"),
  make_option("--gene-table", type = "character", help = "chip_gene_table.tsv"),
  make_option("--rna", type = "character", default = "",
              help = "DESeq2 tables; their union is the universe"),
  make_option("--species", type = "character", default = "mm"),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--min-gene-set", type = "integer", default = 10L),
  make_option("--max-terms", type = "integer", default = 12L),
  make_option("--min-base-mean", type = "double", default = 10),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$common) || is.null(opt$`gene-table`) || is.null(opt$outdir)) {
  stop("missing --common, --gene-table or --outdir")
}

orgdb_name <- if (opt$species == "homo") "org.Hs.eg.db" else "org.Mm.eg.db"
kegg_organism <- if (opt$species == "homo") "hsa" else "mmu"
if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
    !requireNamespace(orgdb_name, quietly = TRUE)) {
  stop("clusterProfiler or ", orgdb_name, " is missing")
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
source(file.path(script_dir, "..", "enrichment", "enrich_lib.R"))

rna_paths <- flag_values("--rna")
if (length(rna_paths) == 0) {
  stop("missing --rna: the universe has to be declared")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
}

universe <- universe_from_rna(rna_paths, opt$`min-base-mean`)$expressed
if (length(universe) < 100) {
  stop("the universe has ", length(universe), " genes; check --rna")
}

common <- read_tsv(opt$common)
if (nrow(common) == 0) stop("common_genes.tsv is empty")
gene_table <- read_tsv(opt$`gene-table`)

key_of <- function(frame) {
  unique(frame[, c("mark", "comparison", "region_type")])
}

# Distinct gene lists per key. The common list keeps the concordant rows only;
# the control list is every gene of the same key that carries a differential
# category, which is the list the region-set stage already enriched.
lists <- list()
for (key_index in seq_len(nrow(key_of(common)))) {
  key <- key_of(common)[key_index, ]
  same <- function(frame) {
    frame$mark == key$mark & frame$comparison == key$comparison &
      frame$region_type == key$region_type
  }
  common_rows <- common[same(common) & common$concordant == "true", , drop = FALSE]
  control_rows <- gene_table[
    same(gene_table) & gene_table$category %in% c("diff_up", "diff_down") &
      gene_table$expressed == "true", , drop = FALSE
  ]
  lists[[key_index]] <- list(
    key = key,
    common = sort(unique(common_rows$gene)),
    all_diff = sort(unique(control_rows$gene))
  )
}

run_lists <- function(which) {
  rows <- list()
  for (entry in lists) {
    genes <- entry[[which]]
    if (length(genes) < opt$`min-gene-set`) next
    label <- sprintf("%s__%s__%s__%s", entry$key$mark, entry$key$comparison,
                     entry$key$region_type, which)
    frame <- enrich_one(genes, universe, label, orgdb, kegg_organism)
    if (is.null(frame)) next
    frame$mark <- entry$key$mark
    frame$comparison <- entry$key$comparison
    frame$region_type <- entry$key$region_type
    frame$list_name <- which
    rows[[length(rows) + 1]] <- frame
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

common_table <- run_lists("common")
control_table <- run_lists("all_diff")
if (is.null(common_table) && is.null(control_table)) {
  stop("no gene list reached --min-gene-set; nothing to enrich")
}

columns <- c("set", "mark", "comparison", "region_type", "list_name", "ontology",
             "term_id", "term", "n_submitted", "n_query", "n_genes", "n_term",
             "overlap", "odds_ratio", "odds_corrected", "p", "FDR", "genes")
write_one <- function(frame, name) {
  if (is.null(frame)) return(invisible(NULL))
  write.table(frame[, columns], file.path(opt$outdir, name), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}
write_one(common_table, "go_kegg_common.tsv")
write_one(control_table, "go_kegg_all_diff.tsv")

# Side-by-side view: same term, two gene lists. Terms absent from one list keep
# NA rather than being dropped, so "the term disappeared" is visible.
if (!is.null(common_table) && !is.null(control_table)) {
  key_columns <- c("mark", "comparison", "region_type", "ontology", "term_id", "term")
  compare <- merge(
    common_table[, c(key_columns, "overlap", "n_query", "odds_ratio", "FDR")],
    control_table[, c(key_columns, "overlap", "n_query", "odds_ratio", "FDR")],
    by = key_columns, all = TRUE, suffixes = c("_common", "_all_diff")
  )
  compare <- compare[order(compare$mark, compare$comparison, compare$region_type,
                           compare$FDR_common, compare$FDR_all_diff), , drop = FALSE]
  write.table(compare, file.path(opt$outdir, "term_compare.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}

# Network edges: significant terms of the common list only.
edges <- term_edges(common_table, opt$alpha, opt$`max-terms`)
if (!is.null(edges)) {
  write.table(edges, file.path(opt$outdir, "gene_term_edges.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
} else {
  write.table(
    data.frame(mark = character(), comparison = character(), region_type = character(),
               ontology = character(), term_id = character(), term = character(),
               gene = character(), FDR = numeric(), odds_ratio = numeric()),
    file.path(opt$outdir, "gene_term_edges.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA"
  )
}

# Table form of the network the figure draws: one file for the edges with
# weight = -log10(FDR), one for the nodes with the degree of the written edges,
# so every line and every node can be traced back to a row here.
network_columns <- c("gene", "term_id", "ontology", "mark", "comparison",
                     "region_type", "weight", "FDR")
node_columns <- c("mark", "comparison", "region_type", "id", "label", "type",
                  "ontology", "fdr", "degree")
if (!is.null(edges) && nrow(edges) > 0) {
  network <- data.frame(
    gene = edges$gene, term_id = edges$term_id, ontology = edges$ontology,
    mark = edges$mark, comparison = edges$comparison,
    region_type = edges$region_type,
    weight = -log10(pmax(edges$FDR, 1e-320)), FDR = edges$FDR,
    stringsAsFactors = FALSE
  )
  write.table(network[, network_columns], file.path(opt$outdir, "network_edges.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  term_degree <- aggregate(
    list(degree = network$gene),
    by = network[, c("mark", "comparison", "region_type", "term_id")],
    FUN = function(genes) length(unique(genes))
  )
  term_label <- unique(network[, c("mark", "comparison", "region_type", "term_id", "FDR")])
  term_label$term <- edges$term[match(
    paste(term_label$mark, term_label$comparison, term_label$region_type, term_label$term_id),
    paste(edges$mark, edges$comparison, edges$region_type, edges$term_id)
  )]
  term_label$ontology <- edges$ontology[match(
    paste(term_label$mark, term_label$comparison, term_label$region_type, term_label$term_id),
    paste(edges$mark, edges$comparison, edges$region_type, edges$term_id)
  )]
  term_nodes <- merge(term_label, term_degree,
                      by = c("mark", "comparison", "region_type", "term_id"))
  term_nodes <- data.frame(
    mark = term_nodes$mark, comparison = term_nodes$comparison,
    region_type = term_nodes$region_type, id = term_nodes$term_id,
    label = term_nodes$term, type = "term", ontology = term_nodes$ontology,
    fdr = term_nodes$FDR, degree = term_nodes$degree, stringsAsFactors = FALSE
  )
  gene_degree <- aggregate(
    list(degree = network$term_id),
    by = network[, c("mark", "comparison", "region_type", "gene")],
    FUN = function(terms) length(unique(terms))
  )
  gene_nodes <- data.frame(
    mark = gene_degree$mark, comparison = gene_degree$comparison,
    region_type = gene_degree$region_type, id = gene_degree$gene,
    label = gene_degree$gene, type = "gene", ontology = NA_character_,
    fdr = NA_real_, degree = gene_degree$degree, stringsAsFactors = FALSE
  )
  nodes <- rbind(term_nodes[, node_columns], gene_nodes[, node_columns])
  nodes <- nodes[order(nodes$mark, nodes$comparison, nodes$region_type, nodes$type,
                       -nodes$degree), , drop = FALSE]
  write.table(nodes, file.path(opt$outdir, "network_nodes.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
} else {
  write.table(
    data.frame(gene = character(), term_id = character(), ontology = character(),
               mark = character(), comparison = character(), region_type = character(),
               weight = numeric(), FDR = numeric()),
    file.path(opt$outdir, "network_edges.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.table(
    data.frame(mark = character(), comparison = character(), region_type = character(),
               id = character(), label = character(), type = character(),
               ontology = character(), fdr = numeric(), degree = integer()),
    file.path(opt$outdir, "network_nodes.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA"
  )
}

# Per key and ontology counts, so the figure caption can state how many terms
# were tested rather than only how many were drawn.
summary_rows <- list()
for (table in list(common_table, control_table)) {
  if (is.null(table)) next
  keys <- unique(table[, c("mark", "comparison", "region_type", "list_name", "ontology")])
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    selected <- table[
      table$mark == key$mark & table$comparison == key$comparison &
        table$region_type == key$region_type & table$list_name == key$list_name &
        table$ontology == key$ontology,
    ]
    if (nrow(selected) == 0) next
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      mark = key$mark, comparison = key$comparison, region_type = key$region_type,
      list_name = key$list_name, ontology = key$ontology,
      n_genes_query = selected$n_query[1], n_genes_mapped = selected$n_genes[1],
      n_terms_tested = nrow(selected),
      n_terms_fdr = sum(!is.na(selected$FDR) & selected$FDR < opt$alpha),
      min_fdr = suppressWarnings(min(selected$FDR, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
}
if (length(summary_rows) > 0) {
  summary_frame <- do.call(rbind, summary_rows)
  summary_frame$min_fdr[!is.finite(summary_frame$min_fdr)] <- NA
  write.table(summary_frame, file.path(opt$outdir, "gene_term_summary.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}

n_common <- if (is.null(common_table)) 0 else nrow(common_table)
n_control <- if (is.null(control_table)) 0 else nrow(control_table)
n_edges <- if (is.null(edges)) 0 else nrow(edges)
cat(sprintf(
  "[common-enrichment] %d keys, %d common rows, %d all-diff rows, %d term-gene edges\n",
  length(lists), n_common, n_control, n_edges
))
