#!/usr/bin/env Rscript

# Gene-term network of the common genes, one page per mark x comparison x region
# type. The page limit is not decoration: 12 terms against 25 genes is what fits
# on a page and still carries legible labels, and the table next to the figure
# holds the rest.
#
# A second figure, enrichment_bar.pdf, ranks the same keys by enrichment
# strength: the bipartite page shows which genes carry which term but cannot
# compare terms to one another, and the KEGG rows drown in the BP half there.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--edges", type = "character"),
  make_option("--nodes", type = "character"),
  make_option("--enrich", type = "character", default = ""),
  make_option("--alpha", type = "numeric", default = 0.05),
  make_option("--max_terms", type = "integer", default = 12L),
  make_option("--max_genes", type = "integer", default = 40L),
  # Hyphenated options need an explicit dest: optparse would otherwise store
  # them under a dotted name and the lookup below would take NULL.
  make_option("--max-bar-terms", type = "integer", default = 12L,
              dest = "max_bar_terms"),
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
source(file.path(script_dir, "network_lib.R"))

edges <- read_tsv_plain(opt$edges)
nodes <- if (!is.null(opt$nodes) && file.exists(opt$nodes)) read_tsv_plain(opt$nodes) else NULL
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
figure <- file.path(opt$outdir, "network.pdf")

# Comparison strings carry the mark twice ("H3K27ac|H3K27ac_TGFb12h_vs_...");
# the title reads better with the prefix stripped once.
short <- function(value, mark) {
  parts <- strsplit(value, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(value)
  strip <- function(x) sub(paste0("^", mark, "_"), "", x)
  sprintf("%s vs %s", strip(parts[1]), strip(parts[2]))
}

# ---------------------------------------------------------------------------
# Enrichment strength, one page per network key. All significant terms rank in
# a single horizontal bar plot (GO-BP and KEGG share the axis), bar = -log10
# (FDR), fill = term type with a legend on the right, light gridlines and the
# overlap count at the bar end. Drawn before the network section so the
# empty-network early exit still leaves the file behind.
# ---------------------------------------------------------------------------
enrich_figure <- file.path(opt$outdir, "enrichment_bar.pdf")

if (nzchar(opt$enrich) && file.exists(opt$enrich)) {
  enrich <- read_tsv_plain(opt$enrich)
  enrich$FDR <- suppressWarnings(as.numeric(enrich$FDR))
  enrich$overlap <- suppressWarnings(as.integer(enrich$overlap))
  enrich <- enrich[!is.na(enrich$FDR) & enrich$FDR < opt$alpha, , drop = FALSE]

  # Okabe-Ito palette: distinguishable under red-green colour blindness, and
  # ontologies absent from the table never reach the legend.
  ontology_colours <- c(BP = "#E69F00", MF = "#56B4E9", CC = "#009E73",
                        KEGG = "#56B4E9")
  ontology_titles <- c(BP = "Biological process", MF = "Molecular function",
                       CC = "Cellular component", KEGG = "KEGG")

  draw_enrich_bars <- function(part, page_title, n_rows) {
    n <- nrow(part)
    heights <- -log10(pmax(part$FDR, 1e-300))
    fills <- ontology_colours[part$ontology]
    x_max <- max(heights) * 1.1
    # Right margin carries the "Term type" legend outside the panel.
    graphics::par(mar = c(4.4, 16.5, 3.2, 10.0))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, x_max), ylim = c(0.5, n_rows + 0.5))
    # Full grid behind the bars: vertical lines at the axis ticks, horizontal
    # lines at every row boundary, like the reference layout. The row count is
    # fixed per page so a page with few terms keeps the same bar thickness as
    # a full one; the bars hang from the top.
    graphics::abline(v = graphics::axTicks(1), col = "grey90")
    graphics::abline(h = seq(0.5, n_rows + 0.5, by = 1), col = "grey90")
    # Bars occupy 70% of each row so adjacent terms stay visually separate.
    mid <- n_rows - n + seq_len(n)
    graphics::rect(0, mid - 0.35, heights, mid + 0.35, col = fills, border = NA)
    graphics::axis(1, cex.axis = 0.8)
    graphics::box(col = "grey40")
    graphics::text(-x_max * 0.03, mid, labels = trim_text(part$term, 48),
                   pos = 2, xpd = NA, cex = 0.62)
    graphics::text(heights, mid, labels = as.character(part$overlap), pos = 4,
                   cex = 0.55, col = "grey30")
    graphics::title(xlab = expression(-log[10](FDR)), line = 2.9, cex.lab = 0.9)
    graphics::mtext(page_title, side = 3, line = 1.4, adj = 0, cex = 0.9,
                    font = 2, xpd = NA)
    shown <- unique(part$ontology)
    graphics::legend(x = x_max * 1.03, y = n_rows * 0.6, xpd = NA, bty = "n",
                     legend = ontology_titles[shown],
                     fill = ontology_colours[shown], title = "Term type",
                     cex = 0.85)
  }

  grDevices::pdf(enrich_figure, width = 9.6, height = 8.4, onefile = TRUE)
  if (nrow(enrich) == 0) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, sprintf("no enriched term below FDR %g", opt$alpha), cex = 0.9)
  } else {
    bar_keys <- unique(enrich[, c("mark", "comparison", "region_type")])
    bar_keys <- bar_keys[order(bar_keys$mark, bar_keys$comparison,
                               bar_keys$region_type), , drop = FALSE]
    for (index in seq_len(nrow(bar_keys))) {
      key <- bar_keys[index, ]
      page <- enrich[enrich$mark == key$mark &
                       enrich$comparison == key$comparison &
                       enrich$region_type == key$region_type, , drop = FALSE]
      part <- page[order(page$FDR), , drop = FALSE]
      part <- part[seq_len(min(nrow(part), opt$max_bar_terms)), , drop = FALSE]
      if (nrow(part) == 0) next
      # barplot draws the first group at the bottom: reverse so the strongest
      # term ends up on top.
      part <- part[nrow(part):1, , drop = FALSE]
      draw_enrich_bars(
        part,
        sprintf("%s  %s  (%s, top %d of %d terms, genes tested: %d)",
                key$mark, short(key$comparison, key$mark), key$region_type,
                nrow(part), nrow(page),
                max(suppressWarnings(as.numeric(page$n_query)), na.rm = TRUE)),
        n_rows = opt$max_bar_terms)
    }
  }
  grDevices::dev.off()
  cat(sprintf("[enrichment_bar] -> %s\n", enrich_figure))
}

if (nrow(edges) == 0 || is.null(nodes) || nrow(nodes) == 0) {
  grDevices::pdf(figure, width = 7.2, height = 7.6)
  graphics::plot.new()
  graphics::text(0.5, 0.5, "no significant term shared by the common genes", cex = 0.9)
  grDevices::dev.off()
  cat("[network] empty network; wrote a placeholder page\n")
  quit(save = "no")
}

keys <- unique(edges[, c("mark", "comparison", "region_type")])

grDevices::pdf(figure, width = 7.2, height = 7.6, onefile = TRUE)
for (index in seq_len(nrow(keys))) {
  key <- keys[index, ]
  page_edges <- edges[
    edges$mark == key$mark & edges$comparison == key$comparison &
      edges$region_type == key$region_type, , drop = FALSE
  ]
  page_nodes <- nodes[
    nodes$mark == key$mark & nodes$comparison == key$comparison &
      nodes$region_type == key$region_type, , drop = FALSE
  ]
  if (nrow(page_edges) == 0) next

  # The per-ontology budget matches term_edges: BP and KEGG are separate claims
  # and a shared cap would let the BP half fill the page again.
  terms <- page_nodes[page_nodes$type == "term", , drop = FALSE]
  terms <- do.call(rbind, lapply(split(terms, terms$ontology), function(part) {
    part <- part[order(-part$degree, part$fdr), , drop = FALSE]
    head(part, opt$max_terms)
  }))
  terms <- terms[order(-terms$degree, terms$fdr), , drop = FALSE]
  genes <- page_nodes[page_nodes$type == "gene", , drop = FALSE]
  genes <- genes[order(-genes$degree), , drop = FALSE]
  genes <- head(genes, opt$max_genes)

  page_edges <- page_edges[page_edges$term_id %in% terms$id & page_edges$gene %in% genes$id, , drop = FALSE]
  if (nrow(page_edges) == 0) next

  term_degrees <- stats::setNames(terms$degree, terms$id)
  gene_degrees <- stats::setNames(genes$degree, genes$id)
  left <- node_positions(term_degrees, 0.97, 0.03)
  right <- node_positions(gene_degrees, 0.97, 0.03)

  # Term labels carry the ontology, because "BP" and "KEGG" terms with the same
  # wording are different claims. The prefix goes to text, not labels: labels
  # stay the term ids the edges are looked up by.
  ontology <- stats::setNames(terms$ontology, terms$id)
  term_names <- stats::setNames(terms$label, terms$id)
  left$text <- trim_text(sprintf("[%s] %s", ontology[left$labels],
                                 term_names[left$labels]), 44)

  edge_list <- data.frame(
    left = page_edges$term_id, right = page_edges$gene, weight = page_edges$weight,
    stringsAsFactors = FALSE
  )
  draw_bipartite(
    left, right, edge_list,
    title = sprintf("%s  %s  (%s)", key$mark, short(key$comparison, key$mark), key$region_type),
    subtitle = sprintf("%d terms x %d genes drawn, edge width = -log10(FDR)",
                       nrow(terms), nrow(genes)),
    max_left_labels = opt$max_terms,
    max_right_labels = opt$max_genes,
    cex_left = 0.62,
    cex_right = 0.5,
    xlim = c(-0.5, 1.1)
  )
}
grDevices::dev.off()
cat(sprintf("[network] %d page(s) over %d key(s) -> %s\n", nrow(keys), nrow(keys), figure))
