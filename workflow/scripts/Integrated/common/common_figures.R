#!/usr/bin/env Rscript

# Figures for the common genes of the joint analysis.
#
# Everything is base graphics: the analysis environments carry r-base and
# optparse only, and 0_4_plots.R already draws its figures this way, so the new
# stage does not add a dependency for four plots.
#
# The correlation printed under each scatter is the genome-wide rho, not the
# one inside the common set. The common set is filtered on binding and on
# expression, so a correlation computed inside it is inflated by construction
# (0_3_tests.R keeps the same number under h2_correlation_selective and marks it
# "BIAS REFERENCE"). Both values live in common_stats.tsv; the figure cites the
# unbiased one.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--common", type = "character", help = "common_genes.tsv"),
  make_option("--expression", type = "character", help = "expression_matrix.tsv"),
  make_option("--signal", type = "character", help = "signal_matrix.tsv"),
  make_option("--stats", type = "character", help = "common_stats.tsv"),
  make_option("--sample-map", type = "character", help = "sample_map.tsv"),
  make_option("--rna", type = "character", action = "list", default = list(),
              help = "DESeq2 tables, for the all-gene expression background"),
  make_option("--heatmap-rows", type = "integer", default = 100L),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$common) || is.null(opt$expression) || is.null(opt$outdir)) {
  stop("missing --common, --expression or --outdir")
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
}

as_num <- function(x) suppressWarnings(as.numeric(x))

common <- read_tsv(opt$common)
expression <- read_tsv(opt$expression)
signal <- if (!is.null(opt$signal)) read_tsv(opt$signal) else NULL
stats <- if (!is.null(opt$stats)) read_tsv(opt$stats) else NULL

if (nrow(common) == 0) {
  stop("common_genes.tsv is empty; there is nothing to plot")
}

# --- sample layout ---------------------------------------------------------
sample_columns <- setdiff(colnames(expression), "gene")
# "TGFb12h_rep1" -> "TGFb12h"; the suffix after the last underscore is the
# replicate, so a time point that itself contains one still parses.
timepoint_of <- function(column) sub("_[^_]+$", "", column)
timepoints <- unique(vapply(sample_columns, timepoint_of, character(1)))
samples_of <- function(timepoint) {
  sample_columns[vapply(sample_columns, timepoint_of, character(1)) == timepoint]
}
signal_columns <- if (!is.null(signal)) setdiff(
  colnames(signal), c("gene", "mark", "comparison", "region_type", "n_regions",
                      "chip_direction", "concordant")
) else character(0)

# run id -> time point, so the DESeq2 tables can be summarised by time point
# instead of pooled across the two contrasts.
run_timepoint <- NULL
if (!is.null(opt$`sample-map`) && file.exists(opt$`sample-map`)) {
  mapping <- read_tsv(opt$`sample-map`)
  mapping <- mapping[mapping$omics == "RNAseq", , drop = FALSE]
  run_timepoint <- setNames(as.character(mapping$timepoint), as.character(mapping$sample))
}

expression_of <- function(gene) {
  row <- expression[expression$gene == gene, , drop = FALSE]
  if (nrow(row) == 0) return(rep(NA_real_, length(sample_columns)))
  vapply(sample_columns, function(column) as_num(row[[column]])[1], numeric(1))
}

# --- per-gene expression curves -------------------------------------------
# log2(x + 1) keeps the two- and three-fold changes readable on one axis and
# avoids dropping genes with a zero in one replicate.
gene_curves <- function(genes) {
  if (length(genes) == 0) return(NULL)
  curves <- t(vapply(genes, expression_of, numeric(length(sample_columns))))
  colnames(curves) <- sample_columns
  log2(curves + 1)
}

# Background: median expression of every gene in the DESeq2 tables, one value
# per time point. The common set is filtered on both axes, so without this
# reference a reader cannot tell whether the up group rises above the
# transcriptome or only above the down group.
background_by_timepoint <- function() {
  if (length(opt$rna) == 0 || is.null(run_timepoint)) return(NULL)
  values <- list()
  for (path in opt$rna) {
    frame <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    for (column in intersect(names(run_timepoint), colnames(frame))) {
      timepoint <- run_timepoint[[column]]
      x <- as_num(frame[[column]])
      values[[timepoint]] <- c(values[[timepoint]], x[!is.na(x)])
    }
  }
  if (length(values) == 0) return(NULL)
  vapply(values, function(x) median(log2(x + 1)), numeric(1))
}

combos <- unique(common[, c("mark", "comparison")])
combos <- combos[order(combos$mark, combos$comparison), , drop = FALSE]
region_levels <- c("promoter", "distal")

# "H3K27ac_TGFb12h_vs_H3K27ac_veh" -> "TGFb12h vs veh": the full name does not
# fit a panel title once the mark and the region are added.
short_comparison <- function(comparison, mark) {
  parts <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(comparison)
  drop_mark <- function(value) sub(paste0("^", mark, "_"), "", value)
  sprintf("%s vs %s", drop_mark(parts[1]), drop_mark(parts[2]))
}

panel_label <- function(mark, comparison) {
  sprintf("%s %s", mark, short_comparison(comparison, mark))
}

subset_of <- function(mark, comparison, region_type) {
  common[common$mark == mark & common$comparison == comparison &
           common$region_type == region_type, , drop = FALSE]
}

# --- figure 1: ChIP change against RNA change ------------------------------
scatter_panel <- function(frame, label) {
  if (nrow(frame) < 5) {
    plot.new()
    title(main = label, cex.main = 0.9)
    text(0.5, 0.5, sprintf("only %d genes", nrow(frame)), col = "grey40", cex = 0.9)
    return(invisible(NULL))
  }
  plot(
    frame$chip_lfc_mean, frame$rna_log2FC,
    pch = ifelse(frame$concordant == "true", 21, 1),
    bg = ifelse(frame$chip_direction == "up", "steelblue", "indianred"),
    col = ifelse(frame$concordant == "true",
                 ifelse(frame$chip_direction == "up", "steelblue", "indianred"), "grey50"),
    cex = 0.9, xlab = "ChIP log2 fold change (mean over regions)",
    ylab = "RNA log2 fold change", main = label, las = 1, bty = "n"
  )
  abline(h = 0, v = 0, lty = 2, col = "grey60")
  entry <- stats[
    stats$mark == frame$mark[1] & stats$comparison == frame$comparison[1] &
      stats$region_type == frame$region_type[1] & stats$set == "common_concordant",
    , drop = FALSE
  ]
  if (nrow(entry) == 1) {
    legend("topleft", bty = "n", cex = 0.75, pch = c(16, 16, 1),
           pt.bg = c("steelblue", "indianred", NA), col = c("steelblue", "indianred", "grey50"),
           legend = c(
             sprintf("concordant up (n=%d)", sum(frame$concordant == "true" & frame$chip_direction == "up")),
             sprintf("concordant down (n=%d)", sum(frame$concordant == "true" & frame$chip_direction == "down")),
             sprintf("discordant (n=%d)", sum(frame$concordant != "true"))
           ))
    mtext(
      sprintf("rho (all genes) = %.2f, n = %d;  rho (common) = %.2f, n = %d",
              as_num(entry$spearman_all_genes), as_num(entry$n_all_genes),
              as_num(entry$spearman), as_num(entry$n)),
      side = 1, line = 4.3, cex = 0.72, adj = 0
    )
    # Same sentence as common_stats.tsv$note and docs/STATISTICS.md section 9:
    # the common-set rho describes a set both variables were used to cut, so it
    # cannot stand in for a test of the relationship.
    mtext(
      "rho (common) is computed after both filters; cite rho (all genes)",
      side = 1, line = 5.3, cex = 0.68, adj = 0, col = "grey30"
    )
  }
}

pdf(file.path(opt$outdir, "scatter_change.pdf"), width = 10, height = 4.6)
# Extra bottom margin: the two rho values and the selection-bias line sit under
# the x-axis label.
par(mfrow = c(1, 2), mar = c(6.2, 4.4, 2.6, 1))
for (index in seq_len(nrow(combos))) {
  mark <- combos$mark[index]
  comparison <- combos$comparison[index]
  for (region_type in region_levels) {
    frame <- subset_of(mark, comparison, region_type)
    scatter_panel(frame, sprintf("%s (%s)", panel_label(mark, comparison), region_type))
  }
}
invisible(dev.off())

# --- figure 2: expression across the time course ---------------------------
background <- background_by_timepoint()

line_panel <- function(frame, label) {
  if (nrow(frame) == 0) {
    plot.new()
    title(main = label, cex.main = 0.9)
    text(0.5, 0.5, "no genes", col = "grey40")
    return(invisible(NULL))
  }
  groups <- list(
    up = frame$gene[frame$concordant == "true" & frame$chip_direction == "up"],
    down = frame$gene[frame$concordant == "true" & frame$chip_direction == "down"]
  )
  curves <- lapply(groups, gene_curves)
  limits <- range(
    unlist(lapply(curves, function(curve) if (is.null(curve)) NA else range(curve, na.rm = TRUE))),
    background, na.rm = TRUE
  )
  plot(NA, xlim = c(0.8, length(timepoints) + 0.2), ylim = limits,
       xaxt = "n", xlab = "", ylab = "log2(normalized counts + 1)",
       main = label, las = 1, bty = "n")
  axis(1, at = seq_along(timepoints), labels = timepoints)
  if (!is.null(background)) {
    for (index in seq_along(timepoints)) {
      value <- background[timepoints[index]]
      if (is.na(value)) next
      segments(index - 0.35, value, index + 0.35, value, col = "grey50", lty = 2)
    }
  }
  colours <- c(up = "steelblue", down = "indianred")
  for (name in names(curves)) {
    curve <- curves[[name]]
    if (is.null(curve) || nrow(curve) == 0) next
    positions <- seq_along(timepoints)
    # One value per gene per time point (mean over that time point's
    # replicates). The group is drawn as a median line with an interquartile
    # band: the earlier per-gene hairlines numbered in the hundreds and buried
    # the trend they were meant to reveal.
    gene_means <- vapply(timepoints, function(timepoint) {
      rowMeans(curve[, samples_of(timepoint), drop = FALSE], na.rm = TRUE)
    }, numeric(nrow(curve)))
    # A group holding exactly one gene makes vapply return a plain vector, and
    # apply() below then aborts with "dim(X) must have a positive length"; the
    # explicit nrow restores the gene x time point shape for that case.
    gene_means <- matrix(gene_means, nrow = nrow(curve))
    if (nrow(curve) > 1) {
      q25 <- apply(gene_means, 2, quantile, probs = 0.25, na.rm = TRUE)
      q75 <- apply(gene_means, 2, quantile, probs = 0.75, na.rm = TRUE)
      polygon(c(positions, rev(positions)), c(q25, rev(q75)),
              col = adjustcolor(colours[[name]], alpha.f = 0.15), border = NA)
    }
    median_curve <- apply(gene_means, 2, median, na.rm = TRUE)
    lines(positions, median_curve, col = colours[[name]], lwd = 2.4)
    points(positions, median_curve, pch = 21, bg = colours[[name]], cex = 1.3)
  }
  legend("topleft", bty = "n", cex = 0.72, lwd = c(2.4, 2.4, 1), lty = c(1, 1, 2),
         col = c("steelblue", "indianred", "grey50"),
         legend = c(sprintf("ChIP up (n=%d, median, IQR band)", length(groups$up)),
                    sprintf("ChIP down (n=%d, median, IQR band)", length(groups$down)),
                    "all genes (median)"))
}

pdf(file.path(opt$outdir, "line_expression.pdf"), width = 10, height = 4.6)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1))
for (index in seq_len(nrow(combos))) {
  mark <- combos$mark[index]
  comparison <- combos$comparison[index]
  for (region_type in region_levels) {
    line_panel(subset_of(mark, comparison, region_type),
               sprintf("%s (%s)", panel_label(mark, comparison), region_type))
  }
}
invisible(dev.off())

# --- figure 3: per-sample expression heatmap -------------------------------
heatmap_page <- function(frame, label) {
  # heatmap() calls layout() internally, which overrides mfrow; resetting it per
  # page keeps one heatmap on one page.
  par(mfrow = c(1, 1))
  frame <- frame[frame$concordant == "true", , drop = FALSE]
  if (nrow(frame) < 5) {
    plot.new()
    title(main = label)
    text(0.5, 0.5, "fewer than 5 concordant genes", col = "grey40")
    return(invisible(NULL))
  }
  n_total <- nrow(frame)
  frame <- frame[order(-abs(frame$chip_lfc_mean)), , drop = FALSE]
  shown <- min(n_total, opt$`heatmap-rows`)
  if (n_total > opt$`heatmap-rows`) {
    frame <- frame[seq_len(opt$`heatmap-rows`), , drop = FALSE]
  }
  matrix_values <- t(vapply(frame$gene, expression_of, numeric(length(sample_columns))))
  colnames(matrix_values) <- sample_columns
  rownames(matrix_values) <- frame$gene
  matrix_values <- log2(matrix_values + 1)
  centred <- t(apply(matrix_values, 1, function(row) {
    spread <- sd(row, na.rm = TRUE)
    if (is.na(spread) || spread == 0) return(rep(0, length(row)))
    (row - mean(row, na.rm = TRUE)) / spread
  }))
  centred[is.na(centred)] <- 0
  # Cluster genes (rows) so co-regulated blocks line up; samples (columns) keep
  # the fixed veh -> 2h -> 12h order, so the time course stays readable.
  row_dend <- NA
  if (nrow(centred) > 2) {
    row_dist <- dist(centred)
    if (any(is.finite(row_dist) & row_dist > 0)) {
      row_dend <- as.dendrogram(hclust(row_dist, method = "average"))
    }
  }
  palette <- colorRampPalette(c("#6baed6", "#f7f7f7", "#f4a582"))(51)
  side <- ifelse(frame$chip_direction == "up", "#4292c6", "#ef6c4d")
  heatmap(
    centred, Rowv = row_dend, Colv = NA, scale = "none", col = palette,
    margins = c(7, 9), cexRow = 0.45, cexCol = 0.8,
    RowSideColors = side, labRow = rownames(centred), labCol = colnames(centred),
    main = sprintf("%s: z-scored expression (top %d of %d concordant genes)",
                   label, shown, n_total)
  )
  legend("topright", inset = c(-0.02, 0), bty = "n", cex = 0.7, fill = c("#4292c6", "#ef6c4d"),
         legend = c("ChIP up", "ChIP down"), xpd = TRUE)
}

pdf(file.path(opt$outdir, "heatmap_expression.pdf"), width = 8, height = 9)
par(mar = c(4, 4, 3, 1))
for (index in seq_len(nrow(combos))) {
  mark <- combos$mark[index]
  comparison <- combos$comparison[index]
  for (region_type in region_levels) {
    heatmap_page(subset_of(mark, comparison, region_type),
                 sprintf("%s %s", panel_label(mark, comparison), region_type))
  }
}
invisible(dev.off())

# --- figure 4: peak signal per time point ----------------------------------
# One page per comparison: a time point that only appears in the 12 h contrast
# would otherwise be drawn with the 2 h gene set, and the boxes would no longer
# describe one set of regions.
signal_page <- function(mark, comparison) {
  frame <- signal[signal$mark == mark & signal$comparison == comparison, , drop = FALSE]
  for (region_type in region_levels) {
    selected <- frame[frame$region_type == region_type, , drop = FALSE]
    label <- sprintf("%s (%s)", panel_label(mark, comparison), region_type)
    if (nrow(selected) == 0) {
      plot.new()
      title(main = label, cex.main = 0.9)
      text(0.5, 0.5, "no regions", col = "grey40")
      next
    }
    columns <- intersect(signal_columns, colnames(selected))
    # Only the time points this contrast measured keep a box.
    here <- timepoints[timepoints %in% vapply(columns, timepoint_of, character(1))]
    values_by_direction <- list(
      up = lapply(here, function(timepoint) {
        columns_of <- columns[vapply(columns, timepoint_of, character(1)) == timepoint]
        if (length(columns_of) == 0) return(numeric(0))
        apply(selected[selected$chip_direction == "up", columns_of, drop = FALSE], 1,
              function(row) mean(as_num(row), na.rm = TRUE))
      }),
      down = lapply(here, function(timepoint) {
        columns_of <- columns[vapply(columns, timepoint_of, character(1)) == timepoint]
        if (length(columns_of) == 0) return(numeric(0))
        apply(selected[selected$chip_direction == "down", columns_of, drop = FALSE], 1,
              function(row) mean(as_num(row), na.rm = TRUE))
      })
    )
    all_values <- unlist(values_by_direction)
    all_values <- log2(all_values[!is.na(all_values)] + 1)
    if (length(all_values) == 0) {
      plot.new()
      title(main = label, cex.main = 0.9)
      text(0.5, 0.5, "no signal values", col = "grey40")
      next
    }
    limits <- range(all_values)
    plot(NA, xlim = c(0.6, length(here) + 0.4), ylim = limits, xaxt = "n",
         xlab = "", ylab = "log2(normalized counts + 1)", main = label, las = 1, bty = "n")
    axis(1, at = seq_along(here), labels = here)
    colours <- c(up = "steelblue", down = "indianred")
    for (name in names(values_by_direction)) {
      for (index in seq_along(here)) {
        values <- values_by_direction[[name]][[index]]
        values <- log2(values[!is.na(values)] + 1)
        if (length(values) == 0) next
        offset <- if (name == "up") -0.16 else 0.16
        boxplot(values, at = index + offset, add = TRUE, axes = FALSE,
                col = adjustcolor(colours[[name]], alpha.f = 0.35),
                border = colours[[name]], outline = FALSE, width = 0.28)
        points(index + offset + runif(length(values), -0.05, 0.05), values,
               pch = 16, cex = 0.4, col = adjustcolor(colours[[name]], alpha.f = 0.5))
      }
    }
    n_up <- sum(selected$chip_direction == "up")
    n_down <- sum(selected$chip_direction == "down")
    legend("topleft", bty = "n", cex = 0.72,
           fill = c(adjustcolor("steelblue", alpha.f = 0.35),
                    adjustcolor("indianred", alpha.f = 0.35)),
           border = c("steelblue", "indianred"),
           legend = c(sprintf("ChIP up (n=%d)", n_up), sprintf("ChIP down (n=%d)", n_down)))
  }
}

if (!is.null(signal) && nrow(signal) > 0) {
  signal_combos <- unique(signal[, c("mark", "comparison")])
  signal_combos <- signal_combos[order(signal_combos$mark, signal_combos$comparison), , drop = FALSE]
  pdf(file.path(opt$outdir, "box_signal.pdf"), width = 10, height = 4.6)
  par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1))
  for (index in seq_len(nrow(signal_combos))) {
    signal_page(signal_combos$mark[index], signal_combos$comparison[index])
  }
  invisible(dev.off())
}

cat(sprintf("[common] figures written to %s\n", opt$outdir))
