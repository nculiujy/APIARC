#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--stats", type = "character", help = "stats_summary.tsv"),
  make_option("--dose", type = "character", help = "dose_response.tsv"),
  make_option("--gene-table", type = "character", help = "chip_gene_table.tsv"),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$stats) || is.null(opt$outdir)) stop("missing --stats or --outdir")

as_num <- function(x) suppressWarnings(as.numeric(x))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

stats <- read.delim(opt$stats, check.names = FALSE, stringsAsFactors = FALSE)
gene <- if (!is.null(opt$`gene-table`) && file.exists(opt$`gene-table`)) {
  read.delim(opt$`gene-table`, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  NULL
}
dose <- if (!is.null(opt$dose) && file.exists(opt$dose)) {
  read.delim(opt$dose, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  NULL
}

family_of <- function(f) stats$analysis[sub(":.*$", "", stats$analysis) == f]
stratum_of <- function(family, stratum) {
  stats$analysis == sprintf("%s:%s", family, stratum)
}

# "H3K27ac|H3K27ac_TGFb12h_vs_H3K27ac_veh" -> "H3K27ac TGFb12h vs veh": the raw
# comparison string repeats the mark and does not fit the axis.
short_label <- function(comparison) {
  parts <- strsplit(comparison, "|", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(comparison)
  mark <- parts[1]
  halves <- strsplit(parts[2], "_vs_", fixed = TRUE)[[1]]
  if (length(halves) != 2) return(comparison)
  drop_mark <- function(value) sub(paste0("^", mark, "_"), "", value)
  sprintf("%s %s vs %s", mark, drop_mark(halves[1]), drop_mark(halves[2]))
}

forest <- function(frame, title, xlab = "odds ratio (95% CI, log scale)") {
  if (nrow(frame) == 0) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "no test available", col = "grey40")
    return(invisible(NULL))
  }
  frame <- frame[order(frame$effect_size), , drop = FALSE]
  low <- as_num(frame$ci_low)
  high <- as_num(frame$ci_high)
  finite <- is.finite(low) & is.finite(high) & is.finite(as_num(frame$effect_size))
  frame <- frame[finite, , drop = FALSE]
  low <- low[finite]
  high <- high[finite]
  if (nrow(frame) == 0) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "confidence intervals unavailable", col = "grey40")
    return(invisible(NULL))
  }
  or <- as_num(frame$effect_size)
  labels <- vapply(frame$comparison, short_label, character(1))
  layout_limits <- range(c(low, high, 1), na.rm = TRUE)
  plot(
    NA,
    xlim = layout_limits, ylim = c(0.5, nrow(frame) + 0.5),
    log = "x", yaxt = "n", xlab = xlab, ylab = "", bty = "n",
    main = title
  )
  abline(v = 1, lty = 2, col = "grey50")
  crossing <- low <= 1 & high >= 1
  segments(low, seq_len(nrow(frame)), high, seq_len(nrow(frame)), col = "grey30")
  points(or, seq_len(nrow(frame)), pch = 21, cex = 1.2,
         bg = ifelse(crossing, NA, "steelblue"), col = "grey20")
  axis(2, at = seq_len(nrow(frame)), labels = labels, las = 1, cex.axis = 0.7)
}

distribution <- function(frame, title) {
  categories <- c("diff_down", "diff_up")
  values <- lapply(categories, function(category) {
    as_num(frame$rna_log2FC[frame$category == category])
  })
  names(values) <- categories
  values <- lapply(values, function(v) v[!is.na(v)])
  values <- values[vapply(values, length, integer(1)) > 0]
  if (length(values) == 0) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "no differential genes", col = "grey40")
    return(invisible(NULL))
  }
  boxplot(
    values, outline = FALSE, col = c("#4C72B0", "#C44E52")[seq_along(values)],
    ylab = "RNA log2 fold change", main = title, las = 1
  )
  abline(h = 0, lty = 2, col = "grey60")
  for (i in seq_along(values)) {
    points(jitter(rep(i, length(values[[i]])), amount = 0.15), values[[i]],
           pch = 16, cex = 0.25, col = adjustcolor("grey20", alpha.f = 0.35))
  }
}

dose_panel <- function(frame, title) {
  if (is.null(frame) || nrow(frame) == 0) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "no decile data", col = "grey40")
    return(invisible(NULL))
  }
  plot(
    frame$mean_chip_log2FC, frame$mean_rna_log2FC,
    type = "n", xlab = "mean ChIP log2 fold change (decile)",
    ylab = "mean RNA log2 fold change", main = title, las = 1, bty = "n"
  )
  abline(h = 0, v = 0, lty = 2, col = "grey60")
  segments(frame$mean_chip_log2FC, frame$ci_low,
           frame$mean_chip_log2FC, frame$ci_high, col = "grey40")
  lines(frame$mean_chip_log2FC, frame$mean_rna_log2FC, col = "steelblue", lwd = 2)
  points(frame$mean_chip_log2FC, frame$mean_rna_log2FC, pch = 21, bg = "steelblue", cex = 1.1)
}

correlation_panel <- function(frame, title) {
  if (is.null(frame) || nrow(frame) < 10) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "too few genes", col = "grey40")
    return(invisible(NULL))
  }
  chip <- as_num(frame$log2FC)
  rna <- as_num(frame$rna_log2FC)
  keep <- !is.na(chip) & !is.na(rna)
  chip <- chip[keep]
  rna <- rna[keep]
  plot(
    chip, rna, pch = 16, cex = 0.35,
    col = adjustcolor("grey20", alpha.f = 0.4),
    xlab = "ChIP log2 fold change", ylab = "RNA log2 fold change",
    main = title, las = 1, bty = "n"
  )
  abline(h = 0, v = 0, lty = 2, col = "grey60")
  if (length(chip) > 2) {
    fit <- suppressWarnings(lm(rna ~ chip))
    abline(fit, col = "indianred", lwd = 2)
  }
}

pdf(file.path(opt$outdir, "fig4_evidence.pdf"), width = 10, height = 8)
if (!is.null(gene)) {
  par(mfrow = c(2, 2), mar = c(4, 9, 3, 1))
  forest(stats[stratum_of("h3_overlap", "all"), , drop = FALSE],
         "H3: differential binding x differential expression", "odds ratio (95% CI)")
  distribution(gene[gene$region_type == "promoter", , drop = FALSE],
               "H2: RNA change by binding state (promoter)")
  dose_panel(dose[dose$region_type == "promoter", , drop = FALSE],
             "H2: dose response (promoter)")
  correlation_panel(gene[gene$region_type == "promoter", , drop = FALSE],
                    "H2: ChIP vs RNA, no direction filter (promoter)")
} else {
  forest(stats[stratum_of("h3_overlap", "all"), , drop = FALSE],
         "H3: differential binding x differential expression", "odds ratio (95% CI)")
}
invisible(dev.off())

# One figure with the three region strata side by side: as three separate
# files each held only a handful of rows and duplicated the forest already in
# fig4_evidence.
pdf(file.path(opt$outdir, "fig4_forest.pdf"), width = 13, height = 5)
par(mfrow = c(1, 3), mar = c(4, 11, 3, 1))
for (stratum in c("all", "promoter", "distal")) {
  forest(stats[stratum_of("h3_overlap", stratum), , drop = FALSE],
         sprintf("H3: overlap odds ratio (%s)", stratum))
}
invisible(dev.off())

if (!is.null(gene)) {
  for (stratum in c("all", "promoter", "distal")) {
    frame <- gene[gene$region_type == stratum, , drop = FALSE]
    if (nrow(frame) == 0) next
    pdf(file.path(opt$outdir, sprintf("fig4_distribution_%s.pdf", stratum)),
        width = 7, height = 6)
    distribution(frame, sprintf("RNA change by binding state (%s)", stratum))
    invisible(dev.off())
  }
}

if (!is.null(dose) && nrow(dose) > 0) {
  groups <- split(dose, with(dose, interaction(mark, comparison, region_type,
                                               drop = TRUE, sep = "|")))
  # A comparison with no differential binding splits every decile at the same
  # ~0 ChIP fold change, so the panel is a vertical smear with no information;
  # drop those and keep only panels whose x actually spans a range.
  informative <- vapply(groups, function(frame) {
    diff(range(as_num(frame$mean_chip_log2FC), na.rm = TRUE)) > 1e-4
  }, logical(1))
  groups <- groups[informative]
  if (length(groups) > 0) {
    dose_label <- function(frame) {
      mark <- frame$mark[1]
      halves <- strsplit(frame$comparison[1], "_vs_", fixed = TRUE)[[1]]
      drop_mark <- function(value) sub(paste0("^", mark, "_"), "", value)
      short <- if (length(halves) == 2)
        sprintf("%s vs %s", drop_mark(halves[1]), drop_mark(halves[2]))
      else
        frame$comparison[1]
      sprintf("%s %s (%s)", mark, short, frame$region_type[1])
    }
    # One row per two panels: a fixed 6-inch page left each row about 0.85 inch,
    # less than the default margins, and plot.new() aborted the whole stage with
    # "figure margins too large".
    rows <- ceiling(length(groups) / 2)
    pdf(file.path(opt$outdir, "fig4_dose_response.pdf"),
        width = 10, height = max(4, 2.6 * rows))
    par(mfrow = c(rows, 2), mar = c(4, 4, 3, 1))
    for (frame in groups) {
      dose_panel(frame, dose_label(frame))
    }
    invisible(dev.off())
  }
}

cat(sprintf("[plots] figures written to %s\n", opt$outdir))
