#!/usr/bin/env Rscript

# Bipartite network drawing shared by the pathway figure and the TF-gene figure.
#
# Both networks are one hop deep and small after the filters, so a two-column
# layout is easier to read than a force-directed one: the left column is the
# regulator (term or TF), the right column the genes, and the vertical order on
# each side is by degree. base graphics is not a limitation of taste here - the
# statistics environment ships r-base without igraph or ggraph.

read_tsv_plain <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
}

# One column of nodes, ordered by degree, spread evenly over the drawing area.
node_positions <- function(degrees, top, bottom) {
  order_index <- order(-degrees, names(degrees))
  labels <- names(degrees)[order_index]
  values <- degrees[order_index]
  n <- length(labels)
  positions <- if (n == 0) numeric(0) else if (n == 1) mean(c(bottom, top)) else
    seq(top, bottom, length.out = n)
  list(labels = labels, degree = values, y = positions)
}

# Labels are thinned rather than overlapped: with 40 genes on one side every
# third name is already more than a page can hold legibly.
thin_labels <- function(labels, positions, max_labels) {
  n <- length(labels)
  # An absent limit must not abort the whole page: reduce to "label everything".
  if (length(max_labels) != 1 || is.na(max_labels) || n <= max_labels) return(rep(TRUE, n))
  keep <- unique(round(seq(1, n, length.out = max_labels)))
  seq_len(n) %in% keep
}

# Term wording runs well past what the left margin holds, and the full string is
# in the table beside the figure; cutting it keeps the page readable.
trim_text <- function(values, width = 34) {
  too_long <- !is.na(values) & nchar(values) > width
  values[too_long] <- paste0(substr(values[too_long], 1, width - 3), "...")
  values
}

draw_bipartite <- function(left, right, edge_list, title, subtitle = "",
                           left_col = "#2f6f9f", right_col = "#a5561c",
                           label_left = TRUE, max_left_labels = 12,
                           max_right_labels = 40, cex_left = 0.7,
                           cex_right = 0.55, xlim = c(-0.05, 1.05),
                           edge_alpha = c(0.08, 0.42), edge_lwd = c(0.2, 1.1)) {
  graphics::plot(NA, xlim = xlim, ylim = c(0, 1), axes = FALSE,
                 xlab = "", ylab = "", main = "")
  graphics::mtext(title, side = 3, line = 2.2, adj = 0, cex = 0.86, font = 2)
  if (nzchar(subtitle)) {
    graphics::mtext(subtitle, side = 3, line = 0.7, adj = 0, cex = 0.7, col = "grey30")
  }

  left_y <- stats::setNames(left$y, left$labels)
  right_y <- stats::setNames(right$y, right$labels)
  left_text <- if (!is.null(left$text)) left$text else left$labels
  right_text <- if (!is.null(right$text)) right$text else right$labels

  if (nrow(edge_list) > 0) {
    weight <- edge_list$weight
    scaled <- if (max(weight) > min(weight)) {
      (weight - min(weight)) / (max(weight) - min(weight))
    } else {
      rep(0.5, length(weight))
    }
    for (i in seq_len(nrow(edge_list))) {
      from <- left_y[[edge_list$left[i]]]
      to <- right_y[[edge_list$right[i]]]
      if (is.null(from) || is.null(to)) next
      graphics::segments(0.02, from, 0.98, to,
                         col = grDevices::adjustcolor(
                           "#4a4a4a",
                           alpha.f = edge_alpha[1] +
                             (edge_alpha[2] - edge_alpha[1]) * scaled[i]
                         ),
                         lwd = edge_lwd[1] + (edge_lwd[2] - edge_lwd[1]) * scaled[i])
    }
  }

  left_cex <- 0.5 + 0.14 * sqrt(left$degree / max(1, max(left$degree)))
  graphics::points(rep(0, length(left$y)), left$y, pch = 21,
                   bg = left_col, col = left_col, cex = left_cex)
  graphics::points(rep(1, length(right$y)), right$y, pch = 21,
                   bg = right_col, col = right_col,
                   cex = 0.45 + 0.12 * sqrt(right$degree / max(1, max(right$degree))))

  if (label_left) {
    keep <- thin_labels(left$labels, left$y, max_left_labels)
    graphics::text(-0.035, left$y[keep], left_text[keep], pos = 2,
                   cex = cex_left, col = left_col, xpd = NA)
  }
  keep <- thin_labels(right$labels, right$y, max_right_labels)
  graphics::text(1.035, right$y[keep], right_text[keep], pos = 4,
                 cex = cex_right, col = right_col, xpd = NA)
}

# degrees of a set of ids in one edge column
degree_of <- function(values, universe) {
  counts <- table(factor(values, levels = universe))
  as.integer(counts)
}
