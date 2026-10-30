#!/usr/bin/env Rscript

bh_adjust <- function(p) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  keep <- !is.na(p)
  out[keep] <- p.adjust(p[keep], method = "BH")
  out
}

fisher_exact_ci <- function(a, b, c, d) {
  test <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  estimate <- unname(test$estimate)
  ci <- unname(test$conf.int)
  list(p = test$p.value, or = estimate, ci_low = ci[1], ci_high = ci[2])
}

cliff_delta <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  n <- length(x)
  m <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  u <- sum(ranks[seq_len(n)]) - n * (n + 1) / 2
  2 * u / (n * m) - 1
}

bootstrap_ci <- function(values, statistic = mean, n = 2000, alpha = 0.05) {
  values <- values[!is.na(values)]
  if (length(values) < 3) return(c(NA_real_, NA_real_))
  draws <- replicate(n, statistic(sample(values, replace = TRUE)))
  quantile(draws, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
}

permutation_p <- function(observed, null, alternative = "two.sided") {
  null <- null[!is.na(null)]
  if (length(null) == 0) return(NA_real_)
  if (alternative == "greater") {
    (sum(null >= observed) + 1) / (length(null) + 1)
  } else if (alternative == "less") {
    (sum(null <= observed) + 1) / (length(null) + 1)
  } else {
    (sum(abs(null) >= abs(observed)) + 1) / (length(null) + 1)
  }
}

new_stats_table <- function() {
  data.frame(
    analysis = character(), comparison = character(), test = character(),
    n1 = integer(), n2 = integer(), statistic = numeric(),
    effect_size = numeric(), ci_low = numeric(), ci_high = numeric(),
    p = numeric(), note = character(), stringsAsFactors = FALSE
  )
}

add_row <- function(table, analysis, comparison, test, n1 = NA_integer_, n2 = NA_integer_,
                    statistic = NA_real_, effect_size = NA_real_, ci_low = NA_real_,
                    ci_high = NA_real_, p = NA_real_, note = "") {
  rbind(
    table,
    data.frame(
      analysis = analysis, comparison = comparison, test = test,
      n1 = as.integer(n1), n2 = as.integer(n2),
      statistic = as.numeric(statistic), effect_size = as.numeric(effect_size),
      ci_low = as.numeric(ci_low), ci_high = as.numeric(ci_high),
      p = as.numeric(p), note = note, stringsAsFactors = FALSE
    )
  )
}

STATS_COLUMNS <- c(
  "analysis", "comparison", "test", "n1", "n2", "statistic",
  "effect_size", "ci_low", "ci_high", "p", "FDR", "note"
)

write_stats <- function(table, path, family = NULL) {
  table$FDR <- NA_real_
  if (nrow(table) > 0) {
    groups <- if (is.null(family)) factor(rep("all", nrow(table))) else factor(table[[family]])
    for (level in levels(groups)) {
      index <- which(groups == level)
      table$FDR[index] <- bh_adjust(table$p[index])
    }
  }
  table <- table[, STATS_COLUMNS, drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(table)
}
