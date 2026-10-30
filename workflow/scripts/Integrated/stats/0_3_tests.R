#!/usr/bin/env Rscript

locate_helpers <- function() {
  if (exists("snakemake")) {
    return(file.path(snakemake@scriptdir, "..", "..", "stats", "stats_common.R"))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script <- if (length(file_arg) == 1) sub("^--file=", "", file_arg) else
    "workflow/scripts/Integrated/stats/0_3_tests.R"
  file.path(dirname(script), "..", "..", "stats", "stats_common.R")
}
source(normalizePath(locate_helpers()))

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--gene-table", type = "character",
              help = "chip_gene_table.tsv from step 2"),
  make_option("--rna", type = "character", action = "list", default = list(),
              help = "comparison=DESeq2 table, one per contrast; the universe is their union"),
  make_option("--mode", type = "character", default = "differential",
              help = "differential | descriptive"),
  # optparse keeps the dash in the option name, so every dashed flag needs an
  # explicit dest; opt$min_group and opt$neg-control would not resolve.
  make_option("--neg-control", type = "character", action = "list", default = list(),
              dest = "neg_control",
              help = "negative_control.tsv files; path encodes mark/comparison"),
  make_option("--permutations", type = "integer", default = 10000L),
  make_option("--bootstrap", type = "integer", default = 2000L),
  make_option("--seed", type = "integer", default = 20261023L),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--rna-lfc", type = "double", default = 0, dest = "rna_lfc"),
  make_option("--min-base-mean", type = "double", default = 10,
              dest = "min_base_mean"),
  make_option("--min-group", type = "integer", default = 5L, dest = "min_group"),
  make_option("--outdir", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$`gene-table`) || length(opt$rna) == 0 || is.null(opt$outdir)) {
  stop("missing --gene-table, --rna or --outdir")
}
set.seed(opt$seed)

as_num <- function(x) suppressWarnings(as.numeric(x))
na_if_blank <- function(x) ifelse(x == "", NA, x)

gene <- read.delim(opt$`gene-table`, check.names = FALSE, stringsAsFactors = FALSE)
# The tests are genome-wide, so the universe comes from every contrast; a
# missing value is filled from the table of the *same* comparison, never from
# another time point.
rna_frames <- list()
rna_contrasts <- list()
for (item in opt$rna) {
  parts <- strsplit(item, "=", fixed = TRUE)[[1]]
  label <- if (length(parts) > 1) parts[1] else ""
  path <- if (length(parts) > 1) paste(parts[-1], collapse = "=") else item
  frame <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  frame$SYMBOL <- trimws(as.character(frame$SYMBOL))
  frame <- frame[!duplicated(frame$SYMBOL), , drop = FALSE]
  rna_frames[[length(rna_frames) + 1]] <- frame
  if (nzchar(label)) {
    rna_contrasts[[label]] <- list(
      lfc = setNames(as_num(frame$log2FoldChange), frame$SYMBOL),
      padj = setNames(as_num(frame$padj), frame$SYMBOL)
    )
  }
}
rna <- do.call(rbind, rna_frames)
rna$SYMBOL <- trimws(as.character(rna$SYMBOL))
# A gene can be missing from one contrast (DESeq2 drops rows without a padj),
# so rows with a usable padj win the de-duplication.
rna <- rna[order(is.na(suppressWarnings(as.numeric(rna$padj)))), , drop = FALSE]
rna <- rna[!duplicated(rna$SYMBOL), , drop = FALSE]

rna_by_symbol <- list()
for (i in seq_len(nrow(rna))) {
  symbol <- trimws(as.character(rna$SYMBOL[i]))
  if (!is.na(symbol) && nzchar(symbol) && is.null(rna_by_symbol[[symbol]])) {
    rna_by_symbol[[symbol]] <- rna[i, ]
  }
}

expressed <- function() {
  bm <- as_num(rna$baseMean)
  !is.na(bm) & bm >= opt$min_base_mean & nzchar(trimws(as.character(rna$SYMBOL)))
}
expressed_symbols <- unique(trimws(as.character(rna$SYMBOL[expressed()])))

category_score <- function(category) {
  ifelse(category %in% c("diff_up", "diff_down"), 3L,
         ifelse(category == "discordant", 2L,
                ifelse(category %in% c("static", "bound"), 1L, 0L)))
}

collapse_genes <- function(frame) {
  if (nrow(frame) == 0) return(frame)
  magnitude <- abs(as_num(frame$log2FC))
  magnitude[is.na(magnitude)] <- 0
  frame <- frame[order(frame$gene, -category_score(frame$category), -magnitude), , drop = FALSE]
  frame[!duplicated(frame$gene), , drop = FALSE]
}

symbols <- trimws(as.character(rna$SYMBOL))
rna_lfc <- setNames(as_num(rna$log2FoldChange), symbols)
rna_padj <- setNames(as_num(rna$padj), symbols)
rna_lfc <- rna_lfc[!duplicated(names(rna_lfc))]
rna_padj <- rna_padj[!duplicated(names(rna_padj))]

marks <- unique(gene$mark)
comparisons <- unique(gene$comparison)
strata <- c("all", "promoter", "distal")
stats <- new_stats_table()
dose_rows <- list()

for (mark in marks) {
  for (comparison in comparisons) {
    slice <- gene[gene$mark == mark & gene$comparison == comparison, , drop = FALSE]
    if (nrow(slice) == 0) next
    label <- sprintf("%s|%s", mark, comparison)

    for (stratum in strata) {
      frame <- if (stratum == "all") slice else slice[slice$region_type == stratum, , drop = FALSE]
      frame <- collapse_genes(frame)
      frame <- frame[frame$gene %in% expressed_symbols, , drop = FALSE]
      if (nrow(frame) < opt$min_group) next

      frame$rna_log2FC <- as_num(frame$rna_log2FC)
      frame$rna_padj <- as_num(frame$rna_padj)
      frame$log2FC <- as_num(frame$log2FC)
      contrast <- rna_contrasts[[comparison]]
      if (!is.null(contrast)) {
        missing_lfc <- is.na(frame$rna_log2FC)
        if (any(missing_lfc)) {
          frame$rna_log2FC[missing_lfc] <- unname(contrast$lfc[frame$gene[missing_lfc]])
        }
        missing_padj <- is.na(frame$rna_padj)
        if (any(missing_padj)) {
          frame$rna_padj[missing_padj] <- unname(contrast$padj[frame$gene[missing_padj]])
        }
      }

      differential <- identical(opt$mode, "differential") &&
        any(frame$category %in% c("diff_up", "diff_down"))

      exposure <- if (differential) {
        frame$category %in% c("diff_up", "diff_down")
      } else {
        frame$category %in% c("bound", "static")
      }
      outcome <- !is.na(frame$rna_padj) & frame$rna_padj < opt$alpha &
        !is.na(frame$rna_log2FC) & abs(frame$rna_log2FC) >= opt$rna_lfc
      exposure_label <- if (differential) "differential" else "bound"

      n_exposed <- sum(exposure)
      n_unexposed <- sum(!exposure)
      n_outcome <- sum(outcome)
      if (n_exposed >= opt$min_group && n_unexposed >= opt$min_group &&
          n_outcome >= opt$min_group) {
        a <- sum(outcome & exposure)
        b <- sum(!outcome & exposure)
        cc <- sum(outcome & !exposure)
        dd <- sum(!outcome & !exposure)
        fit <- fisher_exact_ci(a, b, cc, dd)
        stats <- add_row(
          stats,
          analysis = sprintf("h3_overlap:%s", stratum),
          comparison = label,
          test = sprintf("fisher_%s_vs_rna_de", exposure_label),
          n1 = n_exposed, n2 = n_unexposed,
          statistic = fit$or, effect_size = fit$or,
          ci_low = fit$ci_low, ci_high = fit$ci_high, p = fit$p,
          note = sprintf(
            "a=%d b=%d c=%d d=%d; universe=%d genes with a %s peak and RNA baseMean>=%g",
            a, b, cc, dd, nrow(frame), mark, opt$min_base_mean
          )
        )

        n_perm <- opt$permutations
        if (n_perm > 0) {
          null_or <- replicate(n_perm, {
            drawn <- sample.int(nrow(frame), n_exposed)
            a2 <- sum(outcome[drawn])
            c2 <- n_outcome - a2
            b2 <- n_exposed - a2
            d2 <- (nrow(frame) - n_exposed) - c2
            ((a2 + 0.5) * (d2 + 0.5)) / ((b2 + 0.5) * (c2 + 0.5))
          })
          stats <- add_row(
            stats,
            analysis = sprintf("h3_permutation:%s", stratum),
            comparison = label,
            test = sprintf("label_permutation_%d", n_perm),
            n1 = n_exposed, n2 = n_unexposed,
            statistic = fit$or, effect_size = fit$or,
            p = permutation_p(log(fit$or + 1e-9), log(null_or + 1e-9)),
            note = "empirical null by reshuffling the exposure labels within the universe"
          )
        }
      }

      if (differential) {
        reference <- frame$rna_log2FC[frame$category == "static"]
        reference <- reference[!is.na(reference)]
        for (group in c("diff_up", "diff_down")) {
          values <- frame$rna_log2FC[frame$category == group]
          values <- values[!is.na(values)]
          if (length(values) < opt$min_group || length(reference) < opt$min_group) next

          wt <- suppressWarnings(wilcox.test(values, reference))
          delta <- cliff_delta(values, reference)
          ci <- bootstrap_ci(
            values,
            statistic = function(v) cliff_delta(v, reference),
            n = opt$bootstrap, alpha = opt$alpha
          )
          stats <- add_row(
            stats,
            analysis = sprintf("h2_distribution:%s", stratum),
            comparison = label,
            test = sprintf("wilcoxon_%s_vs_static", group),
            n1 = length(values), n2 = length(reference),
            statistic = unname(wt$statistic), effect_size = delta,
            ci_low = unname(ci[1]), ci_high = unname(ci[2]), p = wt$p.value,
            note = sprintf(
              "median RNA log2FC %.3f vs %.3f; Cliff's delta > 0 means the %s genes sit higher",
              median(values), median(reference), group
            )
          )

          ks <- suppressWarnings(ks.test(values, reference))
          stats <- add_row(
            stats,
            analysis = sprintf("h2_ks:%s", stratum),
            comparison = label,
            test = sprintf("ks_%s_vs_static", group),
            n1 = length(values), n2 = length(reference),
            statistic = unname(ks$statistic), p = ks$p.value,
            note = "two-sample Kolmogorov-Smirnov on RNA log2FC; sensitive to any shape change"
          )
        }

        usable <- !is.na(frame$log2FC) & !is.na(frame$rna_log2FC)
        dose <- frame[usable, , drop = FALSE]
        if (nrow(dose) >= 100 && length(unique(dose$log2FC)) >= 10) {
          breaks <- unique(quantile(dose$log2FC, probs = seq(0, 1, 0.1), na.rm = TRUE))
          if (length(breaks) >= 4) {
            decile <- cut(dose$log2FC, breaks = breaks, include.lowest = TRUE, labels = FALSE)
            for (bin in sort(unique(decile))) {
              values <- dose$rna_log2FC[decile == bin]
              ci <- bootstrap_ci(values, n = opt$bootstrap, alpha = opt$alpha)
              dose_rows[[length(dose_rows) + 1]] <- data.frame(
                mark = mark, comparison = comparison, region_type = stratum,
                decile = bin, n = length(values),
                mean_chip_log2FC = mean(dose$log2FC[decile == bin]),
                mean_rna_log2FC = mean(values),
                ci_low = unname(ci[1]), ci_high = unname(ci[2]),
                stringsAsFactors = FALSE
              )
            }
            trend <- suppressWarnings(cor.test(decile, dose$rna_log2FC, method = "spearman"))
            stats <- add_row(
              stats,
              analysis = sprintf("h2_dose:%s", stratum),
              comparison = label,
              test = "spearman_decile_trend",
              n1 = nrow(dose), n2 = length(unique(decile)),
              statistic = unname(trend$estimate), p = trend$p.value,
              note = "RNA log2FC as a function of ChIP log2FC decile; monotonic trend supports a dose response"
            )
          }
        }

        if (nrow(dose) >= 50) {
          rho <- suppressWarnings(cor(dose$log2FC, dose$rna_log2FC, method = "spearman"))
          n_perm_corr <- min(opt$permutations, 5000L)
          null_rho <- replicate(n_perm_corr, {
            suppressWarnings(cor(rank(sample(dose$log2FC)), rank(dose$rna_log2FC)))
          })
          ci <- bootstrap_ci(
            seq_len(nrow(dose)),
            statistic = function(i) {
              suppressWarnings(cor(dose$log2FC[i], dose$rna_log2FC[i], method = "spearman"))
            },
            n = min(opt$bootstrap, 500L), alpha = opt$alpha
          )
          stats <- add_row(
            stats,
            analysis = sprintf("h2_correlation:%s", stratum),
            comparison = label,
            test = sprintf("spearman_all_genes_permutation_%d", n_perm_corr),
            n1 = nrow(dose),
            statistic = rho, effect_size = rho,
            ci_low = unname(ci[1]), ci_high = unname(ci[2]),
            p = permutation_p(rho, null_rho),
            note = "no direction pre-filter; the whole analysable gene set enters the test"
          )

          same_direction <- (dose$log2FC > 0 & dose$rna_log2FC > 0) |
            (dose$log2FC < 0 & dose$rna_log2FC < 0)
          if (sum(same_direction) >= 20) {
            biased <- suppressWarnings(cor(
              dose$log2FC[same_direction], dose$rna_log2FC[same_direction],
              method = "pearson"
            ))
            stats <- add_row(
              stats,
              analysis = sprintf("h2_correlation_selective:%s", stratum),
              comparison = label,
              test = "pearson_same_direction_only",
              n1 = sum(same_direction),
              statistic = biased, effect_size = biased, p = NA_real_,
              note = "BIAS REFERENCE, not evidence: reproduces the previous same-direction pre-filter"
            )
          }
        }

        up <- frame$rna_log2FC[frame$category == "diff_up"]
        down <- frame$rna_log2FC[frame$category == "diff_down"]
        static <- reference
        up <- up[!is.na(up)]
        down <- down[!is.na(down)]
        if (length(up) >= opt$min_group && length(static) >= opt$min_group) {
          fit <- fisher_exact_ci(
            sum(up > 0), sum(up <= 0),
            sum(static > 0), sum(static <= 0)
          )
          stats <- add_row(
            stats,
            analysis = sprintf("h2_direction:%s", stratum),
            comparison = label,
            test = "fisher_diff_up_positive_fraction",
            n1 = length(up), n2 = length(static),
            statistic = fit$or, effect_size = fit$or,
            ci_low = fit$ci_low, ci_high = fit$ci_high, p = fit$p,
            note = sprintf(
              "P(RNA up | ChIP up)=%.3f vs P(RNA up | static)=%.3f",
              mean(up > 0), mean(static > 0)
            )
          )
        }
        if (length(down) >= opt$min_group && length(static) >= opt$min_group) {
          fit <- fisher_exact_ci(
            sum(down < 0), sum(down >= 0),
            sum(static < 0), sum(static >= 0)
          )
          stats <- add_row(
            stats,
            analysis = sprintf("h2_direction:%s", stratum),
            comparison = label,
            test = "fisher_diff_down_negative_fraction",
            n1 = length(down), n2 = length(static),
            statistic = fit$or, effect_size = fit$or,
            ci_low = fit$ci_low, ci_high = fit$ci_high, p = fit$p,
            note = sprintf(
              "P(RNA down | ChIP down)=%.3f vs P(RNA down | static)=%.3f",
              mean(down < 0), mean(static < 0)
            )
          )
        }

        moved <- frame$category %in% c("diff_up", "diff_down") &
          !is.na(frame$rna_log2FC) & !is.na(frame$log2FC)
        if (sum(moved) >= 10) {
          chip_up <- frame$log2FC[moved] > 0
          rna_up <- frame$rna_log2FC[moved] > 0
          table_up <- table(factor(chip_up, c(FALSE, TRUE)), factor(rna_up, c(FALSE, TRUE)))
          mcnemar <- suppressWarnings(mcnemar.test(table_up))
          stats <- add_row(
            stats,
            analysis = sprintf("h2_direction:%s", stratum),
            comparison = label,
            test = "mcnemar_up_down_symmetry",
            n1 = sum(moved),
            statistic = unname(mcnemar$statistic), p = mcnemar$p.value,
            note = "discordant (up/down) pairs should be symmetric if binding and expression are unrelated"
          )
        }
      } else {
        unbound <- setdiff(expressed_symbols, frame$gene)
        unbound_lfc <- unname(rna_lfc[unbound])
        unbound_lfc <- unbound_lfc[!is.na(unbound_lfc)]
        bound <- frame$rna_log2FC[!is.na(frame$rna_log2FC)]
        if (length(bound) >= opt$min_group && length(unbound_lfc) >= opt$min_group) {
          wt <- suppressWarnings(wilcox.test(bound, unbound_lfc))
          delta <- cliff_delta(bound, unbound_lfc)
          ci <- bootstrap_ci(
            bound,
            statistic = function(v) cliff_delta(v, unbound_lfc),
            n = opt$bootstrap, alpha = opt$alpha
          )
          stats <- add_row(
            stats,
            analysis = sprintf("h1_bound_vs_unbound:%s", stratum),
            comparison = label,
            test = "wilcoxon_bound_vs_unbound",
            n1 = length(bound), n2 = length(unbound_lfc),
            statistic = unname(wt$statistic), effect_size = delta,
            ci_low = unname(ci[1]), ci_high = unname(ci[2]), p = wt$p.value,
            note = "descriptive mode: expression shift of bound genes against the rest of the transcriptome"
          )
        }
      }
    }
  }
}

for (path in opt$neg_control) {
  if (!file.exists(path)) next
  control <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  comparison_dir <- basename(dirname(path))
  mark_dir <- basename(dirname(dirname(path)))
  for (i in seq_len(nrow(control))) {
    stats <- add_row(
      stats,
      analysis = "negative_control",
      comparison = sprintf("%s|%s", mark_dir, comparison_dir),
      test = as.character(control$control[i]),
      statistic = as_num(control$max_false_positive[i]),
      effect_size = as_num(control$median_false_positive[i]),
      p = NA_real_,
      note = sprintf(
        "status=%s threshold=%s %s",
        control$status[i], format(control$threshold[i]), control$note[i]
      )
    )
  }
}

stats$family <- sub(":.*$", "", stats$analysis)
stats$family[stats$analysis == "negative_control"] <- "negative_control"

# The outdir is the shared 0_stats folder: it also holds regions.tsv,
# region_sets.tsv and the BED sets that enrichment, motif and enhancer read.
# Wiping it deleted those upstream outputs in the middle of a pipeline run, so
# the three tables written here are simply overwritten.
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
write_stats(stats, file.path(opt$outdir, "stats_summary.tsv"), family = "family")

dose_table <- if (length(dose_rows) > 0) {
  do.call(rbind, dose_rows)
} else {
  data.frame(
    mark = character(), comparison = character(), region_type = character(),
    decile = integer(), n = integer(), mean_chip_log2FC = numeric(),
    mean_rna_log2FC = numeric(), ci_low = numeric(), ci_high = numeric()
  )
}
write.table(
  dose_table, file.path(opt$outdir, "dose_response.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

families <- data.frame(
  family = c(
    "h3_overlap", "h3_permutation", "h2_distribution", "h2_ks", "h2_dose",
    "h2_correlation", "h2_correlation_selective", "h2_direction",
    "h1_bound_vs_unbound", "negative_control"
  ),
  description = c(
    "Fisher exact test of ChIP differential binding against RNA differential expression",
    "empirical null of the same odds ratio by reshuffling exposure labels",
    "Wilcoxon rank-sum of RNA log2FC, differential vs static, with Cliff's delta",
    "Kolmogorov-Smirnov of the same two distributions",
    "trend of RNA log2FC across ChIP log2FC deciles",
    "genome-wide Spearman correlation without direction pre-filter",
    "selection-biased reference computed on same-direction genes only; never cite",
    "enrichment of the expected direction among differential genes",
    "descriptive mode: bound vs unbound expression shift",
    "negative controls; no p value, status column in note"
  ),
  stringsAsFactors = FALSE
)
write.table(
  families, file.path(opt$outdir, "analysis_families.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat(sprintf(
  "[stats] %d rows, %d families -> %s\n",
  nrow(stats), length(unique(stats$family)), opt$outdir
))
