#!/usr/bin/env Rscript

# Shared enrichment helpers.
#
# Two drivers use this file: go_kegg.R (one gene list per region set) and
# common_enrichment.R (common genes against all differential genes). The odds
# ratio lives here because the earlier version was wrong in the same way in both
# places: it took the number of submitted symbols as the gene-list size, while
# clusterProfiler drops symbols it cannot map and divides GeneRatio by that
# mapped count. The mismatch produced 63,887 rows whose overlap exceeded the
# term size, and odds ratios of -1e18.

# --flags that repeat -------------------------------------------------------
# optparse 1.7.5 builds `action = "list"` options with an empty metavar, and
# where it does parse them only the last value survives: `--rna a --rna b` came
# through as "b", so the universe was built from a single contrast. The repeated
# flags are therefore read from the command line directly.
flag_values <- function(flag, argv = commandArgs(trailingOnly = TRUE)) {
  starts <- which(argv == flag)
  if (length(starts) == 0) return(character())
  values <- character()
  for (start in starts) {
    index <- start + 1L
    while (index <= length(argv) && !grepl("^--", argv[[index]])) {
      values <- c(values, argv[[index]])
      index <- index + 1L
    }
  }
  unique(values[nzchar(values)])
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  num <- vapply(parts, function(p) as.numeric(p[1]), numeric(1))
  den <- vapply(parts, function(p) as.numeric(p[2]), numeric(1))
  list(num = num, den = den)
}

# 2x2 table for over-representation:
#   a = overlap, b = query genes outside the term, c = background genes in the
#   term, d = the rest of the background.
#
# A zero cell makes the ratio either infinite or exactly zero, so every cell is
# nudged by 0.5 (Haldane-Anscombe) and the row is flagged; downstream tables can
# then distinguish a corrected ratio from a measured one.
odds_ratio_table <- function(overlap, query_size, term_size, background_size) {
  a <- overlap
  b <- query_size - overlap
  c <- term_size - overlap
  d <- background_size - term_size - b
  invalid <- a <= 0 | b < 0 | c < 0 | d < 0
  corrected <- !invalid & (b == 0 | c == 0 | d == 0)
  a[corrected] <- a[corrected] + 0.5
  b[corrected] <- b[corrected] + 0.5
  c[corrected] <- c[corrected] + 0.5
  d[corrected] <- d[corrected] + 0.5
  value <- (a * d) / (b * c)
  value[invalid] <- NA_real_
  list(odds_ratio = value, corrected = corrected)
}

# --universe of the two drivers, from the DESeq2 tables ---------------------
# Only the two columns the universe needs are kept: the contrasts were run on
# different sample sets, so the tables do not share column names and cannot be
# row-bound whole.
read_expression <- function(path) {
  frame <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(
    SYMBOL = trimws(as.character(frame$SYMBOL)),
    baseMean = suppressWarnings(as.numeric(frame$baseMean)),
    stringsAsFactors = FALSE
  )
}

# Kept identical for every set: the union with the query was what made the
# background depend on the set being tested.
universe_from_rna <- function(paths, min_base_mean = 10) {
  paths <- unique(paths)
  if (length(paths) == 0) stop("no DESeq2 table was given for the universe")
  rna <- do.call(rbind, lapply(paths, read_expression))
  rna <- rna[nzchar(rna$SYMBOL) & !is.na(rna$SYMBOL), , drop = FALSE]
  if (nrow(rna) == 0) stop("none of the DESeq2 tables carries a SYMBOL column")
  # A symbol appears once per contrast; it is expressed when at least one
  # contrast measured it above the cut-off.
  best <- tapply(rna$baseMean, rna$SYMBOL, function(value) {
    if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE)
  })
  list(
    all = names(best),
    expressed = names(best)[!is.na(best) & best >= min_base_mean]
  )
}

to_entrez <- function(symbols, orgdb) {
  if (length(symbols) == 0) return(character())
  mapped <- suppressWarnings(tryCatch(
    clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb),
    error = function(e) NULL
  ))
  if (is.null(mapped) || nrow(mapped) == 0) return(character())
  unique(as.character(mapped$ENTREZID))
}

frame_of <- function(result, kind, label, n_query) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) return(NULL)
  frame <- as.data.frame(result)
  gene_ratio <- parse_ratio(frame$GeneRatio)
  bg_ratio <- parse_ratio(frame$BgRatio)
  overlap <- gene_ratio$num
  query_size <- gene_ratio$den      # symbols clusterProfiler actually mapped
  term_size <- bg_ratio$num
  background_size <- bg_ratio$den
  ratio <- odds_ratio_table(overlap, query_size, term_size, background_size)
  data.frame(
    set = label, ontology = kind, term_id = frame$ID, term = frame$Description,
    n_query = n_query, n_genes = query_size, n_term = term_size,
    overlap = overlap,
    odds_ratio = ratio$odds_ratio, odds_corrected = ratio$corrected,
    p = frame$pvalue, FDR = frame$p.adjust, genes = frame$geneID,
    stringsAsFactors = FALSE
  )
}

# enrichGO + enrichKEGG for one gene list against a fixed universe.
#
# The query is intersected with the universe first. clusterProfiler counts every
# input symbol in GeneRatio's denominator, including symbols the background does
# not have, which makes the 2x2 table impossible (d < 0) for any set larger than
# the universe - the H3K4me1 static set has 27,822 genes against 16,359
# expressed, and every one of its rows came out all-NA.
enrich_one <- function(genes, universe_genes, label, orgdb, kegg_organism) {
  n_submitted <- length(genes)
  if (!is.null(universe_genes)) genes <- intersect(genes, universe_genes)
  n_query <- length(genes)
  if (n_query == 0) return(NULL)

  go <- tryCatch(
    clusterProfiler::enrichGO(
      gene = genes, OrgDb = orgdb, keyType = "SYMBOL", ont = "BP",
      universe = universe_genes, pAdjustMethod = "BH",
      pvalueCutoff = 1, qvalueCutoff = 1, readable = FALSE
    ),
    error = function(e) NULL
  )
  rows <- list(frame_of(go, "BP", label, n_query))

  # enrichKEGG counts on Entrez ids; handing it symbols made every set end in
  # "No gene can be mapped".
  kegg <- tryCatch(
    {
      entrez <- to_entrez(genes, orgdb)
      if (length(entrez) == 0) stop("no Entrez id for this set")
      clusterProfiler::enrichKEGG(
        gene = entrez, organism = kegg_organism, keyType = "kegg",
        universe = to_entrez(universe_genes, orgdb), pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1
      )
    },
    error = function(e) NULL
  )
  rows[["kegg"]] <- frame_of(kegg, "KEGG", label, n_query)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  frame <- do.call(rbind, rows)
  frame$n_submitted <- n_submitted
  frame
}

# term -> gene edges of the significant terms, for the bipartite network.
#
# `max_terms` is per ontology, not in total: with one shared budget the BP half
# of the table (hundreds of significant terms) always wins and KEGG never
# reaches the figure. Within one ontology the list is also thinned by gene-set
# overlap, because "blood vessel morphogenesis", "angiogenesis" and "regulation
# of angiogenesis" are three rows of one claim and would otherwise fill the
# page with the same edges.
term_edges <- function(table, alpha, max_terms, ontology_filter = NULL,
                       jaccard_cut = 0.6) {
  if (is.null(table) || nrow(table) == 0) return(NULL)
  selected <- table
  if (!is.null(ontology_filter)) {
    selected <- selected[selected$ontology %in% ontology_filter, , drop = FALSE]
  }
  selected <- selected[!is.na(selected$FDR) & selected$FDR < alpha, , drop = FALSE]
  if (nrow(selected) == 0) return(NULL)
  selected <- selected[order(selected$FDR), , drop = FALSE]

  genes_of <- function(value) {
    genes <- strsplit(as.character(value), "/", fixed = TRUE)[[1]]
    genes[nzchar(genes)]
  }

  pick <- function(part) {
    part <- part[order(part$FDR), , drop = FALSE]
    sets <- lapply(part$genes, genes_of)
    kept <- integer()
    for (i in seq_len(nrow(part))) {
      if (length(kept) >= max_terms) break
      if (length(kept) > 0) {
        similar <- max(vapply(kept, function(j) {
          shared <- length(intersect(sets[[i]], sets[[j]]))
          total <- length(union(sets[[i]], sets[[j]]))
          if (total == 0) 0 else shared / total
        }, numeric(1)))
        if (similar > jaccard_cut) next
      }
      kept <- c(kept, i)
    }
    part[kept, , drop = FALSE]
  }

  keys <- unique(selected[, c("mark", "comparison", "region_type")])
  out <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, ]
    subset <- selected[
      selected$mark == key$mark & selected$comparison == key$comparison &
        selected$region_type == key$region_type, , drop = FALSE
    ]
    chosen <- do.call(rbind, lapply(
      split(subset, subset$ontology), pick
    ))
    if (is.null(chosen) || nrow(chosen) == 0) next
    chosen <- chosen[order(chosen$FDR), , drop = FALSE]
    for (j in seq_len(nrow(chosen))) {
      genes <- genes_of(chosen$genes[j])
      if (length(genes) == 0) next
      out[[length(out) + 1]] <- data.frame(
        mark = key$mark, comparison = key$comparison, region_type = key$region_type,
        ontology = chosen$ontology[j], term_id = chosen$term_id[j],
        term = chosen$term[j], gene = genes, FDR = chosen$FDR[j],
        odds_ratio = chosen$odds_ratio[j], stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}
