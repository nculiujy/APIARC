#!/usr/bin/env Rscript

.required <- c("optparse", "ggplot2", "dplyr", "RColorBrewer", "pheatmap", "ggrepel")
.missing <- .required[!vapply(.required, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing) > 0) {
  stop(
    "missing R package(s): ", paste(.missing, collapse = ", "),
    "\nThey must come from the pinned conda environment; run snakemake with --use-conda.",
    call. = FALSE
  )
}

library(optparse)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)



option_list <- list(
  make_option(c("--log"), type="character",
              help="optional HISAT2 log; only used for the alignment-composition panel"),
  make_option(c("--deg"), type="character", help="analysiscsvfile"),
  make_option(c("--groups"), type="character", default=NULL,
              help="optional sample=group pairs, comma separated; colours the PCA"),
  make_option(c("--outdir_table"), type="character", help=""),
  make_option(c("--outdir_plot"), type="character", help=""),
  make_option(c("--bed"), type="character", help="bedfile")
)
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$deg) | is.null(opt$outdir_table) | is.null(opt$outdir_plot) | is.null(opt$bed)){
  stop("--deg --outdir_table --outdir_plot --bed are required")
}

if(!dir.exists(opt$outdir_table)) dir.create(opt$outdir_table, recursive=TRUE)
if(!dir.exists(opt$outdir_plot)) dir.create(opt$outdir_plot, recursive=TRUE)

bed_path <- opt$bed


log_lines <- if (!is.null(opt$log) && file.exists(opt$log)) readLines(opt$log) else character()

sample_idx <- grep("^={5,}\\s*\\S+\\s*={5,}$", log_lines)
sample_names <- sub("^=+\\s*(\\S+)\\s*=+$", "\\1", log_lines[sample_idx])

result <- data.frame()

for (i in seq_along(sample_idx)) {
  start <- sample_idx[i]
  end <- ifelse(i < length(sample_idx), sample_idx[i+1] - 1, length(log_lines))
  block <- log_lines[start:end]
  

  is_paired <- any(grepl("were paired", block))
  is_unpaired <- any(grepl("were unpaired", block))
  
  if (is_paired) {

    unique_pct <- as.numeric(sub(".*\\(([0-9.]+)%\\) aligned concordantly exactly 1 time.*", "\\1",
                                 grep("aligned concordantly exactly 1 time", block, value=TRUE)))
    multi_pct <- as.numeric(sub(".*\\(([0-9.]+)%\\) aligned concordantly >1 times.*", "\\1",
                                grep("aligned concordantly >1 times", block, value=TRUE)))
    unmapped_pct <- as.numeric(sub(".*\\(([0-9.]+)%\\) aligned concordantly 0 times.*", "\\1",
                                   grep("aligned concordantly 0 times", block, value=TRUE)))

    unmapped_pct <- unmapped_pct[1]
  } else if (is_unpaired) {

    unique_pct <- as.numeric(sub(".*\\(([0-9.]+)%\\) aligned exactly 1 time.*", "\\1",
                                 grep("aligned exactly 1 time", block, value=TRUE)))
    multi_pct <- as.numeric(sub(".*\\(([0-9.]+)%\\) aligned >1 times.*", "\\1",
                                grep("aligned >1 times", block, value=TRUE)))
    unmapped_pct <- as.numeric(sub(".*\\(([0-9.]+)%\\) aligned 0 times.*", "\\1",
                                   grep("aligned 0 times", block, value=TRUE)))
  } else {

    unique_pct <- 0
    multi_pct <- 0
    unmapped_pct <- 0
  }
  
  if(length(unique_pct) == 0 || is.na(unique_pct)) unique_pct <- 0
  if(length(multi_pct) == 0 || is.na(multi_pct)) multi_pct <- 0
  if(length(unmapped_pct) == 0 || is.na(unmapped_pct)) unmapped_pct <- 0
  
  if (length(sample_names) >= i) {
      result <- rbind(result,
                      data.frame(
                        sample = sample_names[i],
                        type = c("Unique", "Multiple", "Unmapped"),
                        percent = c(unique_pct, multi_pct, unmapped_pct)
                      )
      )
  }
}

if (nrow(result) > 0) {
    result$type <- factor(result$type, levels = c("Unique", "Multiple", "Unmapped"))

    pdf(file.path(opt$outdir_plot, "hisat2_result_triple_bar.pdf"), width=9, height=5)
    p <- ggplot(result, aes(x = sample, y = percent, fill = type)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.6) +
      geom_text(aes(label = sprintf("%.2f%%", percent)),
                position = position_dodge(width = 0.6), vjust = -0.3, size = 5) +
      scale_fill_manual(values = c("Unique" = "#4daf4a", "Multiple" = "#377eb8", "Unmapped" = "#e41a1c")) +
      ylim(0, 100) +
      labs(
        title = "Hisat2 Mapping Statistics",
        y = "Percentage (%)",
        fill = "Type"
      ) +
      theme_bw() +
      theme(
        axis.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(p)
    dev.off()
}


deg.data <- tryCatch({
  read.csv(opt$deg, row.names = NULL)
}, error = function(e) {
  data.frame()
})

if (nrow(deg.data) == 0) {
  stop("DEG result is empty: ", opt$deg,
       "\nDESeq2 produced no rows; the contrast is undefined.", call. = FALSE)
}
bed <- read.table(bed_path, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
bed_names <- c("chr", "start", "end", "gene_id", "score", "strand", "symbol")
colnames(bed) <- bed_names[seq_len(ncol(bed))]

log2FC_cutoff <- 1
padj_cutoff <- 0.05

up_genes <- dplyr::filter(deg.data, log2FoldChange > log2FC_cutoff, padj < padj_cutoff)
down_genes <- dplyr::filter(deg.data, log2FoldChange < -log2FC_cutoff, padj < padj_cutoff)

matched_up <- merge(up_genes[, c("gene_id", "SYMBOL")], bed, by="gene_id")
matched_down <- merge(down_genes[, c("gene_id", "SYMBOL")], bed, by="gene_id")
matched_up <- matched_up[, c("chr", "start", "end", "gene_id", "SYMBOL")]
matched_down <- matched_down[, c("chr", "start", "end", "gene_id", "SYMBOL")]

write.csv(up_genes, file=file.path(opt$outdir_table, "UP_genes_name.csv"), row.names=FALSE, quote=FALSE)
write.csv(down_genes, file=file.path(opt$outdir_table, "DOWN_genes_name.csv"), row.names=FALSE, quote=FALSE)
write.table(matched_up, file=file.path(opt$outdir_table, "UP_gene_names.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(matched_down, file=file.path(opt$outdir_table, "DOWN_gene_names.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


deg.data_nodup <- deg.data[!duplicated(deg.data$SYMBOL), ]
reserved <- c("SYMBOL", "gene_id", "baseMean", "log2FoldChange", "lfcSE",
              "stat", "pvalue", "padj")
count_columns <- setdiff(colnames(deg.data_nodup), reserved)
if (length(count_columns) == 0) {
  stop("no normalised count columns found in ", opt$deg, call. = FALSE)
}
expr_mat <- as.matrix(deg.data_nodup[, count_columns, drop = FALSE])
rownames(expr_mat) <- deg.data_nodup$SYMBOL
n_sample <- ncol(expr_mat)
mycol <- brewer.pal(min(n_sample, 12), "Set3")
if(n_sample > 12) mycol <- colorRampPalette(brewer.pal(12, "Set3"))(n_sample)
pdf(file.path(opt$outdir_plot, "sample_boxplot.pdf"), width=7, height=5)
boxplot(log2(expr_mat + 1), 
        main = "Sample expression distribution", 
        ylab = "log2(normalized counts + 1)", 
        col = mycol,
        xaxt = "n")
axis(1, at=1:n_sample, labels=FALSE)
labels <- colnames(expr_mat)
x <- 1:n_sample
y <- par("usr")[3] - 0.02 * (par("usr")[4] - par("usr")[3])
text(x, y, labels, srt=45, adj=1, cex=0.7, xpd=TRUE)
dev.off()


sample_dists <- dist(t(log2(expr_mat + 1)))
sample_dist_matrix <- as.matrix(sample_dists)
pdf(file.path(opt$outdir_plot, "sample_dist_heatmap.pdf"), width=7, height=6)
pheatmap(sample_dist_matrix, main = "Sample distance heatmap")
dev.off()


# A gene with the same count in every sample has zero variance and makes
# prcomp's scaling abort; PCA only needs the genes that vary.
log_mat <- log2(expr_mat + 1)
variable <- apply(log_mat, 1, var) > 0
if (sum(variable) < 2) {
  stop("no genes with variance across samples in ", opt$deg, call. = FALSE)
}
pca <- prcomp(t(log_mat[variable, , drop = FALSE]), scale. = TRUE)
pca_df <- data.frame(pca$x)
sample_groups <- NULL
if (!is.null(opt$groups) && nzchar(opt$groups)) {
  pairs <- strsplit(opt$groups, ",")[[1]]
  kv <- do.call(rbind, strsplit(pairs, "=", fixed = TRUE))
  sample_groups <- setNames(kv[, 2], kv[, 1])
}
pca_df$group <- if (is.null(sample_groups)) {
  "all"
} else {
  values <- unname(sample_groups[rownames(pca_df)])
  # Levels are per contrast here, so a sample missing from the metadata shows
  # up as unassigned instead of silently becoming background.
  ifelse(is.na(values), "unassigned", values)
}
pca_df$sample <- rownames(pca_df)
pdf(file.path(opt$outdir_plot, "sample_PCA.pdf"), width=10, height=8)
ggplot(pca_df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel() +
  theme_bw() +
  ggtitle("PCA of samples") +
  theme(strip.text = element_blank()) +
  coord_cartesian(clip = "off") +
  expand_limits(x = range(pca_df$PC1) + c(-1, 1),
                y = range(pca_df$PC2) + c(-1, 1))
dev.off()


deg.data <- deg.data %>% filter(pvalue >= 1e-30)
need_DEG <- deg.data %>%
  transmute(
    log2FoldChange = log2FoldChange,
    padj = padj,
    symbol = if ("SYMBOL" %in% names(deg.data)) as.character(SYMBOL) else as.character(gene_id)
  )
need_DEG$significance <- factor(
  ifelse(need_DEG$padj < padj_cutoff & abs(need_DEG$log2FoldChange) > log2FC_cutoff,
         ifelse(need_DEG$log2FoldChange > 0, 'UP', 'DOWN'), 'NOT'),
  levels = c('UP', 'DOWN', 'NOT')
)

# The contrast name is the last component of the plot directory
# (.../plots/treatA_vs_control); the title used to carry only the counts.
contrast_label <- basename(opt$outdir_plot)
n_up <- sum(need_DEG$significance == 'UP')
n_down <- sum(need_DEG$significance == 'DOWN')

# Label the most significant changed genes, ties broken by fold change.
sig <- need_DEG %>% filter(significance != 'NOT')
top_labels <- sig %>% arrange(padj, -abs(log2FoldChange)) %>% slice_head(n = 10)

# A handful of low-count outliers reach |log2FC| ~ 20 and squash the real
# signal into the centre; clip to a symmetric window that still holds every
# significant gene so the labelled points stay on the page.
x_edge <- if (nrow(sig) > 0) max(abs(sig$log2FoldChange)) else 2
x_lim <- min(max(x_edge, log2FC_cutoff + 1), 8)

ggplot(need_DEG, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.4) +
  scale_colour_manual(values = c(UP = '#C0392B', DOWN = '#2471A3', NOT = 'grey75'),
                      guide = 'none') +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = 2,
             color = 'grey45', linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = 2,
             color = 'grey45', linewidth = 0.5) +
  geom_text_repel(
    data = top_labels, aes(label = symbol),
    size = 3.2, max.overlaps = 20, box.padding = 0.4, min.segment.length = 0,
    segment.color = 'grey55', segment.size = 0.3
  ) +
  coord_cartesian(xlim = c(-x_lim, x_lim), clip = 'off') +
  labs(
    title = contrast_label,
    subtitle = sprintf('%d up / %d down  (|log2FC| > %g, padj < %g)',
                       n_up, n_down, log2FC_cutoff, padj_cutoff),
    x = expression(log[2]~"fold change"),
    y = expression(-log[10]~"adjusted p-value")
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    plot.subtitle = element_text(hjust = 0.5, color = 'grey30'),
    axis.title = element_text(face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  )
ggsave(file.path(opt$outdir_plot, "mRNA_volcano.pdf"), width = 7, height = 6.5)

cat("[INFO] analysis，", opt$outdir_table, ", plots output to", opt$outdir_plot, "\n")


