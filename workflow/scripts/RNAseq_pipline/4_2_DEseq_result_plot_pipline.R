#!/usr/bin/env Rscript


if(!require("optparse")) install.packages("optparse", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("ggplot2")) install.packages("ggplot2", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("dplyr")) install.packages("dplyr", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("RColorBrewer")) install.packages("RColorBrewer", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("pheatmap")) install.packages("pheatmap", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if(!require("ggrepel")) install.packages("ggrepel", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

library(optparse)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)


# bed_path <- "/home/yxiaobo/RCProj/tfanno/gencode.vM25.mRNA.annotation.bed"

option_list <- list(
  make_option(c("--log"), type="character", help="hisat2 logfile"),
  make_option(c("--deg"), type="character", help="analysiscsvfile"),
  make_option(c("--outdir_table"), type="character", help=""),
  make_option(c("--outdir_plot"), type="character", help=""),
  make_option(c("--bed"), type="character", help="bedfile")
)
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$log) | is.null(opt$deg) | is.null(opt$outdir_table) | is.null(opt$outdir_plot) | is.null(opt$bed)){
  stop("--log --deg --outdir_table --outdir_plot --bed！")
}

if(!dir.exists(opt$outdir_table)) dir.create(opt$outdir_table, recursive=TRUE)
if(!dir.exists(opt$outdir_plot)) dir.create(opt$outdir_plot, recursive=TRUE)

bed_path <- opt$bed


log_lines <- readLines(opt$log)

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
  
  # Ensure variables are of length 1, if they are empty set them to 0
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
  cat("[WARNING] DEG result is empty! Exiting gracefully.\n")
  file.create(file.path(opt$outdir_table, "UP_genes_name.csv"))
  file.create(file.path(opt$outdir_table, "DOWN_genes_name.csv"))
  file.create(file.path(opt$outdir_table, "UP_gene_names.bed"))
  file.create(file.path(opt$outdir_table, "DOWN_gene_names.bed"))
  quit(status = 0)
}
bed <- read.table(bed_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(bed) <- c('chr', 'start', 'end', 'gene_id', 'score', 'strand')

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
expr_mat <- deg.data_nodup[, 8:ncol(deg.data_nodup)]
expr_mat <- expr_mat[, 1:4]
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


pca <- prcomp(t(log2(expr_mat + 1)), scale. = TRUE)
pca_df <- data.frame(pca$x)
pca_df$group <- rep(c("Veh", "DOX"), each = 2)
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
need_DEG <- deg.data[, c("log2FoldChange", "padj")]
colnames(need_DEG) <- c('log2FoldChange', 'padj')
need_DEG$significance <- as.factor(ifelse(
  need_DEG$padj < padj_cutoff & abs(need_DEG$log2FoldChange) > log2FC_cutoff,
  ifelse(need_DEG$log2FoldChange > log2FC_cutoff, 'UP', 'DOWN'),
  'NOT'
))
title <- paste0(' Up :  ', nrow(need_DEG[need_DEG$significance == 'UP',]),
                '\n Down : ', nrow(need_DEG[need_DEG$significance == 'DOWN',]))
g <- ggplot(data = need_DEG, 
            aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.4, size = 1) +
  theme_classic() +
  xlab("log2 ( FoldChange )") +
  ylab("-log10 ( P.adjust )") +
  ggtitle(title) +
  scale_colour_manual(values = c('blue', 'grey', 'red')) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "grey", lwd = 0.8) +
  geom_hline(yintercept = -log10(padj_cutoff), lty = 4, col = "grey", lwd = 0.8) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin = unit(c(2, 2, 2, 2), 'lines'),
        legend.title = element_blank(),
        legend.position = "right")
ggsave(g, filename = file.path(opt$outdir_plot, "mRNA_volcano.pdf"), width = 8, height = 7.5)

cat("[INFO] analysis，", opt$outdir_table, ", plots output to", opt$outdir_plot, "\n")


# Rscript DEseq_result_plot_pipline.R --log /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632/hisat2_result.log --deg /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632/mRNA/DEG_result.csv   --outdir_table /home/yxiaobo/RCProj/RNAseq/database/TF-binding/GSE85632/mRNA --outdir_plot /home/yxiaobo/RCProj/picture/RNAseq/DOX