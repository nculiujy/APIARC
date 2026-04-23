#!/usr/bin/env Rscript
rm(list = ls())

install_and_load <- function(package, repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") {
  if (!requireNamespace(package, quietly = TRUE)) {
    if (package %in% c("org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler", "enrichplot")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = repo)
      }
      BiocManager::install(package, update = FALSE)
    } else {
      install.packages(package, repos = repo, dependencies = TRUE)
    }
  }
  library(package, character.only = TRUE)
}

required_packages <- c(
  "dplyr", "org.Mm.eg.db", "org.Hs.eg.db", "VennDiagram", "grid", "ggplot2",
  "clusterProfiler", "enrichplot", "ggrepel", "yaml", "optparse", "futile.logger"
)
lapply(required_packages, install_and_load)

option_list <- list(
  make_option(c("-r", "--RNAseq"), type = "character"),
  make_option(c("-c", "--ChIPseq"), type = "character"),
  make_option(c("-o", "--output"), type = "character"),
  make_option(c("-p", "--picture"), type = "character"),
  make_option(c("-y", "--yaml"), type = "character"),
  make_option(c("-m", "--matrix"), type = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

RNAseq_gene_file <- opt$RNAseq
ChIPseq_gene_file <- opt$ChIPseq
output_folder <- opt$output
outpic <- opt$picture
yaml_file <- opt$yaml
matrix_file <- opt$matrix

if (!file.exists(yaml_file)) q()
if (!file.exists(matrix_file)) q()
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
if (!dir.exists(outpic)) dir.create(outpic, recursive = TRUE)

rna_file <- file.path(RNAseq_gene_file)
if (!file.exists(rna_file)) q()
RNAseq_gene <- tryCatch({
  read.csv(rna_file, stringsAsFactors = FALSE)
}, error = function(e) {
  data.frame()
})

if (nrow(RNAseq_gene) == 0) {
    cat("[WARNING] RNAseq_gene is empty. Exiting gracefully.\n")
    q(status = 0)
}
RNAseq_gene <- RNAseq_gene[!is.na(RNAseq_gene$padj), ]

yaml_info <- yaml.load_file(yaml_file)





genome_type <- "mm" # default
if (!is.null(yaml_info$RNAseq) && !is.null(yaml_info$RNAseq$species)) {
  genome_type <- yaml_info$RNAseq$species
}

anno_file <- switch(
  genome_type,
  mm = "workflow/resources/anno/RNAseq_anno/GRCm38/gencode.vM25.annotation.bed",
  homo = "workflow/resources/anno/RNAseq_anno/GRCh38/gencode.v44.annotation.bed",
  q()
)
if (!file.exists(anno_file)) {
    cat("Annotation file not found:", anno_file, "\n")
    q()
}
anno <- read.table(anno_file, header = FALSE, stringsAsFactors = FALSE)



chip_dirs <- list.dirs(ChIPseq_gene_file, full.names = TRUE, recursive = FALSE)
chip_csv_files <- file.path(chip_dirs, paste0(basename(chip_dirs), "_promoter_peaks.csv"))
# The directories might be named "H3K4me1_vs_Input", but matrix_list.txt uses "H3K4me1". 
# So we strip "_vs_Input" to match the sample names in matrix_list.txt
names(chip_csv_files) <- sub("_vs_Input", "", basename(chip_dirs))

chip_list <- lapply(names(chip_csv_files), function(sample_name) {
  f <- chip_csv_files[sample_name]
  if (file.exists(f)) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    return(df)
  } else {
    return(NULL)
  }
})
names(chip_list) <- names(chip_csv_files)
chip_list <- Filter(Negate(is.null), chip_list)
if (length(chip_list) == 0) {
    cat("No ChIPseq promoter peaks found.\n")
    q()
}

matrix_data <- read.table(matrix_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(matrix_data) < 2) q()
file_paths <- matrix_data[, 2]
sample_names <- matrix_data[, 1]
mean_matrix_list <- list()
for (i in seq_along(file_paths)) {
  f <- file_paths[i]
  if (!file.exists(f)) next




      

      mat <- read.table(
        f, header = FALSE, stringsAsFactors = FALSE, 
        skip = 1, fill = TRUE, sep="\t"
      )
      

      meta_line <- readLines(f, n = 1)
      meta_json <- sub("^@", "", meta_line)
      meta_info <- tryCatch(jsonlite::fromJSON(meta_json), error = function(e) NULL)
      
      if (is.null(meta_info) || is.null(meta_info$sample_labels) || is.null(meta_info$sample_boundaries)) {
          cat("Warning: Could not parse deeptools metadata from", f, "\n")
          next
      }
      
      sample_labels <- meta_info$sample_labels
      sample_boundaries <- meta_info$sample_boundaries
      
      if (ncol(mat) > 6) {
          bed_info <- mat[, 1:6, drop=FALSE]
          colnames(bed_info) <- c("V1", "V2", "V3", "V4", "V5", "V6")
          

          signal_mat <- mat[, 7:ncol(mat), drop=FALSE]
          signal_mat[] <- lapply(signal_mat, function(x) as.numeric(as.character(x)))
          signal_mat[is.na(signal_mat)] <- 0
          


          mean_mat <- data.frame(matrix(ncol = length(sample_labels), nrow = nrow(signal_mat)))
          colnames(mean_mat) <- sample_labels
          
          for (j in seq_along(sample_labels)) {
              start_col <- sample_boundaries[j] + 1
              end_col <- sample_boundaries[j + 1]
              
              if (start_col > end_col) {
                  mean_mat[, j] <- 0
              } else {
                  sample_cols <- signal_mat[, start_col:end_col, drop = FALSE]
                  mean_mat[, j] <- rowMeans(sample_cols, na.rm = TRUE)
              }
          }
          
          mean_matrix_list[[sample_names[i]]] <- cbind(bed_info, mean_mat)
      }
}

if (length(mean_matrix_list) == 0) {
    cat("Error: No matrix loaded successfully.\n")
    q()
}



ChIP_signal <- mean_matrix_list



RNAseq_gene_symbols <- RNAseq_gene$SYMBOL

chip_gene_symbols_list <- lapply(chip_list, function(df) df$SYMBOL)
chip_all_symbols <- unique(unlist(chip_gene_symbols_list))

gene_intersect <- intersect(RNAseq_gene_symbols, chip_all_symbols)
if (length(gene_intersect) == 0) {
    cat("Error: No common genes (SYMBOL) found between RNAseq and ChIPseq.\n")
    q()
}

RNAseq_gene$gene_id_simple <- RNAseq_gene$SYMBOL
RNAseq_gene <- RNAseq_gene[RNAseq_gene$gene_id_simple %in% gene_intersect, ]

chip_list <- lapply(chip_list, function(df) {
  df$gene_id_simple <- df$SYMBOL
  df[df$gene_id_simple %in% gene_intersect, ]
})

ChIP_signal <- lapply(ChIP_signal, function(df) {


    df
  })





rna_meta <- read.csv("config/RNAseq_metadata.csv", stringsAsFactors = FALSE)
group_map <- setNames(rna_meta$group, rna_meta$sample_name)
if (is.null(group_map) || length(group_map) == 0) {
    cat("group_map not found in RNAseq_metadata.csv\n")
    q()
}

get_group_means <- function(df, group_name, input_pattern = "Input") {





  

  input_cols <- grep("Input", colnames(df), ignore.case = TRUE, value = TRUE)
  ip_cols <- setdiff(colnames(df), c(input_cols, "V1","V2","V3","V4","V5","V6","V7"))
  
  if (length(ip_cols) == 0 || length(input_cols) == 0) {
      cat("Missing IP or Input columns in", group_name, "\n")
      return(NULL)
  }
  
  IP_mean <- rowMeans(df[, ip_cols, drop = FALSE], na.rm=TRUE)
  Input_mean <- rowMeans(df[, input_cols, drop = FALSE], na.rm=TRUE)
  



  

  promoter_df <- chip_list[[group_name]]
  if (is.null(promoter_df)) {
      cat("Cannot find promoter mapping for group:", group_name, "\n")
      return(NULL)
  }
  




  promoter_df$matrix_id <- paste0(promoter_df$seqnames, ":", promoter_df$start - 1, "-", promoter_df$end)
  mapping <- promoter_df[, c("matrix_id", "SYMBOL")]
  
  tmp <- data.frame(
    peak_id = df$V4,
    IP_mean = IP_mean, 
    Input_mean = Input_mean,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  

  tmp <- merge(tmp, mapping, by.x = "peak_id", by.y = "matrix_id", all.x = TRUE)

  tmp <- tmp[!is.na(tmp$SYMBOL), ]
  

  tmp_agg <- aggregate(cbind(IP_mean, Input_mean) ~ SYMBOL, data = tmp, FUN = mean)
  
  colnames(tmp_agg) <- c("gene", paste0(group_name, c("_IP_mean", "_Input_mean")))
  return(tmp_agg)
}

valid_groups <- names(mean_matrix_list)
if (length(valid_groups) == 0) q()
group_results <- list()
for (grp in valid_groups) {
    res <- get_group_means(ChIP_signal[[grp]], grp)
    if (!is.null(res)) {
        group_results[[grp]] <- res
    }
}
group_results <- Filter(Negate(is.null), group_results)
if (length(group_results) < 1) {
    cat("Error: No valid group results to merge.\n")
    q()
}
Result <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), group_results)

P_group_name <- names(group_map)[group_map == "P"]
T_group_name <- names(group_map)[group_map == "T"]
if (length(P_group_name) == 0 || length(T_group_name) == 0) {
    cat("Warning: Could not find P or T groups in yaml. Proceeding without correlation plots.\n")
    q()
}

P_group_name <- P_group_name[1]
T_group_name <- T_group_name[1]

RNAseq_gene$gene_id_simple <- sub("\\..*", "", RNAseq_gene$gene_id)


if (!(P_group_name %in% colnames(RNAseq_gene))) {
    # If the exact sample names are not there, try finding any columns that match P or T
    # Since we need a log2FC from DESeq2, let's just use log2FoldChange directly if present
    if ("log2FoldChange" %in% colnames(RNAseq_gene)) {
        RNAseq_gene$log2FoldChange_calc <- RNAseq_gene$log2FoldChange
    } else {
        cat("log2FoldChange not found in RNAseq output.\n")
        q()
    }
} else {
    RNAseq_gene$log2FoldChange_calc <- log2((RNAseq_gene[[P_group_name]] + 1) / (RNAseq_gene[[T_group_name]] + 1))
}

Result_all <- merge(Result, RNAseq_gene, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)
write.csv(Result_all, file.path(output_folder, "Promoter_gene_all.csv"), row.names = FALSE)

for (grp in valid_groups) {
  P_ip_col <- paste0(grp, "_IP_mean")
  T_ip_col <- paste0(grp, "_Input_mean")
  
  if (!P_ip_col %in% colnames(Result) || !T_ip_col %in% colnames(Result)) next
  
  grp_Result <- Result[!is.na(Result[[T_ip_col]]) & Result[[T_ip_col]] > 0, ]
  if (nrow(grp_Result) == 0) next
  
  grp_Result$log2FC <- log2(grp_Result[[P_ip_col]] / grp_Result[[T_ip_col]])
  
  grp_Result <- grp_Result %>%
    left_join(RNAseq_gene, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(log2FoldChange) & !is.na(padj))
  
  if (nrow(grp_Result) == 0) next
  
  RNAseq_up   <- grp_Result %>% filter(padj < 0.05 & log2FoldChange > 1) %>% pull(gene)  
  RNAseq_down <- grp_Result %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(gene)
  ChIPseq_up   <- grp_Result %>% filter(log2FC > 0.5) %>% pull(gene)
  ChIPseq_down <- grp_Result %>% filter(log2FC < -0.5) %>% pull(gene)
  
  futile.logger::flog.threshold(futile.logger::ERROR)
  
  venn_up <- venn.diagram(
    x = list(
      RNAseq  = RNAseq_up,
      ChIPseq = ChIPseq_up
    ),
    fill = c("#fc8d62", "#66c2a5"),
    alpha = c(0.5, 0.5),
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.8,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    scaled = FALSE,
    main = paste(grp, "Common Up-regulated Genes in Promoter Region: ", length(intersect(RNAseq_up, ChIPseq_up))),
    main.cex = 1.5,
    main.fontface = "bold",
    main.fontfamily = "sans",
    margin = 0.1,
    filename = NULL
  )
  venn_pdf_up <- file.path(outpic, paste0(grp, "_common_up_genes_venn_diagram.pdf"))
  pdf(file = venn_pdf_up, width = 9, height = 6)
  grid.draw(venn_up)
  dev.off()
  
  venn_down <- venn.diagram(
    x = list(
      RNAseq = RNAseq_down,
      ChIPseq = ChIPseq_down
    ),
    fill = c("#66c2a5", "#fc8d62"),
    alpha = c(0.5, 0.5),
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.8,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    scaled = FALSE,
    main = paste(grp, "Common Down-regulated Genes: ", length(intersect(RNAseq_down, ChIPseq_down))),
    main.cex = 1.5,
    main.fontface = "bold",
    main.fontfamily = "sans",
    margin = 0.1,
    filename = NULL
  )
  venn_pdf_down <- file.path(outpic, paste0(grp, "_common_down_genes_venn_diagram.pdf"))
  pdf(file = venn_pdf_down, width = 9, height = 6)
  grid.draw(venn_down)
  dev.off()
  
  # Volcano Plot
  plot_df <- grp_Result %>%
    filter(!is.na(padj) & padj > 0 & -log10(padj) <= 200) %>%
    mutate(group = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )) %>%
    mutate(group = factor(group, levels = c("Up", "NS", "Down")))
  
  topn <- 10
  logfc <- 1
  padj_cutoff <- 0.05
  
  plot_df_nonzero <- plot_df %>% filter(padj != 0)
  top_up <- plot_df_nonzero %>% filter(group == "Up") %>% arrange(padj) %>% slice_head(n = topn)
  top_down <- plot_df_nonzero %>% filter(group == "Down") %>% arrange(padj) %>% slice_head(n = topn)
  top_genes <- bind_rows(top_up, top_down)
  
  volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = log2FoldChange, size = -log10(padj)), alpha = 0.7) +
    geom_point(
      data = top_genes,
      aes(color = log2FoldChange, size = -log10(padj)),
      alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3, box.padding = 0.5, point.padding = 0.8,
      segment.color = "black", show.legend = FALSE, max.overlaps = Inf
    ) +
    scale_color_gradientn(
      colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
      values = scales::rescale(c(min(plot_df$log2FoldChange, na.rm = T), -1, 0, 1, max(plot_df$log2FoldChange, na.rm = T)))
    ) +
    scale_size(range = c(1, 7)) +
    geom_vline(xintercept = c(-logfc, logfc), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    xlim(c(-max(abs(plot_df$log2FoldChange), na.rm = T), max(plot_df$log2FoldChange, na.rm = T))) +
    ylim(c(0, max(-log10(plot_df$padj), na.rm = T) + 2)) +
    labs(
      x = "log2 (fold change)",
      y = "-log10 (padj)",
      title = paste(grp, "Volcano Plot (Common Genes)"),
      subtitle = paste("Top", topn, "up/down regulated genes")
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 13, color = "black"),
      axis.title = element_text(size = 15),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  volcano_pdf <- file.path(outpic, paste0(grp, "_common_genes_volcano_plot.pdf"))
  pdf(file = volcano_pdf, width = 10, height = 7.5)
  print(volcano_plot)
  dev.off()
  
  result_df_same_direction <- grp_Result[
    (grp_Result$log2FoldChange > 0 & grp_Result$log2FC > 0) |
      (grp_Result$log2FoldChange < 0 & grp_Result$log2FC < 0),
    ,]
  result_df_same_direction <- result_df_same_direction[ 
    is.finite(result_df_same_direction$log2FoldChange) & 
      is.finite(result_df_same_direction$log2FC), 
    ,]
  result_df_same_direction$color_group <- "Non-significant"
  result_df_same_direction$color_group[
    result_df_same_direction$log2FoldChange > 1 & result_df_same_direction$log2FC > 0.5
  ] <- "Up"
  result_df_same_direction$color_group[
    result_df_same_direction$log2FoldChange < -1 & result_df_same_direction$log2FC < -0.5
  ] <- "Down"
  color_map <- c("Up" = "#D55E00", "Down" = "#0072B2", "Non-significant" = "grey70")
  
  if (nrow(result_df_same_direction) > 1) {
    fit <- lm(log2FC ~ 0 + log2FoldChange, data=result_df_same_direction)
    k <- coef(fit)[1]
    cor_value <- cor(result_df_same_direction$log2FoldChange, result_df_same_direction$log2FC, method = "pearson")
    cor_label <- sprintf("Pearson~italic(R) == %.2f", cor_value)
    
    cor_pdf <- file.path(outpic, paste0(grp, "_RNA_ChIP_col_plot_same_direction.pdf"))
    pdf(file = cor_pdf, width=5, height=4)
    print(ggplot(result_df_same_direction, aes(x=log2FoldChange, y=log2FC, color=color_group)) +
      geom_point(alpha=0.7, size=2) +
      scale_color_manual(values=color_map) +
      geom_abline(intercept=0, slope=k, color="black", linewidth=1) +
      geom_hline(yintercept=0, linetype="dashed", color="grey") +
      geom_vline(xintercept=0, linetype="dashed", color="grey") +
      theme_classic() +
      labs(x="RNA-seq log2FC", y="ChIP-seq log2FC", title=paste(grp, "Same Direction (Up/Down)")) +
      annotate("text", 
               x = min(result_df_same_direction$log2FoldChange, na.rm=TRUE), 
               y = max(result_df_same_direction$log2FC, na.rm=TRUE), 
               label = cor_label, parse = TRUE, hjust = 0, vjust = 1, size = 5) +
      theme(legend.position="none",
            plot.title = element_text(hjust = 0.5)))
    dev.off()
  }
}


    plot_df <- Result_all %>%
      filter(!is.na(padj) & padj > 0 & -log10(padj) <= 200) %>%
      mutate(group = case_when(
        padj < 0.05 & log2FoldChange_calc > 1  ~ "Up",
        padj < 0.05 & log2FoldChange_calc < -1 ~ "Down",
        TRUE ~ "NS"
      )) %>%
      mutate(group = factor(group, levels = c("Up", "NS", "Down")))
    
    topn <- 10
    plot_df_nonzero <- plot_df %>% filter(padj != 0)
    top_up <- plot_df_nonzero %>% filter(group == "Up") %>% arrange(padj) %>% slice_head(n = topn)
    top_down <- plot_df_nonzero %>% filter(group == "Down") %>% arrange(padj) %>% slice_head(n = topn)
    top_genes <- bind_rows(top_up, top_down)
    
    volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange_calc, y = -log10(padj))) +
      geom_point(aes(color = log2FoldChange_calc, size = -log10(padj)), alpha = 0.7) +
      geom_point(
        data = top_genes,
        aes(color = log2FoldChange_calc, size = -log10(padj)),
        alpha = 0.95,
        show.legend = FALSE
      ) +
      geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = 3, box.padding = 0.5, point.padding = 0.8,
        segment.color = "black", show.legend = FALSE, max.overlaps = Inf
      ) +
      scale_color_gradientn(
        colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
        values = scales::rescale(c(min(plot_df$log2FoldChange_calc, na.rm = T), -1, 0, 1, max(plot_df$log2FoldChange_calc, na.rm = T)))
      ) +
      scale_size(range = c(1, 7)) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      xlim(c(-max(abs(plot_df$log2FoldChange_calc), na.rm = T), max(plot_df$log2FoldChange_calc, na.rm = T))) +
      ylim(c(0, max(-log10(plot_df$padj), na.rm = T) + 2)) +
      labs(
        x = "log2 (fold change)",
        y = "-log10 (padj)",
        title = "Volcano Plot (Common Genes)",
        subtitle = paste("Top", topn, "up/down regulated genes")
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
    ggsave(file.path(outpic, "common_genes_volcano_plot_calc.pdf"), plot = volcano_plot, width = 10, height = 7.5)

cat("Integrated Promoter Analysis Completed.\n")

for (grp in valid_groups) {
  P_ip_col <- paste0(grp, "_IP_mean")
  T_ip_col <- paste0(grp, "_Input_mean")
  
  if (!P_ip_col %in% colnames(Result) || !T_ip_col %in% colnames(Result)) next
  
  grp_Result <- Result[!is.na(Result[[T_ip_col]]) & Result[[T_ip_col]] > 0, ]
  if (nrow(grp_Result) == 0) next
  
  grp_Result$log2FC <- log2(grp_Result[[P_ip_col]] / grp_Result[[T_ip_col]])
  
  grp_Result <- grp_Result %>%
    left_join(RNAseq_gene, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(log2FoldChange) & !is.na(padj))
  
  if (nrow(grp_Result) == 0) next
  
  RNAseq_up   <- grp_Result %>% filter(padj < 0.05 & log2FoldChange > 1) %>% pull(gene)
  RNAseq_down <- grp_Result %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(gene)
  ChIPseq_up   <- grp_Result %>% filter(log2FC > 0.5) %>% pull(gene)
  ChIPseq_down <- grp_Result %>% filter(log2FC < -0.5) %>% pull(gene)
  
  common_up <- intersect(RNAseq_up, ChIPseq_up)
  common_down <- intersect(RNAseq_down, ChIPseq_down)
  
  peak_data <- chip_list[[grp]]
  if (is.null(peak_data)) next
  
  # Write UP peaks
  if (length(common_up) > 0) {
    bed_up <- peak_data %>%
      filter(SYMBOL %in% common_up) %>%
      distinct(seqnames, start, end, SYMBOL, .keep_all = TRUE) %>%
      dplyr::select(seqnames, start, end, SYMBOL)
    colnames(bed_up) <- c("chrom", "chromStart", "chromEnd", "name")
    write.table(bed_up,
                file = file.path(output_folder, paste0(grp, "_common_up_peaks.bed")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Write DOWN peaks
  if (length(common_down) > 0) {
    bed_down <- peak_data %>%
      filter(SYMBOL %in% common_down) %>%
      distinct(seqnames, start, end, SYMBOL, .keep_all = TRUE) %>%
      dplyr::select(seqnames, start, end, SYMBOL)
    colnames(bed_down) <- c("chrom", "chromStart", "chromEnd", "name")
    write.table(bed_down,
                file = file.path(output_folder, paste0(grp, "_common_down_peaks.bed")),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # Write ALL peaks
  bed_all <- peak_data %>%
    distinct(seqnames, start, end, SYMBOL, .keep_all = TRUE) %>%
    dplyr::select(seqnames, start, end, SYMBOL)
  colnames(bed_all) <- c("chrom", "chromStart", "chromEnd", "name")
  write.table(bed_all,
              file = file.path(output_folder, paste0(grp, "_all_peaks.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

