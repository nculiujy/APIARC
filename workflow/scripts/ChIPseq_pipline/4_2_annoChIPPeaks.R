#!/usr/bin/env Rscript

required_packages <- c(
  "optparse", "yaml", "ChIPseeker", "ggplot2", "GenomicFeatures",
  "GenomicRanges", "rtracklayer", "RColorBrewer"
)
missing_packages <- required_packages[!vapply(
  required_packages,
  function(pkg) requireNamespace(pkg, quietly = TRUE),
  logical(1)
)]
if (length(missing_packages) > 0) {
  stop(
    "missing R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them with the pinned environment: ",
    "conda env create -f workflow/envs/ChIPseq_4_peak_result.yaml"
  )
}

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(ChIPseeker)
  library(ggplot2)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(graphics)
  library(RColorBrewer)
})


option_list <- list(
  make_option(c("-p", "--peakfile"), type = "character", help = "Path to the input Peak file (.narrowPeak or .bed)"),
  make_option(c("-o", "--output"), type = "character", help = "Path to the output directory", default = NULL),
  make_option(c("-n", "--name"), type = "character", help = "Prefix name for output files", default = "PeakAnno"),
  make_option(c("-c", "--config"), type = "character", help = "Path to the config.yaml file", default = "config/config.yaml"),
  make_option("--assembly", type = "character", help = "Reference assembly, e.g. GRCh37", default = NULL),
  make_option("--txdb", type = "character", help = "TxDb package to annotate peaks with", default = NULL),
  make_option("--orgdb", type = "character", help = "OrgDb package used for gene symbols", default = NULL)
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (is.null(opt$peakfile)) {
  stop("--peakfile must be provided.")
}


if (is.null(opt$output)) {
  opt$output <- getwd()
  cat("No output path provided, using current working directory:", opt$output, "\n")
}

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

if (!file.exists(opt$peakfile)) {
  stop(
    "peak file not found: ", opt$peakfile,
    " -- MACS2 did not produce a peak file for this sample; refusing to write a placeholder result"
  )
}

if (!file.exists(opt$config)) {
  stop("The config file does not exist.")
}


config <- yaml.load_file(opt$config)


species <- "mm"
if (!is.null(config$species)) {
  species <- config$species
}

cat("Detected species from config:", species, "\n")

# The assembly, not the species, decides which TxDb is correct: annotating
# hg38 peaks with an hg19 TxDb mislabels promoters without failing.
assembly <- NULL
if (!is.null(opt$assembly) && nzchar(opt$assembly)) {
  assembly <- opt$assembly
} else if (!is.null(config$annotation$assembly)) {
  assembly <- config$annotation$assembly
}

ANNOTATION_PAIRS <- list(
  GRCh37 = c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"),
  hg19 = c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"),
  GRCh38 = c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
  hg38 = c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
  GRCm38 = c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
  mm10 = c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
  TAIR = c("TxDb.Athaliana.Biomart.plantsmart28", "org.At.tair.db")
)

SPECIES_PAIRS <- list(
  mm = c("TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"),
  homo = c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"),
  TAIR = c("TxDb.Athaliana.Biomart.plantsmart28", "org.At.tair.db")
)

pair <- NULL
if (!is.null(assembly)) {
  pair <- ANNOTATION_PAIRS[[assembly]]
}
if (is.null(pair)) {
  pair <- SPECIES_PAIRS[[species]]
}
if (is.null(pair)) {
  stop("Unsupported assembly/species: ", ifelse(is.null(assembly), species, assembly))
}

txdb_name <- pair[1]
orgdb_name <- pair[2]
if (!is.null(opt$txdb) && nzchar(opt$txdb)) {
  txdb_name <- opt$txdb
}
if (!is.null(opt$orgdb) && nzchar(opt$orgdb)) {
  orgdb_name <- opt$orgdb
}

cat("Using assembly:", ifelse(is.null(assembly), "(from species)", assembly), "\n")
cat("Using TxDb:", txdb_name, "\n")
cat("Using OrgDb:", orgdb_name, "\n")


if (!requireNamespace(txdb_name, quietly = TRUE)) {
  stop(
    "missing annotation package ", txdb_name,
    "; it must be part of the pinned environment, not installed at run time"
  )
}
if (!requireNamespace(orgdb_name, quietly = TRUE)) {
  stop(
    "missing annotation package ", orgdb_name,
    "; it must be part of the pinned environment, not installed at run time"
  )
}
suppressPackageStartupMessages({
  library(txdb_name, character.only = TRUE)
  library(orgdb_name, character.only = TRUE)
})

txdb <- get(txdb_name)
annoDb_name <- orgdb_name


options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)


if (file.info(opt$peakfile)$size == 0) {
  stop(
    "empty peak file: ", opt$peakfile,
    " -- MACS2 found no region under the configured cutoff. Lower the cutoff or report the sample as failing, but do not emit an empty annotation"
  )
}


peak <- tryCatch({
  readPeakFile(opt$peakfile)
}, error = function(e) {
  stop(
    "cannot read peak file ", opt$peakfile, ": ", conditionMessage(e),
    " -- the file is malformed; fix the upstream step instead of writing a placeholder"
  )
})


peakAnno <- annotatePeak(
  peak,
  tssRegion = c(-2000, 2000),
  TxDb = txdb,
  annoDb = annoDb_name
)


annotationData <- as.data.frame(peakAnno)


total_peaks <- nrow(annotationData)
cat("Total peaks:", total_peaks, "\n")


annotationData$annotation_category <- ifelse(grepl("Promoter", annotationData$annotation), "Promoter",
                                    ifelse(grepl("Intron", annotationData$annotation), "Intron",
                                           ifelse(grepl("Exon", annotationData$annotation), "Exon",
                                                  ifelse(grepl("3' UTR", annotationData$annotation), "3' UTR",
                                                         ifelse(grepl("Downstream", annotationData$annotation) | grepl("Distal Intergenic", annotationData$annotation), "Intergenic", annotationData$annotation)))))


annotationCounts <- table(annotationData$annotation_category)
annotationCounts <- as.data.frame(annotationCounts)
colnames(annotationCounts) <- c("Category", "Count")


annotationCounts <- annotationCounts[annotationCounts$Category %in% c("Promoter", "Intergenic", "Exon", "Intron", "3' UTR"), ]


if(sum(annotationCounts$Count) > 0) {
  annotationCounts$Percentage <- annotationCounts$Count / sum(annotationCounts$Count) * 100
} else {
  annotationCounts$Percentage <- 0
}


pdf_file <- file.path(opt$output, paste0(opt$name, "_peak_pie.pdf"))
pdf(pdf_file, width = 7.5, height = 7.5)

if (nrow(annotationCounts) > 0 && sum(annotationCounts$Count) > 0) {

  labels <- annotationCounts$Category
  values <- annotationCounts$Count
  percentages <- paste(format(round(annotationCounts$Percentage, 2), nsmall = 2), "%", sep = "")
  label_text <- paste(labels, "\n(", values, ", ", percentages, ")", sep = "")
  colors <- brewer.pal(max(3, length(annotationCounts$Category)), "Set3")[1:length(annotationCounts$Category)]

  pie(values, labels = label_text, col = colors, border = "white", main = paste(opt$name, "- Total peaks:", total_peaks))
} else {
  plot.new()
  text(0.5, 0.5, "No peaks matched required categories")
}
dev.off()


peak_anno_file <- file.path(opt$output, paste0(opt$name, "_peak_anno.csv"))
write.csv(annotationData, peak_anno_file, row.names = FALSE)


promoter_peaks <- annotationData[grepl("Promoter", annotationData$annotation_category), ]
promoter_file <- file.path(opt$output, paste0(opt$name, "_promoter_peaks.csv"))
write.csv(promoter_peaks, promoter_file, row.names = FALSE)

cat("Annotation finished successfully. Output saved to:", opt$output, "\n")