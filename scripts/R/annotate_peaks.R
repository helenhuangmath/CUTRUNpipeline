#!/usr/bin/env Rscript
#
# Step 6: annotation of peaks with ChIPseeker (and Homer, if given).
# Called by scripts/06_annotate.sh.
#
# Output:
#   <out>.anno.txt          one row per peak: Anno (Promoter, Intron, ...),
#                           Distance.to.TSS, nearest gene (geneId, Gene.Name)
#   <out>.anno_allgenes.txt every gene linked to a peak; type: 1_nearest,
#                           2_flank (within --flank bp) or 3_homer
#   <out>.anno_summary.txt, plot_anno_pie_<name>.pdf: % of peaks per feature

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

kFeatureLevels <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "TTS",
                    "Non-Coding", "Intergenic")

#' Annotates peaks with ChIPseeker.
#'
#' @param peaks Data frame: Chr, Start (0-based), End, peakN.
#' @param txdb_file TxDb sqlite file.
#' @param tss_region Promoter region around the TSS, e.g. c(-2000, 500).
#' @param flank Distance (bp) for flanking genes.
#' @return List: nearest (one row per peak) and genes (all linked genes).
AnnotateWithChipseeker <- function(peaks, txdb_file, tss_region, flank) {
  txdb <- AnnotationDbi::loadDb(txdb_file)
  ranges <- GenomicRanges::GRanges(
    peaks$Chr, IRanges::IRanges(peaks$Start + 1, peaks$End),
    peakN = peaks$peakN
  )
  anno <- as.data.frame(ChIPseeker::annotatePeak(
    ranges, TxDb = txdb, tssRegion = tss_region, flankDistance = flank,
    addFlankGeneInfo = TRUE, verbose = FALSE
  ))

  # "Promoter (<=1kb)" -> "Promoter", "Intron (ENSMUST..., ...)" -> "Intron"
  anno$Anno <- sub(" \\(.*", "", anno$annotation)
  anno$Anno[anno$Anno %in% c("Distal Intergenic", "Downstream")] <-
    "Intergenic"

  nearest <- data.frame(peakN = anno$peakN, Anno = anno$Anno,
                        Distance.to.TSS = anno$distanceToTSS,
                        geneId = anno$geneId)
  genes <- data.frame(peakN = anno$peakN, geneId = anno$geneId,
                      Distance.to.TSS = anno$distanceToTSS, type = "1_nearest")

  # Flanking genes: "geneA;geneB" in one row -> one row per gene
  flanking <- anno[!is.na(anno$flank_geneIds), ]
  if (nrow(flanking) > 0) {
    ids <- strsplit(flanking$flank_geneIds, ";")
    distances <- strsplit(as.character(flanking$flank_gene_distances), ";")
    genes <- rbind(genes, data.frame(
      peakN = rep(flanking$peakN, lengths(ids)), geneId = unlist(ids),
      Distance.to.TSS = as.numeric(unlist(distances)), type = "2_flank"
    ))
  }
  return(list(nearest = nearest, genes = genes))
}

#' Reads Homer annotatePeaks.pl output.
#'
#' @param file Homer output.
#' @param tss_region Promoter region around the TSS.
#' @return List: nearest (one row per peak) and genes.
ReadHomer <- function(file, tss_region) {
  homer <- as.data.frame(data.table::fread(file, quote = ""))
  colnames(homer)[1] <- "peakN"
  anno <- sub(" \\(.*", "", homer$Annotation)
  anno <- c("promoter-TSS" = "Promoter", "intron" = "Intron", "exon" = "Exon",
            "non-coding" = "Non-Coding", "TTS" = "TTS", "5' UTR" = "5' UTR",
            "3' UTR" = "3' UTR", "Intergenic" = "Intergenic")[anno]
  distance <- homer$`Distance to TSS`
  in_promoter <- distance > tss_region[1] & distance < tss_region[2]
  anno[!is.na(distance) & in_promoter] <- "Promoter"

  nearest <- data.frame(peakN = homer$peakN, Anno = unname(anno),
                        Distance.to.TSS = distance,
                        Gene.Name = homer$`Gene Name`)
  genes <- data.frame(peakN = homer$peakN, geneId = homer$`Nearest Ensembl`,
                      Distance.to.TSS = distance, type = "3_homer",
                      Gene.Name = homer$`Gene Name`)
  return(list(nearest = nearest, genes = genes))
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--peaks", help = "union peaks BED"),
  optparse::make_option("--txdb", default = ""),
  optparse::make_option("--genes", default = "", help = "gene table"),
  optparse::make_option("--homer", default = ""),
  optparse::make_option("--tss-up", dest = "tss_up", type = "integer",
                        default = 2000),
  optparse::make_option("--tss-down", dest = "tss_down", type = "integer",
                        default = 500),
  optparse::make_option("--flank", type = "integer", default = 10000),
  optparse::make_option("--out", help = "output prefix")
)))
tss_region <- c(-opts$tss_up, opts$tss_down)

peaks <- as.data.frame(data.table::fread(
  opts$peaks, header = FALSE, select = 1:4,
  col.names = c("Chr", "Start", "End", "peakN")
))
peaks$Width <- peaks$End - peaks$Start
gene_table <- NULL
if (IsFile(opts$genes)) {
  gene_table <- as.data.frame(data.table::fread(opts$genes))
}
GeneName <- function(ids) {
  if (is.null(gene_table)) return(rep(NA_character_, length(ids)))
  return(gene_table$Gene.Name[match(ids, gene_table$geneId)])
}

# 1. ChIPseeker and/or Homer
anno <- peaks
all_genes <- NULL
if (IsFile(opts$txdb)) {
  chipseeker <- AnnotateWithChipseeker(peaks, opts$txdb, tss_region,
                                       opts$flank)
  chipseeker$nearest$Gene.Name <- GeneName(chipseeker$nearest$geneId)
  anno <- merge(anno, chipseeker$nearest, by = "peakN", all.x = TRUE)
  chipseeker$genes$Gene.Name <- GeneName(chipseeker$genes$geneId)
  all_genes <- chipseeker$genes
}
if (IsFile(opts$homer)) {
  homer <- ReadHomer(opts$homer, tss_region)
  if (is.null(all_genes)) {
    anno <- merge(anno, homer$nearest, by = "peakN", all.x = TRUE)
  } else {
    homer_index <- match(anno$peakN, homer$nearest$peakN)
    anno$Homer.Gene <- homer$nearest$Gene.Name[homer_index]
  }
  all_genes <- rbind(all_genes, homer$genes)
}
anno <- anno[match(peaks$peakN, anno$peakN), ]

# 2. Output
WriteTable(anno, paste0(opts$out, ".anno.txt"))
if (!is.null(all_genes)) {
  # Keep each peak-gene pair once (nearest first)
  all_genes <- all_genes[order(all_genes$peakN, all_genes$type), ]
  duplicated_pair <- duplicated(all_genes[, c("peakN", "Gene.Name")]) &
    !is.na(all_genes$Gene.Name)
  all_genes <- all_genes[!duplicated_pair,
                         c("peakN", "Gene.Name", "geneId", "Distance.to.TSS",
                           "type")]
  WriteTable(all_genes, paste0(opts$out, ".anno_allgenes.txt"))
}

# 3. Percent of peaks in each genomic feature
if ("Anno" %in% colnames(anno) && any(!is.na(anno$Anno))) {
  features <- as.data.frame(table(factor(anno$Anno, levels = kFeatureLevels)))
  colnames(features) <- c("Anno", "n")
  features <- features[features$n > 0, ]
  features$percent <- round(features$n / sum(features$n) * 100, 1)
  features$label <- paste0(features$Anno, " (", features$percent, "%)")
  WriteTable(features, paste0(opts$out, ".anno_summary.txt"))

  plot <- ggplot2::ggplot(features,
                          ggplot2::aes(x = "", y = percent, fill = Anno)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_brewer(palette = "Paired", labels = features$label) +
    ggplot2::labs(title = paste0(basename(opts$out), " (", nrow(anno),
                                 " peaks)")) +
    ggplot2::theme_void(base_size = kFontSize)
  SaveFigure(plot, file.path(dirname(opts$out),
                             paste0("plot_anno_pie_", basename(opts$out),
                                    ".pdf")),
             height = 60)
}
