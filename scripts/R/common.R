# Shared functions of the R steps of the CUT&RUN pipeline: reading tables,
# DESeq2 helpers and plots in a Nature-journal style. Sourced by the other R
# scripts.

# Okabe-Ito colors: distinguishable with color blindness
kConditionPalette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
                       "#D55E00", "#F0E442", "#000000")
kDirectionColors <- c(Up = "#D55E00", Down = "#0072B2", NS = "grey80")
kFontSize <- 7         # pt, as asked by Nature for figure text
kSingleColumn <- 89    # mm, Nature single-column width
kDoubleColumn <- 183   # mm, Nature double-column width

#' Reads a tab-separated table, skipping lines that start with "#".
#'
#' @param file Path of the table, e.g. the sample list.
#' @return A data frame with character columns.
ReadTable <- function(file) {
  lines <- sub("\r$", "", readLines(file, warn = FALSE))
  lines <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
  table <- data.table::fread(text = paste(lines, collapse = "\n"),
                             sep = "\t", colClasses = "character")
  return(as.data.frame(table))
}

#' Tells if a file name is given and the file exists.
#'
#' @param path File path, may be NULL or "".
#' @return TRUE or FALSE.
IsFile <- function(path) {
  return(!is.null(path) && nzchar(path) && file.exists(path))
}

#' Writes a data frame as a tab-separated text file.
#'
#' @param table Data frame.
#' @param file Output file.
WriteTable <- function(table, file) {
  data.table::fwrite(table, file, sep = "\t", na = "NA", quote = FALSE)
  return(invisible(file))
}

#' Gives each condition a color (Okabe-Ito palette).
#'
#' @param conditions Vector of conditions, in plot order.
#' @return Named vector of colors.
ConditionColors <- function(conditions) {
  conditions <- unique(as.character(conditions))
  colors <- rep(kConditionPalette, length.out = length(conditions))
  return(stats::setNames(colors, conditions))
}

#' ggplot2 theme in the style of Nature figures: 7 pt sans-serif text, thin
#' black axes, no grid lines or background.
#'
#' @param base_size Font size in pt.
#' @return A ggplot2 theme.
ThemeNature <- function(base_size = kFontSize) {
  line_width <- 0.25  # mm
  return(
    ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
      ggplot2::theme(
        axis.line = ggplot2::element_line(linewidth = line_width),
        axis.ticks = ggplot2::element_line(linewidth = line_width),
        axis.ticks.length = ggplot2::unit(1, "mm"),
        axis.text = ggplot2::element_text(size = base_size, color = "black"),
        axis.title = ggplot2::element_text(size = base_size),
        legend.key.size = ggplot2::unit(3, "mm"),
        legend.title = ggplot2::element_text(size = base_size),
        legend.text = ggplot2::element_text(size = base_size),
        plot.title = ggplot2::element_text(size = base_size, face = "bold"),
        plot.tag = ggplot2::element_text(size = 8, face = "bold"),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = base_size, face = "bold")
      )
  )
}

#' Saves a ggplot as a vector PDF; sizes in mm (journal column widths).
#'
#' @param plot ggplot or patchwork object.
#' @param file Output PDF.
#' @param width,height Size in mm.
SaveFigure <- function(plot, file, width = kSingleColumn, height = 70) {
  ggplot2::ggsave(file, plot, width = width, height = height, units = "mm",
                  device = grDevices::pdf, useDingbats = FALSE)
  return(invisible(file))
}

#' Variance-stabilized counts (log2 scale) for PCA and heatmaps.
#'
#' Uses DESeq2 vst() for 1000 or more peaks, otherwise the slower
#' varianceStabilizingTransformation(), and log2(normalized + 1) if both fail.
#'
#' @param dds DESeqDataSet with size factors.
#' @return Matrix, peaks x samples.
VstMatrix <- function(dds) {
  transform <- if (nrow(dds) >= 1000) {
    DESeq2::vst
  } else {
    DESeq2::varianceStabilizingTransformation
  }
  vst_counts <- tryCatch(
    SummarizedExperiment::assay(transform(dds, blind = TRUE)),
    error = function(e) log2(DESeq2::counts(dds, normalized = TRUE) + 1)
  )
  return(vst_counts)
}

#' Reads a deepTools computeMatrix output (signal around each region).
#'
#' @param file computeMatrix output (.gz).
#' @return List: regions (names), samples, bin_centers (bp from the center),
#'   values (matrix, regions x (samples x bins)).
ReadSignalMatrix <- function(file) {
  con <- gzfile(file)
  header <- jsonlite::fromJSON(sub("^@", "", readLines(con, n = 1)))
  close(con)
  table <- data.table::fread(
    cmd = paste("gzip -dc", shQuote(file), "| tail -n +2"), header = FALSE
  )
  window <- header$upstream[1] + header$downstream[1]
  n_bins <- window / header[["bin size"]][1]
  bin_centers <- seq(-header$upstream[1], header$downstream[1],
                     length.out = n_bins + 1)
  bin_centers <- (utils::head(bin_centers, -1) + bin_centers[-1]) / 2
  return(list(regions = table[[4]], samples = header$sample_labels,
              bin_centers = bin_centers,
              values = as.matrix(table[, -(1:6)])))
}

#' Mean signal of the replicates of each condition.
#'
#' @param signal Output of ReadSignalMatrix().
#' @param sample_condition Named vector: condition of each sample.
#' @return Named list (one per condition) of matrices, regions x bins.
ConditionSignal <- function(signal, sample_condition) {
  n_bins <- length(signal$bin_centers)
  sample_columns <- split(seq_len(ncol(signal$values)),
                          rep(signal$samples, each = n_bins))
  conditions <- unique(sample_condition[signal$samples])
  means <- lapply(conditions, function(cond) {
    samples <- signal$samples[sample_condition[signal$samples] == cond]
    blocks <- lapply(samples, function(s) signal$values[, sample_columns[[s]]])
    mean_signal <- Reduce(`+`, blocks) / length(blocks)
    rownames(mean_signal) <- signal$regions
    return(mean_signal)
  })
  return(stats::setNames(means, conditions))
}

#' PCA of samples, using the most variable peaks.
#'
#' @param counts Matrix, peaks x samples (e.g. VST counts).
#' @param condition Condition of each sample (column).
#' @param title Plot title.
#' @param n_top Number of most variable peaks to use.
#' @return ggplot object.
PcaPlot <- function(counts, condition, title, n_top = 5000) {
  variance <- apply(counts, 1, stats::var)
  top <- order(variance, decreasing = TRUE)[seq_len(min(n_top, nrow(counts)))]
  pca <- stats::prcomp(t(counts[top, , drop = FALSE]))
  percent <- round(summary(pca)$importance[2, ] * 100, 1)
  pcs <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    sample = colnames(counts),
                    condition = factor(condition, levels = unique(condition)))

  plot <- ggplot2::ggplot(pcs, ggplot2::aes(PC1, PC2, color = condition)) +
    ggplot2::geom_point(size = 1.5) +
    ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 2,
                             color = "grey30", segment.size = 0.2) +
    ggplot2::scale_color_manual(values = ConditionColors(condition)) +
    ggplot2::scale_x_continuous(n.breaks = 4) +
    ggplot2::labs(x = paste0("PC1 (", percent[1], "%)"),
                  y = paste0("PC2 (", percent[2], "%)"), color = NULL,
                  title = title) +
    ThemeNature()
  return(plot)
}

#' Correlation heatmap of samples (PDF).
#'
#' @param counts Matrix, peaks x samples.
#' @param condition Condition of each sample.
#' @param title Plot title.
#' @param file Output PDF.
PlotCorrelation <- function(counts, condition, title, file) {
  annotation <- data.frame(condition = condition, row.names = colnames(counts))
  blues <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  pheatmap::pheatmap(
    stats::cor(counts), annotation_col = annotation,
    annotation_row = annotation,
    annotation_colors = list(condition = ConditionColors(condition)),
    color = blues(100), border_color = NA, fontsize = kFontSize, main = title,
    filename = file, width = kDoubleColumn / 25.4 * 0.6,
    height = kDoubleColumn / 25.4 * 0.5
  )
  return(invisible(file))
}

#' MA plot of a comparison.
#'
#' @param result Data frame with baseMean, log2FoldChange and direction.
#' @param title Plot title.
#' @return ggplot object.
MaPlot <- function(result, title) {
  result <- result[!is.na(result$log2FoldChange), ]
  plot <- ggplot2::ggplot(
    result,
    ggplot2::aes(log10(baseMean + 1), log2FoldChange, color = direction)
  ) +
    ggplot2::geom_point(size = 0.4, stroke = 0) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey40") +
    ggplot2::scale_color_manual(values = kDirectionColors) +
    ggplot2::labs(x = "log10 mean normalized count", y = "log2 fold change",
                  color = NULL, title = title) +
    ThemeNature()
  return(plot)
}

#' Volcano plot; the top significant peaks are labeled with their gene.
#'
#' @param result Data frame with log2FoldChange, padj, direction, peakN and
#'   Gene.Name (optional).
#' @param title Plot title.
#' @param n_label Number of peaks to label.
#' @return ggplot object.
VolcanoPlot <- function(result, title, n_label = 15) {
  result <- result[!is.na(result$padj), ]
  result$label <- result$peakN
  if ("Gene.Name" %in% colnames(result)) {
    has_gene <- !is.na(result$Gene.Name)
    result$label[has_gene] <- result$Gene.Name[has_gene]
  }
  significant <- result[result$direction != "NS", ]
  to_label <- utils::head(significant[order(significant$padj), ], n_label)

  plot <- ggplot2::ggplot(
    result, ggplot2::aes(log2FoldChange, -log10(padj), color = direction)
  ) +
    ggplot2::geom_point(size = 0.4, stroke = 0) +
    ggrepel::geom_text_repel(data = to_label, ggplot2::aes(label = label),
                             size = 2, color = "black", segment.size = 0.2,
                             max.overlaps = 50) +
    ggplot2::scale_color_manual(values = kDirectionColors) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted P", color = NULL,
                  title = title) +
    ThemeNature()
  return(plot)
}

#' Heatmap (z-score of each peak) of significant peaks (PDF).
#'
#' @param counts Matrix, peaks x samples, rows in order of significance.
#' @param condition Condition of each sample.
#' @param title Plot title.
#' @param file Output PDF.
#' @param max_rows Maximum number of peaks shown.
PlotHeatmap <- function(counts, condition, title, file, max_rows = 2000) {
  counts <- counts[seq_len(min(nrow(counts), max_rows)), , drop = FALSE]
  counts <- counts[apply(counts, 1, stats::sd) > 0, , drop = FALSE]
  if (nrow(counts) < 2) {
    return(invisible(NULL))
  }
  annotation <- data.frame(condition = condition, row.names = colnames(counts))
  blue_red <- grDevices::colorRampPalette(
    rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )
  pheatmap::pheatmap(
    t(scale(t(counts))), annotation_col = annotation,
    annotation_colors = list(condition = ConditionColors(condition)),
    color = blue_red(100), border_color = NA, show_rownames = FALSE,
    fontsize = kFontSize, main = paste0(title, " (", nrow(counts), " peaks)"),
    filename = file, width = kSingleColumn / 25.4, height = 120 / 25.4
  )
  return(invisible(file))
}

#' Bar plot of the numbers of Up and Down peaks of each comparison.
#'
#' @param summary Data frame with columns comparison, Up and Down.
#' @param title Plot title.
#' @return ggplot object.
DiffSummaryPlot <- function(summary, title) {
  bars <- rbind(
    data.frame(comparison = summary$comparison, direction = "Up",
               count = summary$Up),
    data.frame(comparison = summary$comparison, direction = "Down",
               count = -summary$Down)
  )
  bars$comparison <- factor(bars$comparison, levels = summary$comparison)

  plot <- ggplot2::ggplot(bars,
                          ggplot2::aes(comparison, count, fill = direction)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25) +
    ggplot2::geom_text(data = bars[bars$count != 0, ],
                       ggplot2::aes(label = abs(count),
                                    vjust = ifelse(count > 0, -0.3, 1.2)),
                       size = 2) +
    ggplot2::scale_fill_manual(values = kDirectionColors) +
    ggplot2::labs(x = NULL, y = "Differential peaks", fill = NULL,
                  title = title) +
    ThemeNature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  return(plot)
}
