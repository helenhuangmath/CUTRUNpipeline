#!/usr/bin/env Rscript
#
# Step 9: Nature-style summary figure of one mark (183 mm wide, 7 pt text,
# panels a-h). Called by scripts/09_figures.sh.
#
#   a  peaks of each sample (MACS2 / FlexPeak)   e  differential peaks
#   b  FRiP of each sample                       f  volcano plots
#   c  genomic features of consensus peaks       g  average signal profiles
#   d  PCA of samples                            h  signal heatmap
#
# Output (Figures/<mark>/): Figure_<mark>_summary.pdf,
#   signal_profile_<mark>.pdf, signal_heatmap_<mark>.pdf

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

kFeatureLevels <- c("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "TTS",
                    "Non-Coding", "Intergenic")

#' a. Number of peaks of each sample, for each peak caller.
PeakNumberPlot <- function(workdir, meta, callers) {
  counts <- expand.grid(sample = meta$sample, caller = callers,
                        stringsAsFactors = FALSE)
  counts$peaks <- mapply(function(sample, caller) {
    file <- file.path(workdir, "Peak", paste0("peak_", sample, ".", caller,
                                              ".bed"))
    if (!file.exists(file)) return(NA_integer_)
    return(length(readLines(file)))
  }, counts$sample, counts$caller)
  counts$sample <- factor(counts$sample, levels = meta$sample)
  counts$caller <- c(macs2 = "MACS2", flexpeak = "FlexPeak")[counts$caller]

  plot <- ggplot2::ggplot(counts,
                          ggplot2::aes(sample, peaks, fill = caller)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.75), width = 0.7) +
    ggplot2::scale_fill_manual(values = c(MACS2 = "grey30",
                                          FlexPeak = "#E69F00")) +
    ggplot2::labs(x = NULL, y = "Peaks", fill = NULL) +
    ThemeNature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  return(plot)
}

#' b. FRiP of each sample, by condition (dots and mean).
FripPlot <- function(size_factors, conditions) {
  size_factors$condition <- factor(size_factors$condition, levels = conditions)
  plot <- ggplot2::ggplot(size_factors,
                          ggplot2::aes(condition, FRiP, color = condition)) +
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.5,
                          linewidth = 0.25, color = "black") +
    ggplot2::geom_point(size = 1.2,
                        position = ggplot2::position_jitter(0.1, seed = 1)) +
    ggplot2::scale_color_manual(values = ConditionColors(conditions)) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::labs(x = NULL, y = "FRiP") +
    ThemeNature() +
    ggplot2::theme(legend.position = "none")
  return(plot)
}

#' c. Genomic features of the consensus peaks of each condition (% stacked).
FeaturePlot <- function(consensus_dir, anno, conditions) {
  features <- do.call(rbind, lapply(conditions, function(cond) {
    ids <- readLines(file.path(consensus_dir,
                               paste0("peak_", cond, "_in_union.txt")))
    return(data.frame(condition = cond,
                      Anno = anno$Anno[match(ids, anno$peakN)]))
  }))
  features <- features[!is.na(features$Anno), ]
  features$Anno <- factor(features$Anno,
                          levels = intersect(kFeatureLevels, features$Anno))
  features$condition <- factor(features$condition, levels = conditions)

  plot <- ggplot2::ggplot(features, ggplot2::aes(condition, fill = Anno)) +
    ggplot2::geom_bar(position = "fill", width = 0.7) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::scale_fill_brewer(palette = "Paired") +
    ggplot2::labs(x = NULL, y = "Consensus peaks", fill = NULL) +
    ThemeNature()
  return(plot)
}

#' f. Volcano plots of all comparisons.
VolcanoPanel <- function(diff_dir, mark) {
  files <- list.files(diff_dir, pattern = paste0("^deseq2_All_", mark, "_"),
                      full.names = TRUE)
  if (length(files) == 0) return(NULL)
  results <- do.call(rbind, lapply(files, function(file) {
    table <- as.data.frame(data.table::fread(file))
    columns <- c("peakN", "Gene.Name", "log2FoldChange", "padj", "direction",
                 "comparison")
    return(table[, intersect(columns, colnames(table))])
  }))
  results$comparison <- gsub("_vs_", "\nvs ", results$comparison)
  plot <- VolcanoPlot(results, NULL, n_label = 0) +
    ggplot2::facet_wrap(~ comparison, nrow = 1) +
    ggplot2::scale_x_continuous(n.breaks = 3) +
    ggplot2::theme(legend.position = "none")
  return(plot)
}

#' Groups peaks by their direction in any comparison.
#'
#' @return Named vector: group of each peak (Up / Down / Not differential).
PeakGroups <- function(diff_dir, mark, peaks) {
  files <- list.files(diff_dir, pattern = paste0("^deseq2_Sig_", mark,
                                                 "_.*_vs_.*\\.txt$"),
                      full.names = TRUE)
  groups <- stats::setNames(rep("Not differential", length(peaks)), peaks)
  for (file in files) {
    table <- as.data.frame(data.table::fread(file))
    groups[table$peakN[table$direction == "Down"]] <- "Down"
  }
  for (file in files) {
    table <- as.data.frame(data.table::fread(file))
    groups[table$peakN[table$direction == "Up"]] <- "Up"
  }
  return(factor(groups, levels = c("Up", "Down", "Not differential")))
}

#' g. Average signal around peak centers, by condition and peak group.
ProfilePlot <- function(signal_by_condition, bin_centers, groups) {
  profiles <- do.call(rbind, lapply(names(signal_by_condition), function(cond) {
    values <- signal_by_condition[[cond]]
    by_group <- lapply(levels(groups), function(group) {
      in_group <- names(groups)[groups == group]
      if (length(in_group) == 0) return(NULL)
      return(data.frame(condition = cond, group = group,
                        position = bin_centers,
                        signal = colMeans(values[in_group, , drop = FALSE])))
    })
    return(do.call(rbind, by_group))
  }))
  profiles$condition <- factor(profiles$condition,
                               levels = names(signal_by_condition))
  profiles$group <- factor(profiles$group, levels = levels(groups))

  plot <- ggplot2::ggplot(profiles,
                          ggplot2::aes(position / 1000, signal,
                                       color = condition)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::facet_wrap(~ group, nrow = 1) +
    ggplot2::scale_color_manual(
      values = ConditionColors(names(signal_by_condition))
    ) +
    ggplot2::labs(x = "Distance to peak center (kb)", y = "Mean signal (CPM)",
                  color = NULL) +
    ThemeNature()
  return(plot)
}

#' h. Heatmap of the signal at differential peaks (all peaks if none).
HeatmapPlot <- function(signal_by_condition, bin_centers, groups,
                        max_peaks = 3000) {
  peaks <- names(groups)[groups != "Not differential"]
  if (length(peaks) == 0) peaks <- names(groups)
  total <- Reduce(`+`, lapply(signal_by_condition, function(m) {
    return(rowMeans(m[peaks, , drop = FALSE]))
  }))
  peaks <- utils::head(peaks[order(total, decreasing = TRUE)], max_peaks)

  conditions <- names(signal_by_condition)
  tiles <- do.call(rbind, lapply(conditions, function(cond) {
    values <- signal_by_condition[[cond]][peaks, , drop = FALSE]
    return(data.frame(condition = cond, peak = rep(peaks, ncol(values)),
                      group = rep(as.character(groups[peaks]), ncol(values)),
                      position = rep(bin_centers, each = nrow(values)),
                      signal = as.vector(values)))
  }))
  tiles$peak <- factor(tiles$peak, levels = rev(peaks))
  tiles$condition <- factor(tiles$condition, levels = conditions)
  cap <- stats::quantile(tiles$signal, 0.99)
  tiles$signal <- pmin(tiles$signal, cap)

  plot <- ggplot2::ggplot(tiles,
                          ggplot2::aes(position / 1000, peak, fill = signal)) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(group ~ condition, scales = "free_y",
                        space = "free_y") +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B") +
    ggplot2::labs(x = "Distance to peak center (kb)", y = NULL,
                  fill = "Signal\n(CPM)") +
    ThemeNature() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank())
  return(plot)
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--mark"),
  optparse::make_option("--workdir"),
  optparse::make_option("--sample-list", dest = "sample_list"),
  optparse::make_option("--peak-callers", dest = "peak_callers",
                        default = "macs2")
)))
mark <- opts$mark
workdir <- opts$workdir
out_dir <- file.path(workdir, "Figures", mark)
diff_dir <- file.path(workdir, "DiffPeak", mark)
count_dir <- file.path(workdir, "ReadCount", mark)

meta <- ReadTable(opts$sample_list)
meta <- meta[meta$mark == mark, ]
conditions <- unique(meta$condition)
sample_condition <- stats::setNames(meta$condition, meta$sample)
anno <- as.data.frame(data.table::fread(file.path(
  workdir, "Annotation", mark, paste0("peak_", mark, "_union.anno.txt")
)))
vst <- as.data.frame(data.table::fread(file.path(count_dir,
                                                 paste0("VST_", mark, ".txt"))))
size_factors <- as.data.frame(data.table::fread(
  file.path(count_dir, paste0("sizeFactors_", mark, ".txt"))
))

signal <- ReadSignalMatrix(file.path(out_dir,
                                     paste0("signal_", mark, ".matrix.gz")))
signal_by_condition <- ConditionSignal(signal, sample_condition)
groups <- PeakGroups(diff_dir, mark, signal$regions)

# Panels
panels <- list(
  a = PeakNumberPlot(workdir, meta, strsplit(opts$peak_callers, " ")[[1]]),
  b = FripPlot(size_factors, conditions),
  c = FeaturePlot(file.path(workdir, "Consensus", mark), anno, conditions),
  d = PcaPlot(as.matrix(vst[, meta$sample]), meta$condition, NULL)
)
summary_file <- file.path(diff_dir, paste0("deseq2_summary_", mark, ".txt"))
if (file.exists(summary_file)) {
  summary <- as.data.frame(data.table::fread(summary_file))
  summary$comparison <- gsub("_vs_", " vs ", summary$comparison)
  panels$e <- DiffSummaryPlot(summary, NULL)
}
panels$f <- VolcanoPanel(diff_dir, mark)
profile <- ProfilePlot(signal_by_condition, signal$bin_centers, groups)
heatmap <- HeatmapPlot(signal_by_condition, signal$bin_centers, groups)
panels <- panels[!vapply(panels, is.null, logical(1))]

# Layout: 3 small panels per row, then profiles and heatmap
top <- patchwork::wrap_plots(panels[intersect(c("a", "b", "c"), names(panels))],
                             nrow = 1, widths = c(1.4, 0.8, 1))
middle_widths <- c(d = 1, e = 0.7, f = 1.8)
middle_panels <- intersect(names(middle_widths), names(panels))
middle <- patchwork::wrap_plots(panels[middle_panels], nrow = 1,
                                widths = middle_widths[middle_panels])
figure <- (top / middle / profile / heatmap) +
  patchwork::plot_layout(heights = c(1, 1, 0.8, 1.6)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = paste(mark, "CUT&RUN summary"),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, face = "bold")
    )
  )
SaveFigure(figure, file.path(out_dir, paste0("Figure_", mark, "_summary.pdf")),
           width = kDoubleColumn, height = 240)
SaveFigure(profile, file.path(out_dir, paste0("signal_profile_", mark, ".pdf")),
           width = kDoubleColumn, height = 55)
SaveFigure(heatmap, file.path(out_dir, paste0("signal_heatmap_", mark, ".pdf")),
           width = kDoubleColumn, height = 120)
message(mark, ": summary figure in ", out_dir)
