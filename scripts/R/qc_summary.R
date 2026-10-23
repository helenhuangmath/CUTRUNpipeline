#!/usr/bin/env Rscript
#
# QC summary of all samples: read pairs, alignment rate, duplicate rate,
# number of peaks, FRiP and fragment sizes.
# Called by scripts/run_all_samples.sh. Output: QC/Summary/

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

#' Reads the alignment stats (step 2) of all samples.
#'
#' @param workdir Pipeline WORKDIR.
#' @param samples Sample IDs.
#' @return Data frame, one row per sample.
ReadAlignStats <- function(workdir, samples) {
  stats <- lapply(samples, function(sample) {
    file <- file.path(workdir, "BAM", paste0(sample, ".align_stats.txt"))
    if (!file.exists(file)) return(data.frame(sample = sample))
    return(as.data.frame(data.table::fread(file, colClasses = "character")))
  })
  return(as.data.frame(data.table::rbindlist(stats, fill = TRUE)))
}

#' Number of cleaned peaks (step 4) of each sample.
#'
#' @param workdir Pipeline WORKDIR.
#' @param samples Sample IDs.
#' @return Integer vector.
CountPeaks <- function(workdir, samples) {
  return(vapply(samples, function(sample) {
    file <- file.path(workdir, "Peak", paste0("peak_", sample, ".final.bed"))
    if (!file.exists(file)) return(NA_integer_)
    return(length(readLines(file)))
  }, integer(1)))
}

#' FRiP of each sample, from the size factor tables of step 7.
#'
#' @param workdir Pipeline WORKDIR.
#' @param marks Marks.
#' @return Data frame: sample, FRiP.
ReadFrip <- function(workdir, marks) {
  tables <- lapply(marks, function(mark) {
    file <- file.path(workdir, "ReadCount", mark,
                      paste0("sizeFactors_", mark, ".txt"))
    if (!file.exists(file)) return(NULL)
    return(as.data.frame(data.table::fread(file))[, c("sample", "FRiP")])
  })
  return(do.call(rbind, tables))
}

#' Peak numbers of MACS2 and FlexPeak and their overlap (step 4), if both ran.
#'
#' @param workdir Pipeline WORKDIR.
#' @param samples Sample IDs.
#' @return Data frame: sample, macs2_peaks, flexpeak_peaks, jaccard; or NULL.
ReadCallerComparison <- function(workdir, samples) {
  files <- file.path(workdir, "Peak", paste0("peak_", samples, ".callers.txt"))
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)
  tables <- lapply(files, function(f) as.data.frame(data.table::fread(f)))
  return(do.call(rbind, tables))
}

#' Bar plot of one QC value of all samples, one panel per mark.
#'
#' @param summary QC summary table.
#' @param column Column to plot.
#' @param y_label Axis label.
#' @param file Output PDF.
PlotQcValue <- function(summary, column, y_label, file) {
  if (!column %in% colnames(summary)) return(invisible(NULL))
  summary$value <- suppressWarnings(as.numeric(summary[[column]]))
  summary <- summary[!is.na(summary$value), ]
  if (nrow(summary) == 0) return(invisible(NULL))

  plot <- ggplot2::ggplot(summary,
                          ggplot2::aes(sample, value, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(~ mark, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = ConditionColors(summary$condition)) +
    ggplot2::labs(x = NULL, y = y_label, fill = NULL) +
    ThemeNature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  SaveFigure(plot, file, width = min(kDoubleColumn, 30 + 6 * nrow(summary)),
             height = 60)
  return(invisible(file))
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--workdir"),
  optparse::make_option("--sample-list", dest = "sample_list"),
  optparse::make_option("--igg-mark", dest = "igg_mark", default = "IgG")
)))
workdir <- opts$workdir
out_dir <- file.path(workdir, "QC", "Summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
meta <- ReadTable(opts$sample_list)

# 1. QC table of all samples
summary <- meta[, c("sample", "mark", "condition", "replicate", "batch")]
summary <- merge(summary, ReadAlignStats(workdir, meta$sample),
                 by = "sample", all.x = TRUE, sort = FALSE)
summary$n_peaks <- CountPeaks(workdir, summary$sample)
frip <- ReadFrip(workdir, setdiff(unique(meta$mark), opts$igg_mark))
if (!is.null(frip)) {
  summary <- merge(summary, frip, by = "sample", all.x = TRUE, sort = FALSE)
}
callers <- ReadCallerComparison(workdir, meta$sample)
if (!is.null(callers)) {
  summary <- merge(summary, callers, by = "sample", all.x = TRUE, sort = FALSE)
}
summary <- summary[match(meta$sample, summary$sample), ]
WriteTable(summary, file.path(out_dir, "sample_QC_summary.txt"))

# 2. Bar plots
summary$sample <- factor(summary$sample, levels = meta$sample)
qc_values <- c(unique_pairs = "Unique read pairs (fragments)",
               align_rate = "Alignment rate (%)",
               dup_rate = "Duplicate rate",
               n_peaks = "Number of peaks",
               FRiP = "FRiP (fraction of fragments in peaks)",
               jaccard = "MACS2 vs FlexPeak overlap (Jaccard)")
for (column in names(qc_values)) {
  PlotQcValue(summary, column, qc_values[[column]],
              file.path(out_dir, paste0("plot_", column, ".pdf")))
}

# 3. Fragment size distribution
fragments <- lapply(seq_len(nrow(meta)), function(i) {
  file <- file.path(workdir, "BAM",
                    paste0(meta$sample[i], ".fragment_size.txt"))
  if (!file.exists(file) || file.size(file) == 0) return(NULL)
  sizes <- as.data.frame(data.table::fread(file, header = FALSE))
  colnames(sizes) <- c("size", "n")
  sizes$fraction <- sizes$n / sum(sizes$n)
  sizes$sample <- meta$sample[i]
  sizes$mark <- meta$mark[i]
  sizes$condition <- meta$condition[i]
  return(sizes[sizes$size <= 1000, ])
})
fragments <- do.call(rbind, fragments)

if (!is.null(fragments)) {
  plot <- ggplot2::ggplot(
    fragments,
    ggplot2::aes(size, fraction, group = sample, color = condition)
  ) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::facet_wrap(~ mark, scales = "free_y") +
    ggplot2::scale_color_manual(values = ConditionColors(fragments$condition)) +
    ggplot2::labs(x = "Fragment size (bp)", y = "Fraction of fragments",
                  color = NULL, title = "Fragment size") +
    ThemeNature()
  SaveFigure(plot, file.path(out_dir, "plot_fragment_size.pdf"),
             width = kDoubleColumn, height = 60)
}
message("QC summary: ", file.path(out_dir, "sample_QC_summary.txt"))
