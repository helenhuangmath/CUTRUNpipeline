#!/usr/bin/env Rscript
#
# Plots of step 5 for one mark: Venn diagrams of replicates and conditions,
# replicate support of the merged peaks, width of the union peaks.
# Called by scripts/05_consensus.sh.

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

#' Venn diagram of 2 to 5 sets of peak IDs.
#'
#' @param sets Named list of vectors of peak IDs.
#' @param title Plot title.
#' @param file Output PDF.
PlotVenn <- function(sets, title, file) {
  if (length(sets) < 2 || length(sets) > 5) {
    return(invisible(NULL))
  }
  tryCatch({
    venn <- VennDiagram::venn.diagram(
      sets, filename = NULL, alpha = 0.5, cex = 1.3, cat.cex = 1.1,
      fill = kConditionPalette[seq_along(sets)], fontfamily = "sans",
      cat.fontfamily = "sans", main.fontfamily = "sans",
      main = title, main.cex = 1.4, margin = 0.1
    )
    grDevices::pdf(file, width = 6, height = 6)
    grid::grid.draw(venn)
    grDevices::dev.off()
  }, error = function(e) message("Venn diagram skipped: ", conditionMessage(e)))
  return(invisible(file))
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--dir", help = "Consensus/<mark> folder"),
  optparse::make_option("--mark"),
  optparse::make_option("--conditions", help = "comma-separated"),
  optparse::make_option("--min-reps", dest = "min_reps", type = "integer",
                        default = 2)
)))
mark <- opts$mark
conditions <- strsplit(opts$conditions, ",")[[1]]
InFile <- function(name) file.path(opts$dir, name)
invisible(futile.logger::flog.threshold(futile.logger::ERROR,
                                        name = "VennDiagramLogger"))

# 1. Replicates of each condition: Venn diagram and replicate support
support <- NULL
for (cond in conditions) {
  merged <- data.table::fread(InFile(paste0("peak_", cond, ".rep_support.bed")),
                              header = FALSE)
  if (nrow(merged) == 0) next
  colnames(merged) <- c("Chr", "Start", "End", "n_reps", "reps")
  support <- rbind(support, data.frame(condition = cond,
                                       n_reps = merged$n_reps))

  reps_of_peak <- strsplit(merged$reps, ",")
  reps <- sort(unique(unlist(reps_of_peak)))
  sets <- lapply(reps, function(rep) {
    return(which(vapply(reps_of_peak, function(x) rep %in% x, logical(1))))
  })
  names(sets) <- reps
  PlotVenn(sets, paste(mark, cond, "replicate peaks"),
           InFile(paste0("plot_Venn_", mark, "_", cond, "_rep_peaks.pdf")))
}

if (!is.null(support)) {
  plot <- ggplot2::ggplot(support,
                          ggplot2::aes(factor(n_reps), fill = condition)) +
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::scale_fill_manual(values = ConditionColors(conditions)) +
    ggplot2::labs(x = "Number of replicates with the peak", y = "Peaks",
                  title = paste0(mark, ": consensus = in >= ", opts$min_reps,
                                 " replicates")) +
    ThemeNature()
  SaveFigure(plot, InFile(paste0("plot_", mark, "_rep_support.pdf")),
             height = 55)
}

# 2. Venn diagram of conditions (union peaks found in each condition)
sets <- lapply(conditions, function(cond) {
  return(readLines(InFile(paste0("peak_", cond, "_in_union.txt"))))
})
names(sets) <- conditions
PlotVenn(sets, paste(mark, "peaks of each condition"),
         InFile(paste0("plot_Venn_", mark, "_AllConds_peak.pdf")))

# 3. Number of peaks and width of union peaks
union <- data.table::fread(InFile(paste0("peak_", mark, "_union.bed")),
                           header = FALSE)
colnames(union) <- c("Chr", "Start", "End", "peakN", "Width", "Strand")
WriteTable(data.frame(mark = mark, peaks = c(conditions, "union"),
                      n = c(lengths(sets), nrow(union))),
           InFile(paste0("Npeaks_", mark, ".txt")))

plot <- ggplot2::ggplot(union, ggplot2::aes(Width)) +
  ggplot2::geom_histogram(bins = 100, fill = "grey40") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(x = "Peak width (bp)",
                title = paste0(mark, ": ", nrow(union), " union peaks, ",
                               "median width ", stats::median(union$Width),
                               " bp")) +
  ThemeNature()
SaveFigure(plot, InFile(paste0("plot_", mark, "_UnionALL_peaks_width.pdf")),
           height = 55)
