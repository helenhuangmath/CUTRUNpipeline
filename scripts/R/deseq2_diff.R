#!/usr/bin/env Rscript
#
# Step 8: differential peaks of one mark with DESeq2.
# Called by scripts/08_diff.sh.
#
# One DESeq2 model is fitted to all samples of the mark, using the raw count
# matrix and the size factors of step 7 (so the results agree with
# NormRC_<mark>.txt). Each comparison of the comparison list is a contrast.
#
# Output (in --outdir), for each comparison <mark>_<treatment>_vs_<reference>:
#   deseq2_All_<...>.txt, deseq2_Sig_<...>.txt   all / significant peaks
#   deseq2_Sig_<...>_Up.bed, _Down.bed           significant peaks as BED
#   plot_MA_, plot_Volcano_, plot_Heatmap_<...>.pdf
# and deseq2_summary_<mark>.txt, plot_SigDiff_Summary_Bar_<mark>.pdf,
# sessionInfo_<mark>.txt

script_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "",
                          grep("^--file=", script_args, value = TRUE)))
source(file.path(script_dir, "common.R"))

#' Builds the DESeq2 design. Terms with one level (e.g. one batch) are
#' dropped; if the design is confounded, ~ condition is used.
#'
#' @param design_text Design, e.g. "~ batch + condition".
#' @param col_data Sample information, one row per sample.
#' @return Design formula.
BuildDesign <- function(design_text, col_data) {
  terms <- all.vars(stats::as.formula(design_text))
  missing_terms <- setdiff(terms, colnames(col_data))
  if (length(missing_terms) > 0) {
    stop("DESIGN uses columns not in the sample list: ", missing_terms)
  }
  terms <- terms[vapply(terms, function(term) {
    return(length(unique(col_data[[term]])) > 1)
  }, logical(1))]
  terms <- c(setdiff(terms, "condition"), intersect(terms, "condition"))
  design <- stats::as.formula(paste("~", paste(c(1, terms), collapse = " + ")))

  model <- stats::model.matrix(design, col_data)
  if (qr(model)$rank < ncol(model)) {
    message("Design is confounded (e.g. batch = condition): use ~ condition")
    design <- ~ condition
  }
  return(design)
}

#' Runs DESeq2. If replicates are very similar the dispersion curve cannot be
#' fitted; then the gene-wise dispersions are used, as DESeq2 suggests.
#'
#' @param dds DESeqDataSet with size factors.
#' @return DESeqDataSet after the Wald test.
RunDeseq <- function(dds) {
  dds <- tryCatch(DESeq2::DESeq(dds), error = function(e) {
    if (!grepl("dispersion", conditionMessage(e))) stop(e)
    message("Dispersion fit failed (replicates very similar): ",
            "use gene-wise dispersion estimates")
    dds <- DESeq2::estimateDispersionsGeneEst(dds)
    DESeq2::dispersions(dds) <- S4Vectors::mcols(dds)$dispGeneEst
    return(DESeq2::nbinomWaldTest(dds))
  })
  return(dds)
}

#' Results of one comparison, with peak information and condition means.
#'
#' @param dds DESeqDataSet after DESeq().
#' @param treatment,reference Conditions compared.
#' @param peaks Peak information (peakN, Chr, Start, ...).
#' @param padj,lfc Significance cutoffs.
#' @return Data frame sorted by padj; direction is Up, Down or NS.
CompareConditions <- function(dds, treatment, reference, peaks, padj, lfc) {
  result <- as.data.frame(DESeq2::results(
    dds, contrast = c("condition", treatment, reference), alpha = padj
  ))
  significant <- !is.na(result$padj) & result$padj < padj &
    abs(result$log2FoldChange) > lfc
  result$direction <- "NS"
  result$direction[significant & result$log2FoldChange > 0] <- "Up"
  result$direction[significant & result$log2FoldChange < 0] <- "Down"
  result$comparison <- paste0(treatment, "_vs_", reference)

  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  for (cond in c(treatment, reference)) {
    in_cond <- dds$condition == cond
    result[[cond]] <- round(rowMeans(norm_counts[, in_cond, drop = FALSE]), 2)
  }

  result <- cbind(peaks[match(rownames(result), peaks$peakN), ], result)
  result <- result[order(result$padj, -abs(result$log2FoldChange)), ]
  return(result)
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--mark"),
  optparse::make_option("--rawrc", help = "raw count matrix of step 7"),
  optparse::make_option("--size-factors", dest = "size_factors"),
  optparse::make_option("--sample-list", dest = "sample_list"),
  optparse::make_option("--comparison-list", dest = "comparison_list",
                        default = ""),
  optparse::make_option("--design", default = "~ condition"),
  optparse::make_option("--min-count", dest = "min_count", type = "double",
                        default = 15),
  optparse::make_option("--padj", type = "double", default = 0.05),
  optparse::make_option("--lfc", type = "double", default = 0.585),
  optparse::make_option("--outdir")
)))
mark <- opts$mark
OutFile <- function(name) file.path(opts$outdir, name)

# 1. Raw count matrix, sample information and size factors of step 7
meta <- ReadTable(opts$sample_list)
meta <- meta[meta$mark == mark, ]
rownames(meta) <- meta$sample

raw_table <- as.data.frame(data.table::fread(opts$rawrc))
peaks <- raw_table[, setdiff(colnames(raw_table), meta$sample)]
raw_counts <- as.matrix(raw_table[, meta$sample, drop = FALSE])
rownames(raw_counts) <- raw_table$peakN

size_factors <- as.data.frame(data.table::fread(opts$size_factors))
size_factors <- stats::setNames(size_factors$size_factor, size_factors$sample)

# 2. Remove peaks with few reads in all samples
keep <- apply(raw_counts, 1, max) >= opts$min_count
message(mark, ": ", sum(keep), " of ", length(keep),
        " peaks have a max count >= ", opts$min_count)
WriteTable(peaks[!keep, ], OutFile(paste0("peak_", mark, "_removed_lowRC.txt")))

# 3. DESeq2 model of all samples of the mark
col_data <- meta[, setdiff(colnames(meta), c("fastq1", "fastq2"))]
col_data[] <- lapply(col_data, function(x) factor(x, levels = unique(x)))
design <- BuildDesign(opts$design, col_data)
message(mark, ": design ", deparse(design))

dds <- DESeq2::DESeqDataSetFromMatrix(raw_counts[keep, , drop = FALSE],
                                      col_data, design)
DESeq2::sizeFactors(dds) <- size_factors[colnames(dds)]

n_coefficients <- ncol(stats::model.matrix(design, col_data))
can_test <- "condition" %in% all.vars(design) &&
  ncol(dds) > n_coefficients
if (can_test) {
  dds <- RunDeseq(dds)
} else {
  message(mark, ": not enough conditions or replicates for DESeq2")
}
vst_counts <- VstMatrix(dds)

# 4. Each comparison
comparisons <- data.frame(mark = character())
if (IsFile(opts$comparison_list)) {
  comparisons <- ReadTable(opts$comparison_list)
}
comparisons <- comparisons[comparisons$mark %in% c(mark, "*"), ]
conditions <- levels(col_data$condition)
summary_rows <- list()

for (i in seq_len(nrow(comparisons))) {
  treatment <- comparisons$treatment[i]
  reference <- comparisons$reference[i]
  name <- paste0(mark, "_", treatment, "_vs_", reference)
  if (!can_test || !all(c(treatment, reference) %in% conditions)) {
    message("Skip ", name, ": condition not found or too few replicates")
    next
  }

  result <- CompareConditions(dds, treatment, reference, peaks, opts$padj,
                              opts$lfc)
  significant <- result[result$direction != "NS", ]
  WriteTable(result, OutFile(paste0("deseq2_All_", name, ".txt")))
  WriteTable(significant, OutFile(paste0("deseq2_Sig_", name, ".txt")))
  for (direction in c("Up", "Down")) {
    bed <- significant[significant$direction == direction,
                       c("Chr", "Start", "End", "peakN", "log2FoldChange")]
    data.table::fwrite(bed, OutFile(paste0("deseq2_Sig_", name, "_",
                                           direction, ".bed")),
                       sep = "\t", col.names = FALSE)
  }

  title <- paste(mark, treatment, "vs", reference)
  SaveFigure(MaPlot(result, title), OutFile(paste0("plot_MA_", name, ".pdf")))
  SaveFigure(VolcanoPlot(result, title),
             OutFile(paste0("plot_Volcano_", name, ".pdf")), height = 80)
  in_comparison <- meta$sample[meta$condition %in% c(treatment, reference)]
  PlotHeatmap(vst_counts[significant$peakN, in_comparison, drop = FALSE],
              meta[in_comparison, "condition"], title,
              OutFile(paste0("plot_Heatmap_", name, ".pdf")))

  summary_rows[[name]] <- data.frame(
    mark = mark, comparison = paste0(treatment, "_vs_", reference),
    tested = sum(!is.na(result$padj)), Up = sum(result$direction == "Up"),
    Down = sum(result$direction == "Down")
  )
  message(sprintf("%s: %d Up, %d Down (padj < %g, |log2FC| > %g)", name,
                  summary_rows[[name]]$Up, summary_rows[[name]]$Down,
                  opts$padj, opts$lfc))
}

# 5. Summary of all comparisons, R and package versions
if (length(summary_rows) > 0) {
  summary <- do.call(rbind, summary_rows)
  WriteTable(summary, OutFile(paste0("deseq2_summary_", mark, ".txt")))
  SaveFigure(DiffSummaryPlot(summary, paste("Differential peaks", mark)),
             OutFile(paste0("plot_SigDiff_Summary_Bar_", mark, ".pdf")),
             width = 30 + 15 * nrow(summary))
}
writeLines(utils::capture.output(utils::sessionInfo()),
           OutFile(paste0("sessionInfo_", mark, ".txt")))
