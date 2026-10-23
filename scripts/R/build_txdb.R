#!/usr/bin/env Rscript
#
# Gene annotation files from a GTF file. Called by reference/build_reference.sh.
#   --out-txdb  : TxDb sqlite for ChIPseeker
#   --out-genes : gene table: geneId, Gene.Name, gene.type, Chr, Start, End,
#                 Strand

#' Reads one attribute of the GTF column 9, e.g. gene_name.
#'
#' @param attributes GTF column 9.
#' @param key Attribute name.
#' @return Character vector, NA where the attribute is missing.
GtfAttribute <- function(attributes, key) {
  pattern <- paste0(".*", key, ' "([^"]+)".*')
  values <- sub(pattern, "\\1", attributes)
  values[!grepl(paste0(key, ' "'), attributes)] <- NA
  return(values)
}

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--gtf"),
  optparse::make_option("--out-txdb", dest = "out_txdb"),
  optparse::make_option("--out-genes", dest = "out_genes")
)))

# 1. TxDb (makeTxDbFromGFF moved from GenomicFeatures to txdbmaker in
# Bioconductor 3.19)
MakeTxDb <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
  txdbmaker::makeTxDbFromGFF
} else {
  GenomicFeatures::makeTxDbFromGFF
}
AnnotationDbi::saveDb(MakeTxDb(opts$gtf, format = "gtf"), opts$out_txdb)

# 2. Gene table from the "gene" lines
gtf <- data.table::fread(cmd = paste("gzip -dcf", shQuote(opts$gtf)),
                         sep = "\t", header = FALSE, quote = "")
gtf <- as.data.frame(gtf[gtf$V3 == "gene", ])
gene_type <- GtfAttribute(gtf$V9, "gene_type")
no_type <- is.na(gene_type)
gene_type[no_type] <- GtfAttribute(gtf$V9[no_type], "gene_biotype")  # Ensembl

genes <- data.frame(geneId = GtfAttribute(gtf$V9, "gene_id"),
                    Gene.Name = GtfAttribute(gtf$V9, "gene_name"),
                    gene.type = gene_type, Chr = gtf$V1, Start = gtf$V4,
                    End = gtf$V5, Strand = gtf$V7)
no_name <- is.na(genes$Gene.Name)
genes$Gene.Name[no_name] <- genes$geneId[no_name]
data.table::fwrite(genes, opts$out_genes, sep = "\t")
message("Saved ", opts$out_txdb, " and ", opts$out_genes, " (", nrow(genes),
        " genes)")
