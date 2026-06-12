#!/usr/bin/env Rscript
# Inspect the dropEst/zUMIs dgecounts.rds structure for GSE131882 (one sample).
# We do NOT trust any external description of the format -- print it ourselves.
suppressMessages({library(Matrix)})
args <- commandArgs(trailingOnly=TRUE)
f <- if (length(args)>=1) args[1] else "data/raw/GSE131882/GSM3823939_control.s1.dgecounts.rds.gz"
cat("Reading:", f, "\n")
obj <- readRDS(gzfile(f))
cat("top-level class:", class(obj), "\n")
cat("top-level names:", paste(names(obj), collapse=", "), "\n\n")
str(obj, max.level=2, list.len=20)
# Drill into umicount if present
if ("umicount" %in% names(obj)) {
  cat("\n--- umicount sublist ---\n")
  cat("names(umicount):", paste(names(obj$umicount), collapse=", "), "\n")
  if ("exon" %in% names(obj$umicount)) {
    cat("names(umicount$exon):", paste(names(obj$umicount$exon), collapse=", "), "\n")
    m <- obj$umicount$exon$all
    cat("class(umicount$exon$all):", class(m), "\n")
    cat("dim (rows x cols):", paste(dim(m), collapse=" x "), "\n")
    cat("rownames head:", paste(head(rownames(m), 8), collapse=", "), "\n")
    cat("colnames head:", paste(head(colnames(m), 4), collapse=", "), "\n")
    cat("total UMIs:", sum(m), "\n")
  }
}
