#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

project_root <- getwd()

read_tsv <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_tsv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

run_edgeR <- function(counts, gene_info, metadata, group_levels, contrasts, out_dir, dataset) {
  metadata$group <- factor(metadata$group, levels = group_levels)
  sample_ids <- metadata$sample_id
  count_matrix <- as.matrix(counts[, sample_ids, drop = FALSE])
  storage.mode(count_matrix) <- "integer"
  rownames(count_matrix) <- gene_info$row_id

  y <- DGEList(counts = count_matrix, samples = metadata, genes = gene_info)
  keep <- filterByExpr(y, group = metadata$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")

  fig_dir <- file.path(project_root, "results/figures", dataset)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(fig_dir, paste0(dataset, "_edgeR_MDS.png")), width = 1800, height = 1500, res = 220)
  group_cols <- as.integer(y$samples$group)
  plotMDS(y, labels = y$samples$sample_title, col = group_cols, main = paste(dataset, "bulk RNA-seq MDS"))
  legend("topright", legend = levels(y$samples$group), col = seq_along(levels(y$samples$group)), pch = 16, cex = 0.75)
  dev.off()

  design <- model.matrix(~ 0 + group, data = y$samples)
  colnames(design) <- sub("^group", "", colnames(design))

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)

  all_results <- list()
  summaries <- list()

  for (contrast_name in names(contrasts)) {
    contrast <- makeContrasts(contrasts = contrasts[[contrast_name]], levels = design)
    qlf <- glmQLFTest(fit, contrast = contrast)
    res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
    res$row_id <- rownames(res)
    res$dataset <- dataset
    res$comparison <- contrast_name
    res$direction <- ifelse(res$logFC > 0, "up", ifelse(res$logFC < 0, "down", "flat"))

    front_cols <- c("dataset", "comparison", "row_id")
    gene_cols <- intersect(c("gene_id", "gene_symbol"), colnames(res))
    stat_cols <- c("logFC", "logCPM", "F", "PValue", "FDR", "direction")
    res <- res[, c(front_cols, gene_cols, stat_cols, setdiff(colnames(res), c(front_cols, gene_cols, stat_cols)))]
    all_results[[contrast_name]] <- res

    summaries[[contrast_name]] <- data.frame(
      dataset = dataset,
      comparison = contrast_name,
      n_samples = ncol(y),
      n_genes_tested = nrow(y),
      n_FDR_0_05 = sum(res$FDR < 0.05),
      n_FDR_0_05_up = sum(res$FDR < 0.05 & res$logFC > 0),
      n_FDR_0_05_down = sum(res$FDR < 0.05 & res$logFC < 0),
      n_FDR_0_10 = sum(res$FDR < 0.10),
      n_FDR_0_10_up = sum(res$FDR < 0.10 & res$logFC > 0),
      n_FDR_0_10_down = sum(res$FDR < 0.10 & res$logFC < 0)
    )
  }

  all_results_df <- do.call(rbind, all_results)
  summary_df <- do.call(rbind, summaries)
  write_tsv_gz(all_results_df, file.path(out_dir, paste0(dataset, "_edgeR_all_results.tsv.gz")))
  write_tsv(summary_df, file.path(out_dir, paste0(dataset, "_edgeR_summary.tsv")))
}

run_gse126005 <- function() {
  counts_path <- file.path(project_root, "data/processed/GSE126005/GSE126005_counts_gene.tsv.gz")
  metadata_path <- file.path(project_root, "metadata/GSE126005/samples.tsv")
  out_dir <- file.path(project_root, "results/tables/GSE126005")

  counts <- read_tsv(counts_path)
  metadata <- read_tsv(metadata_path)
  gene_info <- data.frame(
    row_id = make.unique(counts$gene_id),
    gene_id = counts$gene_id,
    gene_symbol = counts$gene_symbol,
    stringsAsFactors = FALSE
  )

  run_edgeR(
    counts = counts,
    gene_info = gene_info,
    metadata = metadata,
    group_levels = c("normal_control", "peyronies_plaque"),
    contrasts = c(PD_plaque_vs_control = "peyronies_plaque - normal_control"),
    out_dir = out_dir,
    dataset = "GSE126005"
  )
}

run_gse146500 <- function() {
  counts_path <- file.path(project_root, "data/processed/GSE146500/GSE146500_counts_gene_symbol.tsv.gz")
  metadata_path <- file.path(project_root, "metadata/GSE146500/samples.tsv")
  out_dir <- file.path(project_root, "results/tables/GSE146500")

  counts <- read_tsv(counts_path)
  metadata <- read_tsv(metadata_path)
  gene_info <- data.frame(
    row_id = make.unique(counts$gene_symbol),
    gene_symbol = counts$gene_symbol,
    stringsAsFactors = FALSE
  )

  run_edgeR(
    counts = counts,
    gene_info = gene_info,
    metadata = metadata,
    group_levels = c("normal_fibroblast", "peyronies_fibroblast", "peyronies_fibroblast_PC_NV"),
    contrasts = c(
      PF_vs_NF = "peyronies_fibroblast - normal_fibroblast",
      PFPC_vs_PF = "peyronies_fibroblast_PC_NV - peyronies_fibroblast",
      PFPC_vs_NF = "peyronies_fibroblast_PC_NV - normal_fibroblast"
    ),
    out_dir = out_dir,
    dataset = "GSE146500"
  )
}

run_gse126005()
run_gse146500()
