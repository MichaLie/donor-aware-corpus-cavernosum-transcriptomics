#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
counts_file <- file.path(root, "data", "processed", "GSE206528", "pseudobulk", "GSE206528_pseudobulk_counts_by_donor_celltype.tsv.gz")
meta_file <- file.path(root, "data", "processed", "GSE206528", "pseudobulk", "GSE206528_pseudobulk_sample_metadata.tsv")
out_dir <- file.path(root, "results", "tables", "GSE206528_pseudobulk_edger")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_cells <- 50
min_group_n <- 2

counts_dt <- if (grepl("\\.gz$", counts_file)) {
  fread(cmd = paste("gzip -dc", shQuote(counts_file)))
} else {
  fread(counts_file)
}
gene_symbols <- counts_dt$gene_symbol
counts <- as.matrix(counts_dt[, setdiff(names(counts_dt), "gene_symbol"), with = FALSE])
rownames(counts) <- gene_symbols
storage.mode(counts) <- "integer"

meta <- fread(meta_file)
meta <- as.data.frame(meta)
rownames(meta) <- meta$pseudobulk_id

stopifnot(all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), , drop = FALSE]

run_edgeR <- function(cell_type, comparison_name, keep_samples, group_values, group_levels) {
  selected_meta <- meta[keep_samples, , drop = FALSE]
  selected_counts <- counts[, rownames(selected_meta), drop = FALSE]

  selected_meta$comparison_group <- factor(group_values[rownames(selected_meta)], levels = group_levels)
  selected_meta <- selected_meta[!is.na(selected_meta$comparison_group), , drop = FALSE]
  selected_counts <- selected_counts[, rownames(selected_meta), drop = FALSE]
  selected_meta$comparison_group <- droplevels(selected_meta$comparison_group)

  group_table <- table(selected_meta$comparison_group)
  if (length(group_table) != 2 || any(group_table < min_group_n)) {
    return(NULL)
  }

  y <- DGEList(counts = selected_counts, samples = selected_meta, genes = data.frame(gene_symbol = rownames(selected_counts)))
  keep_genes <- filterByExpr(y, group = selected_meta$comparison_group)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)

  design <- model.matrix(~ comparison_group, data = selected_meta)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  res$gene_symbol <- rownames(res)
  res$cell_type <- cell_type
  res$comparison <- comparison_name
  res$n_genes_tested <- nrow(res)
  res$n_samples <- nrow(selected_meta)
  res$group_summary <- paste(names(group_table), as.integer(group_table), sep = "=", collapse = ";")
  res <- res[, c("cell_type", "comparison", "gene_symbol", "logFC", "logCPM", "F", "PValue", "FDR", "n_genes_tested", "n_samples", "group_summary")]
  res
}

all_results <- list()
diagnostics <- list()

cell_types <- sort(unique(meta$major_cell_type))
for (cell_type in cell_types) {
  cell_meta <- meta[meta$major_cell_type == cell_type & meta$n_cells >= min_cells, , drop = FALSE]

  diagnostics[[cell_type]] <- data.frame(
    cell_type = cell_type,
    pseudobulks_total = sum(meta$major_cell_type == cell_type),
    pseudobulks_passing_min_cells = nrow(cell_meta),
    donors_passing = paste(cell_meta$donor_id, collapse = ","),
    min_cells = ifelse(nrow(cell_meta) > 0, min(cell_meta$n_cells), NA_integer_),
    max_cells = ifelse(nrow(cell_meta) > 0, max(cell_meta$n_cells), NA_integer_)
  )

  if (nrow(cell_meta) < 2 * min_group_n) {
    next
  }

  sample_ids <- rownames(cell_meta)

  all_vs_normal_group <- ifelse(cell_meta$primary_group == "normal", "normal", "organic_ED")
  names(all_vs_normal_group) <- sample_ids
  result <- run_edgeR(
    cell_type,
    "organic_ED_vs_normal",
    sample_ids,
    all_vs_normal_group,
    c("normal", "organic_ED")
  )
  if (!is.null(result)) {
    all_results[[paste(cell_type, "organic_ED_vs_normal", sep = "__")]] <- result
  }

  non_dm_meta <- cell_meta[cell_meta$ed_subgroup %in% c("normal", "non_diabetic_ED"), , drop = FALSE]
  non_dm_group <- non_dm_meta$ed_subgroup
  names(non_dm_group) <- rownames(non_dm_meta)
  result <- run_edgeR(
    cell_type,
    "non_diabetic_ED_vs_normal",
    rownames(non_dm_meta),
    non_dm_group,
    c("normal", "non_diabetic_ED")
  )
  if (!is.null(result)) {
    all_results[[paste(cell_type, "non_diabetic_ED_vs_normal", sep = "__")]] <- result
  }

  dm_meta <- cell_meta[cell_meta$ed_subgroup %in% c("normal", "diabetic_ED"), , drop = FALSE]
  dm_group <- dm_meta$ed_subgroup
  names(dm_group) <- rownames(dm_meta)
  result <- run_edgeR(
    cell_type,
    "diabetic_ED_vs_normal",
    rownames(dm_meta),
    dm_group,
    c("normal", "diabetic_ED")
  )
  if (!is.null(result)) {
    all_results[[paste(cell_type, "diabetic_ED_vs_normal", sep = "__")]] <- result
  }
}

diagnostics_df <- rbindlist(diagnostics, fill = TRUE)
fwrite(diagnostics_df, file.path(out_dir, "diagnostics_pseudobulk_edgeR.tsv"), sep = "\t")

if (length(all_results) == 0) {
  stop("No edgeR comparisons could be run with the current thresholds.")
}

results_df <- rbindlist(all_results, fill = TRUE)
fwrite(results_df, file.path(out_dir, "all_pseudobulk_edgeR_results.tsv.gz"), sep = "\t")

summary_df <- results_df[, .(
  genes_tested = unique(n_genes_tested)[1],
  n_samples = unique(n_samples)[1],
  group_summary = unique(group_summary)[1],
  n_FDR_0_05 = sum(FDR < 0.05, na.rm = TRUE),
  n_FDR_0_10 = sum(FDR < 0.10, na.rm = TRUE),
  n_abs_logFC_1_FDR_0_05 = sum(abs(logFC) >= 1 & FDR < 0.05, na.rm = TRUE)
), by = .(cell_type, comparison)]
fwrite(summary_df, file.path(out_dir, "summary_pseudobulk_edgeR_results.tsv"), sep = "\t")

message("Wrote edgeR pseudobulk results to: ", out_dir)
