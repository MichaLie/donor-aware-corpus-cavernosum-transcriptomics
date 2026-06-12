#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

counts_file <- file.path(
  root,
  "data",
  "processed",
  "GSE206528",
  "pseudobulk_mural_subtypes",
  "GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
)
meta_file <- file.path(
  root,
  "data",
  "processed",
  "GSE206528",
  "pseudobulk_mural_subtypes",
  "GSE206528_mural_subtype_pseudobulk_metadata.tsv"
)
out_dir <- file.path(root, "results", "tables", "GSE206528_mural_subtype_edger")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_cells <- 50
min_group_n <- 2

counts_dt <- fread(cmd = paste("gzip -dc", shQuote(counts_file)))
gene_symbols <- counts_dt$gene_symbol
counts <- as.matrix(counts_dt[, setdiff(names(counts_dt), "gene_symbol"), with = FALSE])
rownames(counts) <- gene_symbols
storage.mode(counts) <- "integer"

meta <- as.data.frame(fread(meta_file))
rownames(meta) <- meta$pseudobulk_id
meta <- meta[colnames(counts), , drop = FALSE]

run_edgeR <- function(mural_subtype, comparison_name, selected_meta, group_values, group_levels) {
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
  res$mural_subtype <- mural_subtype
  res$comparison <- comparison_name
  res$n_genes_tested <- nrow(res)
  res$n_samples <- nrow(selected_meta)
  res$group_summary <- paste(names(group_table), as.integer(group_table), sep = "=", collapse = ";")
  res[, c("mural_subtype", "comparison", "gene_symbol", "logFC", "logCPM", "F", "PValue", "FDR", "n_genes_tested", "n_samples", "group_summary")]
}

all_results <- list()
diagnostics <- list()

for (subtype in sort(unique(meta$mural_subtype))) {
  subtype_meta <- meta[meta$mural_subtype == subtype & meta$n_cells >= min_cells, , drop = FALSE]
  diagnostics[[subtype]] <- data.frame(
    mural_subtype = subtype,
    pseudobulks_total = sum(meta$mural_subtype == subtype),
    pseudobulks_passing_min_cells = nrow(subtype_meta),
    donors_passing = paste(subtype_meta$donor_id, collapse = ","),
    min_cells = ifelse(nrow(subtype_meta) > 0, min(subtype_meta$n_cells), NA_integer_),
    max_cells = ifelse(nrow(subtype_meta) > 0, max(subtype_meta$n_cells), NA_integer_)
  )

  if (nrow(subtype_meta) < 2 * min_group_n) next

  all_vs_normal_group <- ifelse(subtype_meta$primary_group == "normal", "normal", "organic_ED")
  names(all_vs_normal_group) <- rownames(subtype_meta)
  result <- run_edgeR(
    subtype,
    "organic_ED_vs_normal",
    subtype_meta,
    all_vs_normal_group,
    c("normal", "organic_ED")
  )
  if (!is.null(result)) {
    all_results[[paste(subtype, "organic_ED_vs_normal", sep = "__")]] <- result
  }

  non_dm_meta <- subtype_meta[subtype_meta$ed_subgroup %in% c("normal", "non_diabetic_ED"), , drop = FALSE]
  non_dm_group <- non_dm_meta$ed_subgroup
  names(non_dm_group) <- rownames(non_dm_meta)
  result <- run_edgeR(
    subtype,
    "non_diabetic_ED_vs_normal",
    non_dm_meta,
    non_dm_group,
    c("normal", "non_diabetic_ED")
  )
  if (!is.null(result)) {
    all_results[[paste(subtype, "non_diabetic_ED_vs_normal", sep = "__")]] <- result
  }

  dm_meta <- subtype_meta[subtype_meta$ed_subgroup %in% c("normal", "diabetic_ED"), , drop = FALSE]
  dm_group <- dm_meta$ed_subgroup
  names(dm_group) <- rownames(dm_meta)
  result <- run_edgeR(
    subtype,
    "diabetic_ED_vs_normal",
    dm_meta,
    dm_group,
    c("normal", "diabetic_ED")
  )
  if (!is.null(result)) {
    all_results[[paste(subtype, "diabetic_ED_vs_normal", sep = "__")]] <- result
  }
}

diagnostics_df <- rbindlist(diagnostics, fill = TRUE)
fwrite(diagnostics_df, file.path(out_dir, "diagnostics_mural_subtype_edgeR.tsv"), sep = "\t")

if (length(all_results) == 0) {
  stop("No refined mural subtype comparisons could be run with current thresholds.")
}

results_df <- rbindlist(all_results, fill = TRUE)
fwrite(results_df, file.path(out_dir, "all_mural_subtype_edgeR_results.tsv.gz"), sep = "\t")

summary_df <- results_df[, .(
  genes_tested = unique(n_genes_tested)[1],
  n_samples = unique(n_samples)[1],
  group_summary = unique(group_summary)[1],
  n_FDR_0_05 = sum(FDR < 0.05, na.rm = TRUE),
  n_FDR_0_10 = sum(FDR < 0.10, na.rm = TRUE),
  n_abs_logFC_1_FDR_0_05 = sum(abs(logFC) >= 1 & FDR < 0.05, na.rm = TRUE)
), by = .(mural_subtype, comparison)]
fwrite(summary_df, file.path(out_dir, "summary_mural_subtype_edgeR_results.tsv"), sep = "\t")

message("Wrote refined mural subtype edgeR results to: ", out_dir)
