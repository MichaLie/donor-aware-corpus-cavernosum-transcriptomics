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
full_results_file <- file.path(
  root,
  "results",
  "tables",
  "GSE206528_mural_subtype_edger",
  "all_mural_subtype_edgeR_results.tsv.gz"
)
out_dir <- file.path(root, "results", "tables", "GSE206528_mural_subtype_edger_lodo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_cells <- 50
min_group_n <- 2
candidate_fdr_threshold <- 0.10

read_table <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)))
  } else {
    fread(path)
  }
}

counts_dt <- read_table(counts_file)
gene_symbols <- counts_dt$gene_symbol
counts <- as.matrix(counts_dt[, setdiff(names(counts_dt), "gene_symbol"), with = FALSE])
rownames(counts) <- gene_symbols
storage.mode(counts) <- "integer"

meta <- as.data.frame(fread(meta_file))
rownames(meta) <- meta$pseudobulk_id
meta <- meta[colnames(counts), , drop = FALSE]

full_results <- read_table(full_results_file)
candidates <- full_results[FDR < candidate_fdr_threshold]
candidate_keys <- unique(candidates[, .(
  mural_subtype,
  comparison,
  gene_symbol,
  full_logFC = logFC,
  full_PValue = PValue,
  full_FDR = FDR
)])

comparison_specs <- list(
  organic_ED_vs_normal = list(group_column = "primary_group", levels = c("normal", "organic_ED")),
  non_diabetic_ED_vs_normal = list(group_column = "ed_subgroup", levels = c("normal", "non_diabetic_ED")),
  diabetic_ED_vs_normal = list(group_column = "ed_subgroup", levels = c("normal", "diabetic_ED"))
)

run_edgeR <- function(selected_meta, selected_counts, group_column, group_levels) {
  selected_meta$comparison_group <- factor(selected_meta[[group_column]], levels = group_levels)
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
  res <- topTags(qlf, n = Inf, sort.by = "none")$table
  res$gene_symbol <- rownames(res)
  res$n_genes_tested <- nrow(res)
  res$n_samples <- nrow(selected_meta)
  res$group_summary <- paste(names(group_table), as.integer(group_table), sep = "=", collapse = ";")
  res
}

lodo_results <- list()
lodo_diagnostics <- list()

for (comparison_name in names(comparison_specs)) {
  spec <- comparison_specs[[comparison_name]]
  for (subtype in sort(unique(meta$mural_subtype))) {
    subtype_candidates <- candidate_keys[mural_subtype == subtype & comparison == comparison_name]
    if (nrow(subtype_candidates) == 0) next

    subtype_meta <- meta[
      meta$mural_subtype == subtype &
        meta$n_cells >= min_cells &
        meta[[spec$group_column]] %in% spec$levels,
      ,
      drop = FALSE
    ]
    subtype_counts <- counts[, rownames(subtype_meta), drop = FALSE]

    for (leave_out_donor in sort(unique(subtype_meta$donor_id))) {
      lodo_meta <- subtype_meta[subtype_meta$donor_id != leave_out_donor, , drop = FALSE]
      lodo_counts <- subtype_counts[, rownames(lodo_meta), drop = FALSE]
      group_values <- lodo_meta[[spec$group_column]]
      group_table <- table(factor(group_values, levels = spec$levels))

      diagnostic <- data.frame(
        mural_subtype = subtype,
        comparison = comparison_name,
        leave_out_donor = leave_out_donor,
        runnable = length(group_table[group_table > 0]) == 2 && all(group_table >= min_group_n),
        group_summary = paste(names(group_table), as.integer(group_table), sep = "=", collapse = ";"),
        stringsAsFactors = FALSE
      )

      res <- run_edgeR(lodo_meta, lodo_counts, spec$group_column, spec$levels)
      if (is.null(res)) {
        lodo_diagnostics[[paste(comparison_name, subtype, leave_out_donor, sep = "__")]] <- diagnostic
        next
      }

      diagnostic$runnable <- TRUE
      diagnostic$n_genes_tested <- unique(res$n_genes_tested)[1]
      diagnostic$n_samples <- unique(res$n_samples)[1]
      lodo_diagnostics[[paste(comparison_name, subtype, leave_out_donor, sep = "__")]] <- diagnostic

      res_dt <- as.data.table(res)
      merged <- merge(subtype_candidates, res_dt, by = "gene_symbol", all.x = TRUE)
      merged[, `:=`(
        mural_subtype = subtype,
        comparison = comparison_name,
        leave_out_donor = leave_out_donor,
        lodo_logFC = logFC,
        lodo_logCPM = logCPM,
        lodo_F = get("F"),
        lodo_PValue = PValue,
        lodo_FDR = FDR,
        lodo_n_genes_tested = n_genes_tested,
        lodo_n_samples = n_samples,
        lodo_group_summary = group_summary,
        lodo_gene_tested = !is.na(logFC)
      )]

      keep_cols <- c(
        "mural_subtype", "comparison", "gene_symbol", "leave_out_donor",
        "full_logFC", "full_PValue", "full_FDR",
        "lodo_logFC", "lodo_logCPM", "lodo_F", "lodo_PValue", "lodo_FDR",
        "lodo_n_genes_tested", "lodo_n_samples", "lodo_group_summary", "lodo_gene_tested"
      )
      lodo_results[[paste(comparison_name, subtype, leave_out_donor, sep = "__")]] <- merged[, ..keep_cols]
    }
  }
}

diagnostics_df <- rbindlist(lodo_diagnostics, fill = TRUE)
fwrite(diagnostics_df, file.path(out_dir, "lodo_diagnostics.tsv"), sep = "\t")

if (length(lodo_results) == 0) {
  stop("No leave-one-donor-out candidate comparisons could be run.")
}

lodo_df <- rbindlist(lodo_results, fill = TRUE)
lodo_df[, same_direction := !is.na(lodo_logFC) & sign(lodo_logFC) == sign(full_logFC)]
lodo_df[, lodo_FDR_0_05 := !is.na(lodo_FDR) & lodo_FDR < 0.05]
lodo_df[, lodo_FDR_0_10 := !is.na(lodo_FDR) & lodo_FDR < 0.10]
lodo_df[, lodo_nominal_0_05 := !is.na(lodo_PValue) & lodo_PValue < 0.05]

fwrite(lodo_df, file.path(out_dir, "candidate_lodo_edgeR_results.tsv.gz"), sep = "\t")

summary_df <- lodo_df[, .(
  full_logFC = unique(full_logFC)[1],
  full_PValue = unique(full_PValue)[1],
  full_FDR = unique(full_FDR)[1],
  n_lodo_runnable = .N,
  n_lodo_gene_tested = sum(lodo_gene_tested, na.rm = TRUE),
  n_lodo_same_direction = sum(same_direction, na.rm = TRUE),
  fraction_lodo_same_direction = sum(same_direction, na.rm = TRUE) / .N,
  n_lodo_nominal_0_05 = sum(lodo_nominal_0_05, na.rm = TRUE),
  n_lodo_FDR_0_10 = sum(lodo_FDR_0_10, na.rm = TRUE),
  n_lodo_FDR_0_05 = sum(lodo_FDR_0_05, na.rm = TRUE),
  min_lodo_logFC = suppressWarnings(min(lodo_logFC, na.rm = TRUE)),
  max_lodo_logFC = suppressWarnings(max(lodo_logFC, na.rm = TRUE)),
  leave_out_donors = paste(leave_out_donor, collapse = ",")
), by = .(mural_subtype, comparison, gene_symbol)]

summary_df[is.infinite(min_lodo_logFC), min_lodo_logFC := NA_real_]
summary_df[is.infinite(max_lodo_logFC), max_lodo_logFC := NA_real_]
summary_df[, stable_direction_all_lodo := n_lodo_runnable > 0 & n_lodo_same_direction == n_lodo_runnable]
summary_df[, gene_tested_all_lodo := n_lodo_gene_tested == n_lodo_runnable]
summary_df[, stable_tested_direction_all_lodo := gene_tested_all_lodo & n_lodo_same_direction == n_lodo_runnable]
summary_df[, stable_direction_min_80pct := fraction_lodo_same_direction >= 0.80]
summary_df[, retained_nominal_signal_any_lodo := n_lodo_nominal_0_05 > 0]

setorder(summary_df, mural_subtype, comparison, full_FDR, full_PValue, gene_symbol)
fwrite(summary_df, file.path(out_dir, "candidate_lodo_summary.tsv"), sep = "\t")

overview_df <- summary_df[, .(
  n_candidates = .N,
  n_gene_tested_all_lodo = sum(gene_tested_all_lodo, na.rm = TRUE),
  n_stable_tested_direction_all_lodo = sum(stable_tested_direction_all_lodo, na.rm = TRUE),
  n_stable_direction_all_lodo = sum(stable_direction_all_lodo, na.rm = TRUE),
  n_stable_direction_min_80pct = sum(stable_direction_min_80pct, na.rm = TRUE),
  n_retained_nominal_signal_any_lodo = sum(retained_nominal_signal_any_lodo, na.rm = TRUE),
  median_fraction_same_direction = median(fraction_lodo_same_direction, na.rm = TRUE)
), by = .(mural_subtype, comparison)]
setorder(overview_df, mural_subtype, comparison)
fwrite(overview_df, file.path(out_dir, "candidate_lodo_overview.tsv"), sep = "\t")

message("Wrote refined mural leave-one-donor-out candidate results to: ", out_dir)
