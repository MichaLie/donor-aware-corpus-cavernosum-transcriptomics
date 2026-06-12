#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(limma)
})

project_root <- getwd()

read_tsv <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

read_tsv_gz <- function(path) {
  read.delim(gzfile(path), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

split_genes <- function(x) {
  genes <- unique(trimws(unlist(strsplit(x, ","))))
  genes <- genes[!is.na(genes) & genes != ""]
  toupper(genes)
}

collapse_stat_df <- function(df) {
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "" & !is.na(df$F), ]
  df$signed_sqrt_F <- sign(df$logFC) * sqrt(df$F)
  df <- df[is.finite(df$signed_sqrt_F), ]
  df$gene_symbol <- toupper(df$gene_symbol)
  df <- df[order(abs(df$signed_sqrt_F), decreasing = TRUE), ]
  df <- df[!duplicated(df$gene_symbol), ]
  df$rank_desc <- rank(-df$signed_sqrt_F, ties.method = "average")
  df$rank_abs <- rank(-abs(df$signed_sqrt_F), ties.method = "average")
  df$rank_percentile_up <- 1 - ((df$rank_desc - 1) / max(1, nrow(df) - 1))
  df
}

camera_for_genes <- function(stat_df, genes, min_set_size = 3) {
  universe <- stat_df$gene_symbol
  genes <- intersect(unique(toupper(genes)), universe)
  if (length(genes) < min_set_size) {
    return(NULL)
  }
  statistic <- stat_df$signed_sqrt_F
  names(statistic) <- stat_df$gene_symbol
  statistic <- statistic[order(statistic, decreasing = TRUE)]
  universe <- names(statistic)
  cam <- cameraPR(statistic = statistic, index = match(genes, universe), use.ranks = FALSE)
  list(
    n = length(genes),
    direction = as.character(cam$Direction[1]),
    p_value = as.numeric(cam$PValue[1]),
    genes = paste(sort(genes), collapse = ",")
  )
}

expected_mode <- function(peyronie_comparison, camera_direction) {
  if (peyronie_comparison == "PFPC_vs_PF") {
    if (camera_direction == "Down") return("treatment_opposes_ED_direction")
    if (camera_direction == "Up") return("treatment_matches_ED_direction")
    return("treatment_ambiguous")
  }
  if (camera_direction == "Up") return("disease_same_direction")
  if (camera_direction == "Down") return("disease_opposite_direction")
  "disease_ambiguous"
}

out_dir <- file.path(project_root, "results/tables/robustness_checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gene_sets_path <- file.path(project_root, "metadata/gene_sets/robustness_module_gene_sets.tsv")
if (!file.exists(gene_sets_path)) {
  stop("Missing gene-set file: ", gene_sets_path, call. = FALSE)
}
gene_sets <- read_tsv(gene_sets_path)
gene_sets$gene_vector <- lapply(gene_sets$gene_symbols, split_genes)

ed <- read_tsv_gz(file.path(project_root, "results/tables/GSE206528_mural_subtype_edger/all_mural_subtype_edgeR_results.tsv.gz"))
ed$gene_symbol <- toupper(ed$gene_symbol)

module_rows <- list()
row_i <- 1
for (subtype in unique(ed$mural_subtype)) {
  for (ed_comp in unique(ed$comparison[ed$mural_subtype == subtype])) {
    stat_df <- collapse_stat_df(ed[ed$mural_subtype == subtype & ed$comparison == ed_comp, ])
    if (nrow(stat_df) == 0) next
    for (gs_i in seq_len(nrow(gene_sets))) {
      cam <- camera_for_genes(stat_df, gene_sets$gene_vector[[gs_i]], min_set_size = 3)
      if (is.null(cam)) next
      module_rows[[row_i]] <- data.frame(
        dataset = "GSE206528",
        mural_subtype = subtype,
        comparison = ed_comp,
        program_id = gene_sets$program_id[gs_i],
        program_label = gene_sets$program_label[gs_i],
        program_family = gene_sets$program_family[gs_i],
        source_status = gene_sets$source_status[gs_i],
        n_genes_in_universe = cam$n,
        camera_direction = cam$direction,
        camera_p_value = cam$p_value,
        genes_in_universe = cam$genes,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1
    }
  }
}

module_results <- if (length(module_rows)) do.call(rbind, module_rows) else data.frame()
if (nrow(module_results) > 0) {
  module_results$camera_FDR_global <- p.adjust(module_results$camera_p_value, method = "BH")
  group_key <- paste(module_results$mural_subtype, module_results$comparison, sep = "||")
  module_results$camera_FDR_within_context <- unsplit(
    lapply(split(module_results$camera_p_value, group_key), function(x) p.adjust(x, method = "BH")),
    group_key
  )
  module_results <- module_results[order(module_results$mural_subtype, module_results$comparison, module_results$camera_FDR_within_context), ]
}
write_tsv(module_results, file.path(out_dir, "GSE206528_refined_mural_module_cameraPR_results.tsv"))

ed_sig <- ed[
  ed$mural_subtype == "Pericyte_like" &
    ed$comparison == "diabetic_ED_vs_normal" &
    ed$FDR < 0.05,
]
ed_sig$ed_direction <- ifelse(ed_sig$logFC > 0, "up_in_diabetic_ED", "down_in_diabetic_ED")
candidate_genes_all <- unique(ed_sig$gene_symbol)
candidate_genes_up <- unique(ed_sig$gene_symbol[ed_sig$logFC > 0])

pey_a <- read_tsv_gz(file.path(project_root, "results/tables/GSE126005/GSE126005_edgeR_all_results.tsv.gz"))
pey_b <- read_tsv_gz(file.path(project_root, "results/tables/GSE146500/GSE146500_edgeR_all_results.tsv.gz"))
common_cols <- Reduce(intersect, list(colnames(pey_a), colnames(pey_b)))
pey <- rbind(pey_a[, common_cols], pey_b[, common_cols])
pey$gene_symbol <- toupper(pey$gene_symbol)

loo_rows <- list()
rank_rows <- list()
loo_i <- 1
rank_i <- 1

for (dataset in unique(pey$dataset)) {
  for (pey_comp in unique(pey$comparison[pey$dataset == dataset])) {
    stat_df <- collapse_stat_df(pey[pey$dataset == dataset & pey$comparison == pey_comp, ])
    universe <- stat_df$gene_symbol
    genes <- intersect(candidate_genes_up, universe)
    cam_full <- camera_for_genes(stat_df, genes, min_set_size = 3)
    if (!is.null(cam_full)) {
      loo_rows[[loo_i]] <- data.frame(
        dataset = dataset,
        peyronie_comparison = pey_comp,
        omitted_gene = "NONE_FULL_SIGNATURE",
        n_signature_genes_in_universe = cam_full$n,
        camera_direction = cam_full$direction,
        camera_p_value = cam_full$p_value,
        validation_mode = expected_mode(pey_comp, cam_full$direction),
        remaining_genes = cam_full$genes,
        stringsAsFactors = FALSE
      )
      loo_i <- loo_i + 1
      for (gene in sort(genes)) {
        remaining <- setdiff(genes, gene)
        cam <- camera_for_genes(stat_df, remaining, min_set_size = 3)
        if (is.null(cam)) next
        loo_rows[[loo_i]] <- data.frame(
          dataset = dataset,
          peyronie_comparison = pey_comp,
          omitted_gene = gene,
          n_signature_genes_in_universe = cam$n,
          camera_direction = cam$direction,
          camera_p_value = cam$p_value,
          validation_mode = expected_mode(pey_comp, cam$direction),
          remaining_genes = cam$genes,
          stringsAsFactors = FALSE
        )
        loo_i <- loo_i + 1
      }
    }

    for (gene in sort(candidate_genes_all)) {
      ed_gene <- ed_sig[ed_sig$gene_symbol == gene, ]
      pey_gene <- stat_df[stat_df$gene_symbol == gene, ]
      if (nrow(pey_gene) == 0) {
        rank_rows[[rank_i]] <- data.frame(
          dataset = dataset,
          peyronie_comparison = pey_comp,
          gene_symbol = gene,
          in_peyronie_universe = FALSE,
          ed_logFC = ed_gene$logFC[1],
          ed_FDR = ed_gene$FDR[1],
          ed_direction = ed_gene$ed_direction[1],
          peyronie_logFC = NA,
          peyronie_FDR = NA,
          signed_sqrt_F = NA,
          rank_desc = NA,
          rank_abs = NA,
          rank_percentile_up = NA,
          same_direction_as_ed = NA,
          stringsAsFactors = FALSE
        )
      } else {
        same_dir <- sign(ed_gene$logFC[1]) == sign(pey_gene$logFC[1])
        rank_rows[[rank_i]] <- data.frame(
          dataset = dataset,
          peyronie_comparison = pey_comp,
          gene_symbol = gene,
          in_peyronie_universe = TRUE,
          ed_logFC = ed_gene$logFC[1],
          ed_FDR = ed_gene$FDR[1],
          ed_direction = ed_gene$ed_direction[1],
          peyronie_logFC = pey_gene$logFC[1],
          peyronie_FDR = pey_gene$FDR[1],
          signed_sqrt_F = pey_gene$signed_sqrt_F[1],
          rank_desc = pey_gene$rank_desc[1],
          rank_abs = pey_gene$rank_abs[1],
          rank_percentile_up = pey_gene$rank_percentile_up[1],
          same_direction_as_ed = same_dir,
          stringsAsFactors = FALSE
        )
      }
      rank_i <- rank_i + 1
    }
  }
}

loo_results <- if (length(loo_rows)) do.call(rbind, loo_rows) else data.frame()
if (nrow(loo_results) > 0) {
  group_key <- paste(loo_results$dataset, loo_results$peyronie_comparison, sep = "||")
  loo_results$camera_FDR_within_comparison <- unsplit(
    lapply(split(loo_results$camera_p_value, group_key), function(x) p.adjust(x, method = "BH")),
    group_key
  )
  loo_results <- loo_results[order(loo_results$dataset, loo_results$peyronie_comparison, loo_results$omitted_gene), ]
}
write_tsv(loo_results, file.path(out_dir, "peyronie_leave_one_gene_out_cameraPR.tsv"))

rank_results <- if (length(rank_rows)) do.call(rbind, rank_rows) else data.frame()
rank_results <- rank_results[order(rank_results$dataset, rank_results$peyronie_comparison, rank_results$gene_symbol), ]
write_tsv(rank_results, file.path(out_dir, "peyronie_candidate_gene_rank_table.tsv"))

cat("Wrote", nrow(module_results), "module cameraPR rows\n")
cat("Wrote", nrow(loo_results), "Peyronie leave-one-gene-out rows\n")
cat("Wrote", nrow(rank_results), "Peyronie per-gene rank rows\n")
