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

collapse_stat_by_gene <- function(df) {
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "" & !is.na(df$F), ]
  df$signed_sqrt_F <- sign(df$logFC) * sqrt(df$F)
  split_idx <- split(seq_len(nrow(df)), df$gene_symbol)
  keep <- vapply(split_idx, function(idx) idx[which.max(abs(df$signed_sqrt_F[idx]))], integer(1))
  collapsed <- df[keep, ]
  statistic <- collapsed$signed_sqrt_F
  names(statistic) <- collapsed$gene_symbol
  statistic <- statistic[is.finite(statistic)]
  sort(statistic, decreasing = TRUE)
}

expected_mode <- function(peyronie_comparison, ed_direction, camera_direction) {
  if (peyronie_comparison == "PFPC_vs_PF") {
    if (ed_direction == "up_in_ED" && camera_direction == "Down") return("treatment_opposes_ED_direction")
    if (ed_direction == "down_in_ED" && camera_direction == "Up") return("treatment_opposes_ED_direction")
    if (ed_direction == "up_in_ED" && camera_direction == "Up") return("treatment_matches_ED_direction")
    if (ed_direction == "down_in_ED" && camera_direction == "Down") return("treatment_matches_ED_direction")
    return("treatment_ambiguous")
  }

  if (ed_direction == "up_in_ED" && camera_direction == "Up") return("disease_same_direction")
  if (ed_direction == "down_in_ED" && camera_direction == "Down") return("disease_same_direction")
  if (ed_direction == "up_in_ED" && camera_direction == "Down") return("disease_opposite_direction")
  if (ed_direction == "down_in_ED" && camera_direction == "Up") return("disease_opposite_direction")
  "disease_ambiguous"
}

ed <- read_tsv(file.path(
  project_root,
  "results/tables/GSE206528_pseudobulk_edger_lodo/robust_FDR_0_05_candidates_with_programs.tsv"
))
ed <- ed[!is.na(ed$gene_symbol) & ed$gene_symbol != "", ]

pey_a <- read_tsv_gz(file.path(project_root, "results/tables/GSE126005/GSE126005_edgeR_all_results.tsv.gz"))
pey_b <- read_tsv_gz(file.path(project_root, "results/tables/GSE146500/GSE146500_edgeR_all_results.tsv.gz"))
common_cols <- Reduce(intersect, list(colnames(pey_a), colnames(pey_b)))
pey <- rbind(pey_a[, common_cols], pey_b[, common_cols])

rows <- list()
row_i <- 1
min_set_size <- 3

for (dataset in unique(pey$dataset)) {
  for (pey_comp in unique(pey$comparison[pey$dataset == dataset])) {
    pey_sub <- pey[pey$dataset == dataset & pey$comparison == pey_comp, ]
    stat <- collapse_stat_by_gene(pey_sub)
    universe <- names(stat)

    for (cell_type in unique(ed$cell_type)) {
      for (ed_comp in unique(ed$comparison[ed$cell_type == cell_type])) {
        for (ed_direction in c("up_in_ED", "down_in_ED")) {
          genes <- unique(ed$gene_symbol[
            ed$cell_type == cell_type &
              ed$comparison == ed_comp &
              ed$direction == ed_direction
          ])
          genes <- intersect(genes, universe)
          if (length(genes) < min_set_size) next

          idx <- match(genes, universe)
          cam <- cameraPR(statistic = stat, index = idx, use.ranks = FALSE)
          camera_direction <- as.character(cam$Direction[1])
          mode <- expected_mode(pey_comp, ed_direction, camera_direction)

          rows[[row_i]] <- data.frame(
            dataset = dataset,
            peyronie_comparison = pey_comp,
            cell_type = cell_type,
            ed_comparison = ed_comp,
            ed_direction = ed_direction,
            n_signature_genes_in_universe = length(genes),
            camera_direction = camera_direction,
            camera_p_value = cam$PValue[1],
            validation_mode = mode,
            genes = paste(sort(genes), collapse = ","),
            stringsAsFactors = FALSE
          )
          row_i <- row_i + 1
        }
      }
    }
  }
}

results <- do.call(rbind, rows)
results$camera_FDR <- p.adjust(results$camera_p_value, method = "BH")
results <- results[order(results$camera_FDR, results$camera_p_value), ]

summary <- aggregate(
  camera_FDR ~ dataset + peyronie_comparison + validation_mode,
  data = results,
  FUN = function(x) sum(x < 0.05)
)
names(summary)[names(summary) == "camera_FDR"] <- "n_camera_FDR_0_05"
summary <- summary[order(summary$dataset, summary$peyronie_comparison, summary$validation_mode), ]

out_dir <- file.path(project_root, "results/tables/peyronie_validation")
write_tsv(results, file.path(out_dir, "ed_signature_cameraPR_results.tsv"))
write_tsv(summary, file.path(out_dir, "ed_signature_cameraPR_mode_summary.tsv"))

cat("Wrote", nrow(results), "cameraPR rows\n")
