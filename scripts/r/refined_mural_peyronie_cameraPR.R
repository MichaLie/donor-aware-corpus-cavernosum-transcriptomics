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
  sort(statistic[is.finite(statistic)], decreasing = TRUE)
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

ed <- read_tsv_gz(file.path(project_root, "results/tables/GSE206528_mural_subtype_edger/all_mural_subtype_edgeR_results.tsv.gz"))
ed <- ed[ed$FDR < 0.05 & !is.na(ed$gene_symbol) & ed$gene_symbol != "", ]
ed$direction <- ifelse(ed$logFC > 0, "up_in_ED", "down_in_ED")

pey_a <- read_tsv_gz(file.path(project_root, "results/tables/GSE126005/GSE126005_edgeR_all_results.tsv.gz"))
pey_b <- read_tsv_gz(file.path(project_root, "results/tables/GSE146500/GSE146500_edgeR_all_results.tsv.gz"))
common_cols <- Reduce(intersect, list(colnames(pey_a), colnames(pey_b)))
pey <- rbind(pey_a[, common_cols], pey_b[, common_cols])

rows <- list()
row_i <- 1
min_set_size <- 3

for (dataset in unique(pey$dataset)) {
  for (pey_comp in unique(pey$comparison[pey$dataset == dataset])) {
    stat <- collapse_stat_by_gene(pey[pey$dataset == dataset & pey$comparison == pey_comp, ])
    universe <- names(stat)

    for (subtype in unique(ed$mural_subtype)) {
      for (ed_comp in unique(ed$comparison[ed$mural_subtype == subtype])) {
        for (ed_direction in unique(ed$direction[ed$mural_subtype == subtype & ed$comparison == ed_comp])) {
          genes <- unique(ed$gene_symbol[
            ed$mural_subtype == subtype &
              ed$comparison == ed_comp &
              ed$direction == ed_direction
          ])
          genes <- intersect(genes, universe)
          if (length(genes) < min_set_size) next

          cam <- cameraPR(statistic = stat, index = match(genes, universe), use.ranks = FALSE)
          camera_direction <- as.character(cam$Direction[1])
          rows[[row_i]] <- data.frame(
            dataset = dataset,
            peyronie_comparison = pey_comp,
            mural_subtype = subtype,
            ed_comparison = ed_comp,
            ed_direction = ed_direction,
            n_signature_genes_in_universe = length(genes),
            camera_direction = camera_direction,
            camera_p_value = cam$PValue[1],
            validation_mode = expected_mode(pey_comp, ed_direction, camera_direction),
            genes = paste(sort(genes), collapse = ","),
            stringsAsFactors = FALSE
          )
          row_i <- row_i + 1
        }
      }
    }
  }
}

out_dir <- file.path(project_root, "results/tables/refined_mural_validation")
if (length(rows) == 0) {
  results <- data.frame()
} else {
  results <- do.call(rbind, rows)
  results$camera_FDR <- p.adjust(results$camera_p_value, method = "BH")
  results <- results[order(results$camera_FDR, results$camera_p_value), ]
}

write_tsv(results, file.path(out_dir, "refined_mural_cameraPR_results.tsv"))
cat("Wrote", nrow(results), "refined mural cameraPR rows\n")
