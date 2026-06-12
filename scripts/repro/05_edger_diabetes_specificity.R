#!/usr/bin/env Rscript
# =============================================================================
# 05_edger_diabetes_specificity.R
#
# edgeR-QL exact-engine reproduction of donor-level Pericyte_like pseudobulk DE
# for the DIABETES-SPECIFICITY contrast on GSE206528:
#
#   PRIMARY  : diabetic_ED  vs  non_diabetic_ED   (NEW; the diabetes-associated axis)
#   CONTROL  : diabetic_ED  vs  normal            (reference edgeR contrast)
#
# Imitates scripts/r/gse206528_mural_subtype_edger.R exactly:
#   DGEList -> filterByExpr -> normLibSizes -> estimateDisp -> glmQLFit(robust=TRUE)
#   -> glmQLFTest(coef=2)
#
# Donor is the unit (pericyte_like pseudobulk per donor). No cell-level inference.
# We FORCE-KEEP the 8 core genes through filterByExpr so we always report them
# (filterByExpr is computed on the full matrix, as in the reference script; we
#  then ensure the 8 reporting genes survive so none is silently dropped).
#
# Data provenance: GSE206528 pericyte-like pseudobulk (donor x gene), derived by
# our own pipeline from GEO accession GSE206528. Counts are integer pseudobulk
# sums. No external intermediate files are reused.
#
# Output: results/tables/repro/05_edger_diabetes_specificity.csv
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
})

root <- "/work"  # bind-mounted repo root inside the Docker image

counts_file <- file.path(
  root, "data", "processed", "GSE206528",
  "pseudobulk_mural_subtypes",
  "GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
)
meta_file <- file.path(
  root, "data", "processed", "GSE206528",
  "pseudobulk_mural_subtypes",
  "GSE206528_mural_subtype_pseudobulk_metadata.tsv"
)
out_dir <- file.path(root, "results", "tables", "repro")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "05_edger_diabetes_specificity.csv")

CORE_GENES <- c("IL6", "ADAMTS4", "PLA2G2A", "CYR61", "PI15",
                "SH3BP5", "ERAP2", "RGS16")
SHARED_ED  <- c("IL6", "ADAMTS4", "PLA2G2A", "CYR61", "PI15")
DIAB_AXIS  <- c("SH3BP5", "ERAP2")

MURAL_SUBTYPE <- "Pericyte_like"
min_group_n <- 2L  # matches reference script

# ---------------------------------------------------------------------------
# Load counts + metadata (same loading idiom as the reference script)
# ---------------------------------------------------------------------------
counts_dt <- fread(cmd = paste("gzip -dc", shQuote(counts_file)))
gene_symbols <- counts_dt$gene_symbol
counts <- as.matrix(counts_dt[, setdiff(names(counts_dt), "gene_symbol"), with = FALSE])
rownames(counts) <- gene_symbols
storage.mode(counts) <- "integer"

meta <- as.data.frame(fread(meta_file))
rownames(meta) <- meta$pseudobulk_id
meta <- meta[colnames(counts), , drop = FALSE]

# Restrict to the pericyte-like pseudobulks
peri_meta <- meta[meta$mural_subtype == MURAL_SUBTYPE, , drop = FALSE]
cat("Pericyte_like donor pseudobulks (donor x ed_subgroup x n_cells):\n")
print(peri_meta[, c("donor_id", "ed_subgroup", "diabetes_status", "n_cells")])
cat("\n")

# ---------------------------------------------------------------------------
# edgeR-QL runner (identical pipeline to gse206528_mural_subtype_edger.R)
# ---------------------------------------------------------------------------
run_edgeR <- function(comparison_name, selected_meta, group_values, group_levels,
                      report_genes) {
  selected_meta$comparison_group <- factor(group_values[rownames(selected_meta)],
                                           levels = group_levels)
  selected_meta <- selected_meta[!is.na(selected_meta$comparison_group), , drop = FALSE]
  selected_counts <- counts[, rownames(selected_meta), drop = FALSE]
  selected_meta$comparison_group <- droplevels(selected_meta$comparison_group)

  group_table <- table(selected_meta$comparison_group)
  stopifnot(length(group_table) == 2, all(group_table >= min_group_n))

  y <- DGEList(counts = selected_counts, samples = selected_meta,
               genes = data.frame(gene_symbol = rownames(selected_counts)))
  keep_genes <- filterByExpr(y, group = selected_meta$comparison_group)
  # force-keep the reporting genes so none is silently filtered out
  keep_genes[rownames(y) %in% report_genes] <- TRUE
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)

  design <- model.matrix(~ comparison_group, data = selected_meta)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2)               # coef 2 = second level vs reference
  res <- topTags(qlf, n = Inf, sort.by = "none")$table
  res$gene_symbol <- rownames(res)
  res$comparison <- comparison_name
  res$n_genes_tested <- nrow(res)
  res$n_samples <- nrow(selected_meta)
  res$group_summary <- paste(names(group_table), as.integer(group_table),
                             sep = "=", collapse = ";")
  # reference level is the first level; coef2 = (second level) vs (first level)
  res$reference_level <- group_levels[1]
  res$test_level <- group_levels[2]
  res
}

# ---------------------------------------------------------------------------
# PRIMARY contrast: diabetic_ED vs non_diabetic_ED
#   levels = c("non_diabetic_ED","diabetic_ED") so coef2 = diabetic_ED vs non_diab
#   -> positive logFC = UP in diabetic_ED relative to non-diabetic ED
# ---------------------------------------------------------------------------
prim_meta <- peri_meta[peri_meta$ed_subgroup %in% c("non_diabetic_ED", "diabetic_ED"), , drop = FALSE]
prim_group <- prim_meta$ed_subgroup
names(prim_group) <- rownames(prim_meta)
res_primary <- run_edgeR(
  "diabetic_ED_vs_non_diabetic_ED",
  prim_meta, prim_group,
  c("non_diabetic_ED", "diabetic_ED"),
  CORE_GENES
)

# ---------------------------------------------------------------------------
# CONTROL contrast: diabetic_ED vs normal (reference edgeR contrast)
# ---------------------------------------------------------------------------
ctrl_meta <- peri_meta[peri_meta$ed_subgroup %in% c("normal", "diabetic_ED"), , drop = FALSE]
ctrl_group <- ctrl_meta$ed_subgroup
names(ctrl_group) <- rownames(ctrl_meta)
res_control <- run_edgeR(
  "diabetic_ED_vs_normal",
  ctrl_meta, ctrl_group,
  c("normal", "diabetic_ED"),
  CORE_GENES
)

# ---------------------------------------------------------------------------
# Assemble: core genes only, both contrasts, long format
# ---------------------------------------------------------------------------
cols <- c("comparison", "gene_symbol", "logFC", "logCPM", "F", "PValue", "FDR",
          "reference_level", "test_level", "n_samples", "n_genes_tested",
          "group_summary")
prim_core <- as.data.table(res_primary)[gene_symbol %in% CORE_GENES, ..cols]
ctrl_core <- as.data.table(res_control)[gene_symbol %in% CORE_GENES, ..cols]

# annotate gene class
class_map <- c(setNames(rep("shared-ED", length(SHARED_ED)), SHARED_ED),
               setNames(rep("diabetes-associated-axis", length(DIAB_AXIS)), DIAB_AXIS),
               RGS16 = "stress/IEG")
prim_core[, gene_class := class_map[gene_symbol]]
ctrl_core[, gene_class := class_map[gene_symbol]]

# order genes: shared-ED first, then diabetes axis, then RGS16
gene_order <- c(SHARED_ED, DIAB_AXIS, "RGS16")
prim_core[, gene_symbol := factor(gene_symbol, levels = gene_order)]
ctrl_core[, gene_symbol := factor(gene_symbol, levels = gene_order)]
setorder(prim_core, gene_symbol)
setorder(ctrl_core, gene_symbol)
prim_core[, gene_symbol := as.character(gene_symbol)]
ctrl_core[, gene_symbol := as.character(gene_symbol)]

out_long <- rbindlist(list(prim_core, ctrl_core))
setcolorder(out_long, c("comparison", "gene_symbol", "gene_class",
                        "logFC", "logCPM", "F", "PValue", "FDR",
                        "reference_level", "test_level"))

fwrite(out_long, out_csv)

# ---------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------
cat("================ edgeR-QL results: 8 core genes ================\n\n")
cat(">>> PRIMARY: diabetic_ED vs non_diabetic_ED",
    "(positive logFC = UP in diabetic_ED)\n")
cat("    n_samples =", res_primary$n_samples[1],
    " group_summary =", res_primary$group_summary[1],
    " genes_tested =", res_primary$n_genes_tested[1], "\n")
print(prim_core[, .(gene_symbol, gene_class,
                    logFC = round(logFC, 3),
                    logCPM = round(logCPM, 2),
                    F = round(F, 2),
                    PValue = signif(PValue, 3),
                    FDR = signif(FDR, 3))])

cat("\n>>> CONTROL: diabetic_ED vs normal (reference edgeR contrast)\n")
cat("    n_samples =", res_control$n_samples[1],
    " group_summary =", res_control$group_summary[1],
    " genes_tested =", res_control$n_genes_tested[1], "\n")
print(ctrl_core[, .(gene_symbol, gene_class,
                    logFC = round(logFC, 3),
                    logCPM = round(logCPM, 2),
                    F = round(F, 2),
                    PValue = signif(PValue, 3),
                    FDR = signif(FDR, 3))])

cat("\n---- diabetic-versus-non-diabetic-ED contrast verdict (diabetic_ED vs non_diabetic_ED) ----\n")
for (g in c(SHARED_ED, DIAB_AXIS)) {
  r <- prim_core[gene_symbol == g]
  if (nrow(r) == 0) next
  cat(sprintf("  %-9s class=%-22s logFC=%+.3f  FDR=%.3g  rawP=%.3g\n",
              g, r$gene_class, r$logFC, r$FDR, r$PValue))
}
cat("\nExpectation: SH3BP5 logFC>0, ERAP2 logFC<0 are the largest-magnitude /",
    "lowest-p genes; the 5 shared-ED genes ~0 and non-significant.\n")

cat("\nSaved ->", out_csv, "\n")
