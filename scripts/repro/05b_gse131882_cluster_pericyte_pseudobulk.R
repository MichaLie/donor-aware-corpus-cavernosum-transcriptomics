#!/usr/bin/env Rscript
# GSE131882 (human kidney snRNA-seq, diabetic nephropathy vs control; Wu et al. 2019, Humphreys lab).
# Re-fetched from GEO per-GSM FTP as zUMIs/dropEst dgecounts.rds (umicount$exon$all = Ensembl x cells).
# 3 control donors (s1-s3), 3 diabetes donors (s1-s3) -> matches the 3-vs-3 design in the task.
#
# Goal of this stage:
#  - Load all 6 samples' UMI matrices, QC, merge to one SingleCellExperiment.
#  - Map Ensembl -> symbol (org.Hs.eg.db).
#  - Normalize + cluster (scran/igraph Louvain) and annotate clusters by canonical kidney markers.
#  - Identify the pericyte/mural cluster (RGS5, PDGFRB, NOTCH3).
#  - Export: (a) per-donor pericyte PSEUDOBULK raw counts (gene x 6 donors),
#            (b) cluster marker means table, (c) per-cell pericyte expression for the
#               activation genes (descriptive effect sizes only -- inference is donor-level downstream).
suppressMessages({
  library(Matrix); library(SingleCellExperiment); library(scran); library(scater)
  library(igraph); library(org.Hs.eg.db); library(AnnotationDbi)
})
set.seed(0)
ROOT <- "/work"
RAW  <- file.path(ROOT, "data/raw/GSE131882")
PROC <- file.path(ROOT, "data/processed/GSE131882"); dir.create(PROC, showWarnings=FALSE, recursive=TRUE)
OUT  <- file.path(ROOT, "results/tables/repro");      dir.create(OUT,  showWarnings=FALSE, recursive=TRUE)

samples <- data.frame(
  gsm   = c("GSM3823939","GSM3823940","GSM3823941","GSM3823942","GSM3823943","GSM3823944"),
  file  = c("GSM3823939_control.s1.dgecounts.rds","GSM3823940_control.s2.dgecounts.rds",
            "GSM3823941_control.s3.dgecounts.rds","GSM3823942_diabetes.s1.dgecounts.rds",
            "GSM3823943_diabetes.s2.dgecounts.rds","GSM3823944_diabetes.s3.dgecounts.rds"),
  donor = c("control_s1","control_s2","control_s3","diabetes_s1","diabetes_s2","diabetes_s3"),
  cond  = c("control","control","control","diabetes","diabetes","diabetes"),
  stringsAsFactors=FALSE)

MIN_UMI <- 500   # nuclei QC: keep cells with >=500 UMI (these are snRNA-seq nuclei)
mats <- list()
for (i in seq_len(nrow(samples))) {
  f <- file.path(RAW, samples$file[i])
  obj <- readRDS(f)
  m <- obj$umicount$exon$all              # genes(Ensembl) x cells, dgCMatrix
  keep <- Matrix::colSums(m) >= MIN_UMI
  m <- m[, keep, drop=FALSE]
  colnames(m) <- paste0(samples$donor[i], ":", colnames(m))
  mats[[i]] <- m
  cat(sprintf("%s (%s): %d genes x %d cells passing >=%d UMI (median %d UMI/cell)\n",
              samples$gsm[i], samples$donor[i], nrow(m), ncol(m), MIN_UMI,
              as.integer(median(Matrix::colSums(m)))))
  rm(obj); gc()
}
# union of genes (all should share the same zUMIs reference, but be safe)
all_genes <- sort(unique(unlist(lapply(mats, rownames))))
align <- function(m) {
  full <- Matrix(0, nrow=length(all_genes), ncol=ncol(m), sparse=TRUE,
                 dimnames=list(all_genes, colnames(m)))
  full[rownames(m), ] <- m
  full
}
same_genes <- all(sapply(mats, function(m) identical(rownames(m), rownames(mats[[1]]))))
cat("all samples share identical gene order:", same_genes, "\n")
counts <- if (same_genes) do.call(cbind, mats) else do.call(cbind, lapply(mats, align))
cat("merged counts:", nrow(counts), "genes x", ncol(counts), "cells\n")

# cell metadata
cell_donor <- sub(":.*$", "", colnames(counts))
cd <- samples[match(cell_donor, samples$donor), c("donor","cond","gsm")]
rownames(cd) <- colnames(counts)
sce <- SingleCellExperiment(assays=list(counts=counts), colData=DataFrame(cd))

# Ensembl -> symbol
ens <- rownames(sce)
sym <- mapIds(org.Hs.eg.db, keys=ens, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
rowData(sce)$ensembl <- ens
rowData(sce)$symbol  <- ifelse(is.na(sym), ens, sym)
cat("genes with a symbol:", sum(!is.na(sym)), "/", length(ens), "\n")

# Normalize (library-size factors + log) -- scran pooling clusters then deconvolution
clust0 <- tryCatch(quickCluster(sce, min.size=100), error=function(e) NULL)
sce <- computeSumFactors(sce, clusters=clust0, min.mean=0.1)
sce <- logNormCounts(sce)

# HVGs, PCA, graph clustering (Louvain)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n=2000)
sce <- fixedPCA(sce, subset.row=hvg, rank=30)
g <- buildSNNGraph(sce, k=20, use.dimred="PCA")
cl <- igraph::cluster_louvain(g)$membership
colLabels(sce) <- factor(cl)
cat("clusters:", paste(sort(unique(cl)), collapse=","), " (n=", length(unique(cl)), ")\n")

# marker means per cluster (mean log-norm) for annotation
markers <- list(
  Pericyte_Mural = c("RGS5","PDGFRB","NOTCH3","ACTA2","TAGLN","MYH11","PDGFRA"),
  Endothelial    = c("PECAM1","FLT1","EMCN","KDR"),
  Podocyte       = c("NPHS1","NPHS2","PODXL"),
  Proximal_Tub   = c("LRP2","CUBN","SLC34A1","SLC5A12"),
  LoopHenle      = c("UMOD","SLC12A1"),
  Distal_Coll    = c("SLC12A3","AQP2","GATA3"),
  Immune         = c("PTPRC","CD3E","LYZ","CD68"),
  Fibroblast     = c("COL1A1","COL1A2","DCN")
)
ln <- logcounts(sce)
sym_vec <- rowData(sce)$symbol
mean_by_cluster <- function(genes) {
  genes <- genes[genes %in% sym_vec]
  if (!length(genes)) return(setNames(rep(NA, nlevels(colLabels(sce))), levels(colLabels(sce))))
  idx <- match(genes, sym_vec)
  sub <- ln[idx, , drop=FALSE]
  sapply(split(seq_len(ncol(sub)), colLabels(sce)), function(j) mean(Matrix::rowMeans(sub[, j, drop=FALSE])))
}
mk_tab <- do.call(rbind, lapply(names(markers), function(nm) mean_by_cluster(markers[[nm]])))
rownames(mk_tab) <- names(markers)
cat("\n=== mean log-norm marker score per cluster (rows=lineage) ===\n")
print(round(mk_tab, 2))

# also per-cluster individual mural marker means + cluster sizes + donor composition
mural_genes <- c("RGS5","PDGFRB","NOTCH3","ACTA2","TAGLN")
mural_idx <- match(mural_genes[mural_genes %in% sym_vec], sym_vec)
mural_each <- sapply(split(seq_len(ncol(ln)), colLabels(sce)), function(j)
  Matrix::rowMeans(ln[mural_idx, j, drop=FALSE]))
rownames(mural_each) <- mural_genes[mural_genes %in% sym_vec]
cat("\n=== individual mural marker mean log-norm per cluster ===\n")
print(round(mural_each, 2))
csize <- table(colLabels(sce))
cat("\ncluster sizes:\n"); print(csize)

# pick the pericyte/mural cluster: highest combined RGS5+PDGFRB+NOTCH3, and clearly above endothelial/epithelial
mural_score <- mk_tab["Pericyte_Mural", ]
# require RGS5 and PDGFRB both expressed (mural identity), not an SMC-only or contamination cluster
rgs5  <- mural_each["RGS5", ]; pdgfrb <- mural_each["PDGFRB", ]; notch3 <- mural_each["NOTCH3", ]
peri_combo <- rgs5 + pdgfrb + notch3
cat("\nRGS5+PDGFRB+NOTCH3 combined per cluster:\n"); print(round(sort(peri_combo, decreasing=TRUE), 2))
peri_cluster <- names(which.max(peri_combo))
cat("\n>>> selected pericyte/mural cluster:", peri_cluster,
    " (n=", as.integer(csize[peri_cluster]), "cells)\n")

# donor composition of the pericyte cluster
peri_cells <- colLabels(sce) == peri_cluster
cat("\npericyte-cluster cells per donor:\n")
print(table(sce$donor[peri_cells]))

# ---- PSEUDOBULK: sum raw counts over pericyte cells, per donor ----
peri_sce <- sce[, peri_cells]
pb <- sapply(split(seq_len(ncol(peri_sce)), peri_sce$donor), function(j)
  Matrix::rowSums(counts(peri_sce)[, j, drop=FALSE]))
# collapse to symbol (sum Ensembl ids mapping to same symbol)
pb_df <- as.data.frame(as.matrix(pb))
pb_df$symbol <- sym_vec
pb_sym <- rowsum(pb_df[, samples$donor[samples$donor %in% colnames(pb)], drop=FALSE],
                 group=pb_df$symbol)
write.csv(pb_sym, file.path(PROC, "GSE131882_pericyte_pseudobulk_counts_bySymbol.csv"))
write.csv(as.matrix(pb), file.path(PROC, "GSE131882_pericyte_pseudobulk_counts_byEnsembl.csv"))
cat("\nwrote pericyte pseudobulk (gene x donor) to data/processed/GSE131882/\n")
cat("pseudobulk donors:", paste(colnames(pb), collapse=", "), "\n")
cat("pseudobulk total counts per donor:\n"); print(round(colSums(pb)))

# ---- per-cell pericyte expression of activation genes (descriptive only) ----
act_genes <- c("IL6","ADAMTS4","PLA2G2A","CYR61","CCN1","PI15","SH3BP5","ERAP2","RGS16",
               "NR4A1","NR4A2","NR4A3","MT1A")
present <- act_genes[act_genes %in% sym_vec]
cat("\nactivation/stress genes present in data:", paste(present, collapse=", "), "\n")
cat("activation/stress genes ABSENT:", paste(setdiff(act_genes, present), collapse=", "), "\n")
idx <- match(present, sym_vec)
peri_ln <- ln[idx, peri_cells, drop=FALSE]
percell <- as.data.frame(t(as.matrix(peri_ln)))
colnames(percell) <- present
percell$donor <- peri_sce$donor
percell$cond  <- peri_sce$cond
write.csv(percell, file.path(PROC, "GSE131882_pericyte_percell_activation_lognorm.csv"), row.names=FALSE)
cat("wrote per-cell pericyte activation log-norm matrix (",
    nrow(percell), "pericyte cells )\n")

saveRDS(list(samples=samples, peri_cluster=peri_cluster,
             cluster_sizes=as.matrix(csize), mk_tab=mk_tab, mural_each=mural_each),
        file.path(PROC, "GSE131882_cluster_summary.rds"))
cat("\nSTAGE1 DONE\n")
