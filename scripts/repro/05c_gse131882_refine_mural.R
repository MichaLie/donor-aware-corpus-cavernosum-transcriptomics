#!/usr/bin/env Rscript
# Stage 2: the global k=20 clustering did NOT yield a clean RGS5+/NOTCH3+ pericyte cluster
# (cluster 9 = PDGFRB-high but RGS5~0, NOTCH3~0, and Podocyte-marker high -> mesangial/podocyte mix).
# Refine: rebuild from raw, higher-resolution clustering, and a marker-based pericyte/mesangial-mural
# gate so the pseudobulk is on genuinely mural cells. Report what we find honestly either way.
suppressMessages({
  library(Matrix); library(SingleCellExperiment); library(scran); library(scater)
  library(igraph); library(org.Hs.eg.db); library(AnnotationDbi)
})
set.seed(0)
root <- Sys.getenv("PROJECT_ROOT", unset = "/work")
RAW  <- file.path(root, "data", "raw", "GSE131882")
PROC <- file.path(root, "data", "processed", "GSE131882")
samples <- data.frame(
  file  = c("GSM3823939_control.s1.dgecounts.rds","GSM3823940_control.s2.dgecounts.rds",
            "GSM3823941_control.s3.dgecounts.rds","GSM3823942_diabetes.s1.dgecounts.rds",
            "GSM3823943_diabetes.s2.dgecounts.rds","GSM3823944_diabetes.s3.dgecounts.rds"),
  donor = c("control_s1","control_s2","control_s3","diabetes_s1","diabetes_s2","diabetes_s3"),
  cond  = c("control","control","control","diabetes","diabetes","diabetes"),
  stringsAsFactors=FALSE)
MIN_UMI <- 500
mats <- lapply(seq_len(nrow(samples)), function(i){
  m <- readRDS(file.path(RAW, samples$file[i]))$umicount$exon$all
  m <- m[, Matrix::colSums(m) >= MIN_UMI, drop=FALSE]
  colnames(m) <- paste0(samples$donor[i], ":", colnames(m)); m
})
common <- Reduce(intersect, lapply(mats, rownames))
counts <- do.call(cbind, lapply(mats, function(m) m[common, , drop=FALSE]))
cell_donor <- sub(":.*$","",colnames(counts))
cd <- samples[match(cell_donor, samples$donor), c("donor","cond")]
sce <- SingleCellExperiment(assays=list(counts=counts),
                            colData=DataFrame(cd, row.names=colnames(counts)))
sym <- mapIds(org.Hs.eg.db, keys=rownames(sce), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
rowData(sce)$symbol <- ifelse(is.na(sym), rownames(sce), sym)
cls <- tryCatch(quickCluster(sce, min.size=100), error=function(e) NULL)
sce <- computeSumFactors(sce, clusters=cls, min.mean=0.1)
sce <- sce[, sizeFactors(sce) > 0]                 # drop non-positive size factors
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce); hvg <- getTopHVGs(dec, n=2000)
sce <- fixedPCA(sce, subset.row=hvg, rank=30)
g <- buildSNNGraph(sce, k=10, use.dimred="PCA")    # higher resolution (smaller k)
cl <- igraph::cluster_louvain(g, resolution=2)$membership
colLabels(sce) <- factor(cl)
symv <- rowData(sce)$symbol
gm <- function(gs){ gs<-gs[gs%in%symv]; if(!length(gs)) return(NULL)
  idx<-match(gs,symv); m<-logcounts(sce)[idx,,drop=FALSE]; rownames(m)<-gs
  res<-sapply(split(seq_len(ncol(m)),colLabels(sce)),function(j) Matrix::rowMeans(m[,j,drop=FALSE]))
  if(is.null(dim(res))){ res<-matrix(res,nrow=length(gs),dimnames=list(gs,names(res))) }
  res }
mural  <- gm(c("RGS5","PDGFRB","NOTCH3","GATA3"))      # GATA3 to exclude CD-cells; mesangial=PDGFRB+GATA3-
podo   <- gm(c("NPHS1","NPHS2","PODXL"))
endo   <- gm(c("PECAM1","FLT1","EMCN"))
immune <- gm(c("PTPRC","LYZ","CD68"))
tub    <- gm(c("LRP2","SLC34A1","UMOD","SLC12A1","SLC12A3","AQP2"))
saveRDS(sce, file.path(PROC,"GSE131882_clustered_sce.rds"))
cat("=== refined clusters: mural markers ===\n"); print(round(mural,2))
cat("\npodocyte (mean of NPHS1/NPHS2/PODXL):\n"); print(round(colMeans(podo),2))
cat("endothelial:\n"); print(round(colMeans(endo),2))
cat("immune:\n"); print(round(colMeans(immune),2))
cat("tubular(max):\n"); print(round(apply(tub,2,max),2))
cat("cluster sizes:\n"); print(table(colLabels(sce)))
# pericyte/mesangial-mural gate: PDGFRB-high AND (RGS5 or NOTCH3) present AND podocyte/endo/immune low
pdgfrb <- mural["PDGFRB",]; rgs5 <- mural["RGS5",]; notch3 <- mural["NOTCH3",]
podsc <- colMeans(podo); endosc <- colMeans(endo); immsc <- colMeans(immune); tubsc <- apply(tub,2,max)
score <- pdgfrb + pmax(rgs5,notch3)
ok <- pdgfrb > 0.15 & (rgs5 > 0.05 | notch3 > 0.05) & podsc < 0.5 & endosc < 0.1 & immsc < 0.1 & tubsc < 0.5
cat("\nmural gate score (PDGFRB + max(RGS5,NOTCH3)):\n"); print(round(sort(score,decreasing=TRUE),2))
cat("clusters passing strict mural gate:", paste(names(which(ok)),collapse=", "), "\n")
mural_clusters <- if (any(ok)) names(which(ok)) else names(which.max(score))
cat(">>> mural cluster(s):", paste(mural_clusters,collapse=", "), "\n")
peri_cells <- colLabels(sce) %in% mural_clusters
cat("mural cells per donor:\n"); print(table(sce$donor[peri_cells]))
# pseudobulk by symbol
pbc <- sapply(split(which(peri_cells), sce$donor[peri_cells]), function(j)
  Matrix::rowSums(counts(sce)[, j, drop=FALSE]))
pbdf <- as.data.frame(as.matrix(pbc)); pbdf$symbol <- symv
pbsym <- rowsum(pbdf[,colnames(pbc),drop=FALSE], group=pbdf$symbol)
write.csv(pbsym, file.path(PROC,"GSE131882_mural_pseudobulk_counts_bySymbol_refined.csv"))
cat("wrote refined mural pseudobulk; donors:", paste(colnames(pbc),collapse=", "),"\n")
cat("counts/donor:\n"); print(round(colSums(pbc)))
# per-cell activation log-norm on refined mural cells
act <- c("IL6","ADAMTS4","PLA2G2A","CCN1","CYR61","PI15","SH3BP5","ERAP2","RGS16","NR4A1","NR4A2","NR4A3","MT1A")
present <- act[act %in% symv]
pc <- as.data.frame(t(as.matrix(logcounts(sce)[match(present,symv), peri_cells, drop=FALSE])))
colnames(pc) <- present; pc$donor <- sce$donor[peri_cells]; pc$cond <- sce$cond[peri_cells]
write.csv(pc, file.path(PROC,"GSE131882_mural_percell_activation_lognorm_refined.csv"), row.names=FALSE)
cat("refined mural cells total:", sum(peri_cells), "\n")
cat("activation genes present:", paste(present,collapse=", "), " | absent:", paste(setdiff(act,present),collapse=", "),"\n")
cat("STAGE2 DONE\n")
