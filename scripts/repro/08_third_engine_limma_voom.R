#!/usr/bin/env Rscript
# Third DE engine (limma-voom; the framework underlying muscat/dreamlet) on the pericyte-like
# donor pseudobulk, to confirm SH3BP5/ERAP2 diabetic-versus-non-diabetic-ED contrast across a third statistical method.
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)})
root <- "/work"
cnt <- fread(cmd=paste("gzip -dc", file.path(root,"data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_counts.tsv.gz")))
genes <- cnt$gene_symbol; m <- as.matrix(cnt[,-1]); rownames(m) <- genes
meta <- as.data.frame(fread(file.path(root,"data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_metadata.tsv")))
meta <- meta[meta$mural_subtype=="Pericyte_like",]; rownames(meta) <- meta$pseudobulk_id
m <- m[, rownames(meta)]; storage.mode(m) <- "integer"
grp <- c(diabetic_ED="diab",non_diabetic_ED="nondiab",normal="normal")[meta$ed_subgroup]
meta$cond <- factor(grp, levels=c("normal","nondiab","diab"))
core <- c("IL6","ADAMTS4","PLA2G2A","CYR61","PI15","SH3BP5","ERAP2","RGS16")

y <- DGEList(counts=m, samples=meta)
keep <- filterByExpr(y, group=meta$cond) | rownames(y) %in% core
y <- y[keep,,keep.lib.sizes=FALSE]; y <- normLibSizes(y)
design <- model.matrix(~0+cond, data=meta); colnames(design) <- levels(meta$cond)
v <- voom(y, design)
fit <- lmFit(v, design)
cm <- makeContrasts(diab_vs_nondiab=diab-nondiab, diab_vs_normal=diab-normal, levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm))
cat("=== limma-voom, Pericyte_like donor pseudobulk ===\n")
for (ct in colnames(cm)) {
  tt <- topTable(fit2, coef=ct, number=Inf, sort.by="none")
  cat(sprintf("\n--- %s ---\n", ct))
  sub <- tt[core,c("logFC","P.Value","adj.P.Val")]
  print(round(sub,4))
}
# write the diab-vs-nondiab core result
tt <- topTable(fit2, coef="diab_vs_nondiab", number=Inf, sort.by="none")
out <- tt[core,c("logFC","P.Value","adj.P.Val")]; out$gene <- rownames(out)
fwrite(out, file.path(root,"results/tables/repro/08_limma_voom_diab_vs_nondiab.csv"))
cat("\nsaved -> results/tables/repro/08_limma_voom_diab_vs_nondiab.csv\n")
