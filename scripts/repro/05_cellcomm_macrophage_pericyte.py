#!/usr/bin/env python3
"""
05 — Macrophage -> pericyte inflammatory axis, tested on the archived annotations.

Analysis question: are macrophage-derived TNF/OSM patterns consistent with
pericyte NF-kB/STAT3 activity and the
pericyte activation core. That claim is pure INFERENCE (CellChat ligand-receptor +
TF-activity). We CANNOT measure signaling; we can only ask whether the *expression*
patterns are consistent with it. This script therefore makes two descriptive,
inference-flavoured checks and is explicit that neither is a measurement of signaling:

  (1) Are macrophage ligands (TNF, IL1B, OSM) higher in DIABETIC macrophages than in
      reference macrophages?  Donor-level pseudobulk fraction-of-macrophages-expressing
      and mean log-norm, in GSE206528 (8 donors: 2 diab / 3 nondiab-ED / 3 normal).
      GSE259348 (single diabetic donor) is positioned descriptively against that range.

  (2) Within PERICYTES, does the pericyte-activation "core score" (mean of the shared-ED
      genes IL6/ADAMTS4/PLA2G2A/CYR61/PI15) co-vary with *inferred* TF-target-module
      activity for STAT3 and NF-kB (RELA/NFKB1 targets), versus a negative-control TF
      module (GATA1, an erythroid TF with no expected role here)?
      INFERENCE = TF activity is proxied by the average expression of curated target
      genes (an AUCell/module-score surrogate), NOT a measured TF activity.

DISCIPLINE (anti-pseudoreplication):
  - Donor is the unit. The ONLY inferential numbers are donor-level: per-donor pericyte
    means of (core score, TF modules), then Pearson r across the 8 donors, plus a
    Spearman as a rank check. Cell-level correlations are reported ONLY as descriptive
    effect sizes and explicitly flagged as non-inferential (pseudoreplication).
  - GATA1 negative control must be ~0 / non-significant for the test to be meaningful.

Data re-derivation: GSE206528 cells were re-fetched from GEO and annotated by us
(major_cell_type for Macrophage; mural_subtype Pericyte_like for pericytes).
GSE259348 (GSM8114573, independent diabetic donor) was re-fetched from GEO FTP and
annotated by us in script 04; here pericytes/macrophages are taken from its leiden
clusters by marker expression (RGS5/NOTCH3-high mural = pericyte-like; CD68/CD163-high
= macrophage).
"""
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
import scipy.sparse as sp
from scipy.stats import pearsonr, spearmanr

ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
OUT=f"{ROOT}/results/tables/repro"; os.makedirs(OUT, exist_ok=True)
np.random.seed(0)

# ---------- gene sets ----------
LIGANDS=["TNF","IL1B","OSM"]
CORE=["IL6","ADAMTS4","PLA2G2A","CYR61","PI15"]   # shared-ED pericyte activation core
# TF target modules = AUCell/module-score surrogate for INFERRED TF activity (NOT measured).
# Curated, literature-standard targets; core genes deliberately excluded to avoid circularity.
STAT3_T=["SOCS3","BCL3","MCL1","BCL2L1","VEGFA","MMP9","CCND1","MYC","JUNB","OSMR",
         "IL6R","BATF","HIF1A","TIMP1","SAA1","SAA2","PIM1"]
NFKB_T =["NFKBIA","RELB","TNFAIP3","BIRC3","CXCL1","CXCL2","CXCL3","CXCL8","CCL2","CCL20",
         "ICAM1","VCAM1","SOD2","PTGS2","NFKB2","CXCL5","SELE","TRAF1","IER3"]  # RELA/NFKB1 targets, ligands TNF/IL1B excluded
GATA1_T=["HBB","HBA1","HBA2","ALAS2","SLC4A1","EPB42","SPTA1","AHSP","TAL1","RHAG",
         "EPOR","ANK1","TFRC"]  # negative-control erythroid module
MODULES={"STAT3_targets":STAT3_T,"NFKB_targets":NFKB_T,"GATA1_neg_targets":GATA1_T}

def dense(M):
    return M.toarray() if sp.issparse(M) else np.asarray(M)

def module_score(adata, genes, name):
    g=[x for x in genes if x in adata.var_names]
    sc.tl.score_genes(adata, g, score_name=name, ctrl_size=50, random_state=0, use_raw=False)
    return g

def mean_expr_frac(adata, genes):
    """per-cell mean log-norm of a gene list and fraction expressing each gene."""
    g=[x for x in genes if x in adata.var_names]
    X=dense(adata[:,g].X)
    return pd.DataFrame(X, columns=g, index=adata.obs_names)

rows=[]   # tidy long output

# =========================================================================
# PART A — GSE206528 (8 donors)
# =========================================================================
print("="*70); print("GSE206528"); print("="*70)
a=ad.read_h5ad(f"{ROOT}/data/processed/GSE206528/GSE206528_scanpy_preprocessed.h5ad")
major=pd.read_csv(f"{ROOT}/results/tables/GSE206528_cell_major_celltype_annotations.tsv.gz",sep="\t",index_col=0)
mural=pd.read_csv(f"{ROOT}/results/tables/GSE206528_cell_mural_subtype_annotations.tsv.gz",sep="\t",index_col=0)
a.obs["major_ct"]=major["major_cell_type"].reindex(a.obs_names).values
a.obs["mural_subtype"]=mural["mural_subtype"].reindex(a.obs_names).values
# X is log-norm (verified max ~6.8); score modules on all cells once
for nm,gs in MODULES.items():
    used=module_score(a, gs, nm)
    print(f"  module {nm}: {len(used)} genes used")
score_core=module_score(a, CORE, "core_score")

grp_map={"diabetic_ED":"diabetic","non_diabetic_ED":"nondiab_ED","normal":"normal"}
a.obs["grp"]=a.obs["ed_subgroup"].map(grp_map)

# ---- (1) macrophage ligands, donor-level ----
mac=a[a.obs.major_ct=="Macrophage"].copy()
print(f"\nMacrophages: {mac.n_obs}  by donor:")
print(mac.obs.groupby(['grp','donor_id'],observed=True).size().to_string())
ligX=mean_expr_frac(mac, LIGANDS)
ligX["donor_id"]=mac.obs["donor_id"].values; ligX["grp"]=mac.obs["grp"].values
# per-donor mean log-norm and fraction expressing
don_mean=ligX.groupby("donor_id",observed=True)[LIGANDS].mean()
don_frac=ligX.groupby("donor_id",observed=True)[LIGANDS].apply(lambda d:(d>0).mean())
don_grp=ligX.groupby("donor_id",observed=True)["grp"].first()
print("\n[GSE206528] macrophage ligand per-donor MEAN log-norm:")
print(pd.concat([don_grp.rename('grp'),don_mean.round(3)],axis=1).to_string())
print("\n[GSE206528] macrophage ligand per-donor FRACTION expressing:")
print(pd.concat([don_grp.rename('grp'),don_frac.round(3)],axis=1).to_string())
for lig in LIGANDS:
    for grp in ["diabetic","nondiab_ED","normal"]:
        dvals=don_mean.loc[don_grp[don_grp==grp].index, lig]
        rows.append(dict(dataset="GSE206528",analysis="ligand_donor_mean",metric=lig,
                         group=grp,n_donors=len(dvals),value=round(float(dvals.mean()),4),
                         sd=round(float(dvals.std(ddof=1)) if len(dvals)>1 else np.nan,4)))
# diabetic vs all-reference (nondiab+normal) effect, donor-level (descriptive; n too small for a real test)
for lig in LIGANDS:
    d=don_mean.loc[don_grp[don_grp=="diabetic"].index,lig]
    ref=don_mean.loc[don_grp[don_grp!="diabetic"].index,lig]
    rows.append(dict(dataset="GSE206528",analysis="ligand_diab_vs_ref_logfc",metric=lig,
                     group="diab_vs_ref",n_donors=int(len(d)+len(ref)),
                     value=round(float(d.mean()-ref.mean()),4),sd=np.nan))

# ---- (2) pericyte core vs TF modules, DONOR-level (inference) ----
peri=a[a.obs.mural_subtype=="Pericyte_like"].copy()
print(f"\nPericyte_like cells: {peri.n_obs}  by donor:")
pcount=peri.obs.groupby('donor_id',observed=True).size()
print(pcount.to_string())
pdf=peri.obs[["donor_id","grp","core_score"]+list(MODULES)].copy()
don=pdf.groupby("donor_id",observed=True)[["core_score"]+list(MODULES)].mean()
don_g=pdf.groupby("donor_id",observed=True)["grp"].first()
print("\n[GSE206528] per-donor pericyte module means (n=%d donors):"%len(don))
print(pd.concat([don_g.rename('grp'),don.round(3)],axis=1).to_string())

print("\n[GSE206528] DONOR-LEVEL correlation core_score vs TF module (INFERENTIAL UNIT):")
for mod in MODULES:
    r,p=pearsonr(don["core_score"],don[mod])
    rs,ps=spearmanr(don["core_score"],don[mod])
    print(f"   core ~ {mod:18s}  pearson r={r:+.3f} p={p:.3f}   spearman rho={rs:+.3f} p={ps:.3f}  (n={len(don)})")
    rows.append(dict(dataset="GSE206528",analysis="donor_core_vs_TF_pearson",metric=mod,
                     group="pericyte",n_donors=len(don),value=round(float(r),4),sd=round(float(p),4)))
    rows.append(dict(dataset="GSE206528",analysis="donor_core_vs_TF_spearman",metric=mod,
                     group="pericyte",n_donors=len(don),value=round(float(rs),4),sd=round(float(ps),4)))

# ---- cell-level effect SIZE only (descriptive; NOT inference) ----
print("\n[GSE206528] CELL-LEVEL core~TF Pearson (DESCRIPTIVE effect size ONLY, pseudoreplication if used as p):")
cell_core=peri.obs["core_score"].values
for mod in MODULES:
    r,_=pearsonr(cell_core, peri.obs[mod].values)
    print(f"   core ~ {mod:18s}  cell-level r={r:+.3f}  (n_cells={peri.n_obs})  [descriptive]")
    rows.append(dict(dataset="GSE206528",analysis="cell_core_vs_TF_pearson_DESCRIPTIVE",metric=mod,
                     group="pericyte",n_donors=peri.n_obs,value=round(float(r),4),sd=np.nan))

# =========================================================================
# PART B — GSE259348 (single diabetic donor; descriptive replication)
# =========================================================================
print("\n"+"="*70); print("GSE259348 (single diabetic donor — descriptive)"); print("="*70)
b=ad.read_h5ad(f"{ROOT}/data/processed/GSE259348/GSE259348_annotated_ours.h5ad")
# X is log-norm (normalize_total+log1p in script 04). Identify subsets from leiden by markers.
# pericyte-like mural: clusters with high RGS5/NOTCH3 (0,14); macrophage: CD68/CD163 high (6,7)
PERI_CL={"0","14"}; MAC_CL={"6","7"}
b.obs["subset"]=np.where(b.obs.leiden.astype(str).isin(PERI_CL),"Pericyte_like",
                 np.where(b.obs.leiden.astype(str).isin(MAC_CL),"Macrophage","other"))
print("subset counts:", b.obs.subset.value_counts().to_dict())
for nm,gs in MODULES.items():
    module_score(b, gs, nm)
module_score(b, CORE, "core_score")

bmac=b[b.obs.subset=="Macrophage"]
blig=mean_expr_frac(bmac, LIGANDS)
print("\n[GSE259348] macrophage ligand mean log-norm (single diab donor):")
mac_mean=blig.mean(); mac_frac=(blig>0).mean()
print("  mean:", mac_mean.round(3).to_dict())
print("  frac:", mac_frac.round(3).to_dict())
for lig in LIGANDS:
    rows.append(dict(dataset="GSE259348",analysis="ligand_donor_mean",metric=lig,group="diabetic",
                     n_donors=1,value=round(float(mac_mean[lig]),4),sd=round(float(mac_frac[lig]),4)))

bperi=b[b.obs.subset=="Pericyte_like"]
print(f"\n[GSE259348] Pericyte_like cells: {bperi.n_obs}")
print("[GSE259348] pericyte module means (single donor):")
mods=["core_score"]+list(MODULES)
print(bperi.obs[mods].mean().round(3).to_string())
# only a single donor: cell-level correlation is DESCRIPTIVE within-donor (NOT a disease claim)
print("\n[GSE259348] within-donor CELL-LEVEL core~TF Pearson (DESCRIPTIVE; single donor, NOT inference):")
ccore=bperi.obs["core_score"].values
for mod in MODULES:
    r,_=pearsonr(ccore, bperi.obs[mod].values)
    print(f"   core ~ {mod:18s}  cell-level r={r:+.3f} (n_cells={bperi.n_obs}) [descriptive]")
    rows.append(dict(dataset="GSE259348",analysis="cell_core_vs_TF_pearson_DESCRIPTIVE",metric=mod,
                     group="pericyte",n_donors=bperi.n_obs,value=round(float(r),4),sd=np.nan))

# =========================================================================
out=pd.DataFrame(rows)
out.to_csv(f"{OUT}/05_cellcomm.csv", index=False)
print("\nsaved ->", f"{OUT}/05_cellcomm.csv", "rows:", len(out))
