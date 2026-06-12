#!/usr/bin/env python3
"""
P4b — refine the GSE259348 mural compartment into pericyte-like vs SMC-like
(mirrors the GSE206528 mural refinement), then combine the
GSE259348 pericyte pseudobulk with GSE206528 pericyte pseudobulks and run a
donor-level DESeq2 ~dataset+condition for the 8 non-stress core genes.
Loads the saved annotated object (no re-QC).
"""
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2])); OUT=f"{ROOT}/results/tables/repro"
CORE8=["IL6","ADAMTS4","PLA2G2A","CYR61","PI15","SH3BP5","ERAP2","RGS16"]
PERI=["RGS5","PDGFRB","NOTCH3","KCNJ8","HIGD1B"]; SMC=["ACTA2","MYH11","TAGLN","CNN1","DES"]

a=ad.read_h5ad(f"{ROOT}/data/processed/GSE259348/GSE259348_annotated_ours.h5ad")
mural=a[a.obs.leiden.isin(["0","1","14"])].copy()          # SMC/pericyte mural clusters
print(f"mural cells: {mural.n_obs}")
# re-cluster mural only
sc.pp.highly_variable_genes(mural, n_top_genes=2000)
mv=mural[:,mural.var.highly_variable].copy(); sc.pp.scale(mv,max_value=10)
sc.tl.pca(mv,n_comps=30,random_state=0); sc.pp.neighbors(mv,n_neighbors=15,random_state=0)
sc.tl.leiden(mv,resolution=0.5,random_state=0,flavor="igraph",n_iterations=2,directed=False)
mural.obs["msub"]=mv.obs["leiden"].values
sub=mural.obs.groupby("msub",observed=True)[["sc_Pericyte","sc_SMC"]].mean()
sub["peri_minus_smc"]=sub.sc_Pericyte-sub.sc_SMC
sub["n"]=mural.obs.msub.value_counts()
print("\nmural subclusters (pericyte vs SMC score):"); print(sub.round(2).to_string())
peri_subs=list(sub.index[sub.peri_minus_smc>sub.peri_minus_smc.median()])  # pericyte-leaning half
# ensure we take the genuinely pericyte-leaning subclusters (positive contrast preferred)
if (sub.peri_minus_smc>0).any(): peri_subs=list(sub.index[sub.peri_minus_smc>0])
print("pericyte-like subclusters:", peri_subs)
peri=mural[mural.obs.msub.isin(peri_subs)]
print(f"GSE259348 pericyte-like cells: {peri.n_obs}")
print("  marker check (mean log-norm): PERI",
      {g:round(float(np.asarray(peri[:,g].X.todense()).mean()),2) for g in PERI if g in peri.var_names},
      "SMC", {g:round(float(np.asarray(peri[:,g].X.todense()).mean()),2) for g in SMC if g in peri.var_names})

pb=pd.Series(np.asarray(peri.layers["counts"].sum(0)).ravel(), index=a.var_names, name="GSE259348_DMED")

# combine with GSE206528 pericyte pseudobulk -> DESeq2 ~dataset+condition
cnt=pd.read_csv(f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_counts.tsv.gz",sep="\t",index_col=0)
meta=pd.read_csv(f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_metadata.tsv",sep="\t")
pm=meta[meta.mural_subtype=="Pericyte_like"].set_index("pseudobulk_id")
cols=[c for c in cnt.columns if c in pm.index]
gmap={"diabetic_ED":"diab_ED","non_diabetic_ED":"nondiab_ED","normal":"normal"}
g=cnt.index.intersection(pb.index)
comb=cnt.loc[g,cols].copy(); comb["GSE259348_DMED"]=pb.reindex(g).fillna(0).astype(int)
md=pd.DataFrame(index=comb.columns)
md["dataset"]=["GSE206528"]*len(cols)+["GSE259348"]
md["condition"]=pd.Categorical([gmap[pm.loc[c,"ed_subgroup"]] for c in cols]+["diab_ED"],
                               categories=["normal","nondiab_ED","diab_ED"])
X=comb.T.astype(int); keep=(X.sum(0)>=10)|X.columns.isin(CORE8); X=X.loc[:,keep]
dds=DeseqDataSet(counts=X,metadata=md,design="~dataset+condition",refit_cooks=False,quiet=True); dds.deseq2()
def ctr(x,y):
    st=DeseqStats(dds,contrast=["condition",x,y],quiet=True); st.summary(); return st.results_df.loc[CORE8]
dn=ctr("diab_ED","normal"); dne=ctr("diab_ED","nondiab_ED")
rep=pd.DataFrame({"lfc_DvsNormal":dn.log2FoldChange,"padj_DvsNormal":dn.padj,
                  "lfc_DvsNonDiabED":dne.log2FoldChange,"p_DvsNonDiabED":dne.pvalue}).loc[CORE8]
pd.set_option("display.float_format",lambda x:f"{x:.3g}","display.width",160)
print("\n===== Re-fetched GSE259348 + GSE206528 donor-level ~dataset+condition =====")
print(rep.to_string())
rep.to_csv(f"{OUT}/04_replication_combined_DE.csv")
print(f"\nsaved -> {OUT}/04_replication_combined_DE.csv")
print("NOTE: with a single GSE259348 sample the dataset term absorbs it -> condition effect = GSE206528-only")
print("      (same single-donor replication limitation). Real n=1 replication = directional concordance below.\n")

# ---- honest n=1 replication: directional concordance of per-donor logCPM across all 9 pericyte donors ----
pb.to_frame().to_csv(f"{OUT}/04_GSE259348_pericyte_pseudobulk.csv")
allc = comb.copy()                                   # genes x 9 donors (raw)
lib = allc.sum(0); logcpm = np.log2(allc.div(lib,axis=1)*1e6 + 1)
key=["SH3BP5","ERAP2","PLA2G2A","PI15","IL6","ADAMTS4"]
grp = {c:("diab_ED(206528)" if c.startswith(("DMED")) else
          "nondiab_ED" if c.startswith("ED") else
          "normal" if c.startswith("N") else "DIAB_ED(259348-repro)") for c in allc.columns}
tab = logcpm.loc[key].T; tab.insert(0,"group",[grp[c] for c in tab.index])
pd.set_option("display.float_format",lambda x:f"{x:.2f}","display.width",160)
print("Per-donor pericyte logCPM (9 donors incl. re-fetched GSE259348):")
print(tab.sort_values("group").to_string())
# rank of GSE259348 among all 9 for SH3BP5 (expect high) and ERAP2 (expect low)
for gname,hi in [("SH3BP5",True),("ERAP2",False)]:
    order=logcpm.loc[gname].sort_values(ascending=not hi)
    rank=list(order.index).index("GSE259348_DMED")+1
    print(f"  GSE259348 rank for {gname} ({'highest' if hi else 'lowest'} first): {rank}/9  "
          f"(value {logcpm.loc[gname,'GSE259348_DMED']:.2f})")
print("\nVerdict: n=1 directional concordance — does the independent diabetic donor sit SH3BP5-high / ERAP2-low?")
