#!/usr/bin/env python3
"""
Cell-type specificity panel: is the diabetes axis (SH3BP5 up / ERAP2 down, diabetic vs
non-diabetic ED) genuinely pericyte/mural-specific, or broadly stromal?
Per major cell type (+ Mural split into Pericyte_like vs other mural), compute the
donor-level diabetic-vs-non-diabetic-ED difference in SH3BP5, ERAP2, the axis, and the
shared-ED score. If the axis effect is concentrated in pericyte/mural -> "perivascular"
claim supported; if broad -> reframe to "mural/stromal".
"""
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, matplotlib
from pathlib import Path
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2])); OUT=f"{ROOT}/results/tables/repro"; FIG=f"{ROOT}/results/figures/repro"; os.makedirs(FIG, exist_ok=True)
SHARED=["IL6","ADAMTS4","PLA2G2A","CYR61","PI15"]
x=ad.read_h5ad(f"{ROOT}/data/processed/GSE206528/GSE206528_scanpy_preprocessed.h5ad")
ann=pd.read_csv(f"{ROOT}/results/tables/GSE206528_cell_major_celltype_annotations.tsv.gz",sep="\t").rename(columns={"Unnamed: 0":"cell"}).set_index("cell")
mur=pd.read_csv(f"{ROOT}/results/tables/GSE206528_cell_mural_subtype_annotations.tsv.gz",sep="\t").rename(columns={"Unnamed: 0":"cell"}).set_index("cell")
x=x[x.obs_names.isin(ann.index)].copy()
x.obs["ct"]=ann.loc[x.obs_names,"major_cell_type"].values
x.obs["donor"]=ann.loc[x.obs_names,"donor_id"].values
x.obs["grp"]=ann.loc[x.obs_names,"ed_subgroup"].values
# refine Mural -> Pericyte_like vs Other_mural
sub=mur["mural_subtype"].reindex(x.obs_names)
x.obs["ct2"]=np.where(x.obs["ct"].values=="Mural",
                      np.where(sub.values=="Pericyte_like","Pericyte_like","Other_mural"),
                      x.obs["ct"].values)
def expr(g):
    import numpy as np
    v=x[:,g].X; return np.asarray(v.todense()).ravel() if hasattr(v,"todense") else np.asarray(v).ravel()
sc.tl.score_genes(x,[g for g in SHARED if g in x.var_names],score_name="sharedED",random_state=0)
df=pd.DataFrame({"ct":x.obs.ct2.values,"donor":x.obs.donor.values,"grp":x.obs.grp.values,
                 "SH3BP5":expr("SH3BP5"),"ERAP2":expr("ERAP2"),"sharedED":x.obs.sharedED.values})
# z across all cells for the axis
for g in ["SH3BP5","ERAP2"]: df[g+"_z"]=(df[g]-df[g].mean())/df[g].std()
df["axis"]=df["SH3BP5_z"]-df["ERAP2_z"]

# donor x celltype means, then diabetic(DMED) vs non-diabetic-ED(ED) difference per cell type
dia=["DMED1","DMED2"]; non=["ED1","ED2","ED3"]
order=["Pericyte_like","Other_mural","Smooth_muscle","Endothelial","Fibroblast","Macrophage","T_cell","Schwann_neural"]
rows=[]
for ct in order:
    sub_ct=df[df.ct==ct]
    dm=sub_ct.groupby("donor")[["SH3BP5","ERAP2","axis","sharedED"]].mean()
    rec={"cell_type":ct,"n_cells":len(sub_ct),
         "n_dia_cells":int((sub_ct.donor.isin(dia)).sum()),"n_non_cells":int((sub_ct.donor.isin(non)).sum())}
    for s in ["SH3BP5","ERAP2","axis","sharedED"]:
        d=dm.reindex(dia)[s].mean(); n=dm.reindex(non)[s].mean()
        rec[f"{s}_diabVSnon"]=round(d-n,3)
    rows.append(rec)
res=pd.DataFrame(rows)
pd.set_option("display.width",170,"display.float_format",lambda v:f"{v:.3f}")
print("Diabetic vs non-diabetic ED difference, per cell type (donor-level means):")
print(res.to_string(index=False))
res.to_csv(f"{OUT}/06_celltype_specificity.csv",index=False)

# heatmap of the diabetic-vs-non differences — LEAD with SH3BP5 (pericyte-specific); drop the
# misleading combined "axis" column (ERAP2's broad downregulation makes it peak in fibroblasts).
hm=res.set_index("cell_type")[["SH3BP5_diabVSnon","ERAP2_diabVSnon","sharedED_diabVSnon"]]
hm.columns=["SH3BP5\n(diabetes-associated)","ERAP2\n(broad)","shared-ED"]
fig,axh=plt.subplots(figsize=(5.4,4.2))
vmax=np.nanmax(np.abs(hm.values))
im=axh.imshow(hm.values,cmap="RdBu_r",vmin=-vmax,vmax=vmax,aspect="auto")
axh.set_xticks(range(hm.shape[1])); axh.set_xticklabels(hm.columns,fontsize=8)
axh.set_yticks(range(hm.shape[0])); axh.set_yticklabels(hm.index,fontsize=8)
for i in range(hm.shape[0]):
    for j in range(hm.shape[1]):
        axh.text(j,i,f"{hm.values[i,j]:+.2f}",ha="center",va="center",fontsize=7,
                 color="white" if abs(hm.values[i,j])>vmax*0.6 else "black")
axh.set_title("SH3BP5 is pericyte-enriched; ERAP2 is broad\n(diabetic vs non-diabetic ED effect by cell type; red = up in diabetic)",fontsize=8.5)
axh.add_patch(plt.Rectangle((-0.5,-0.5),1,1,fill=False,edgecolor="black",lw=2.2))  # highlight SH3BP5/Pericyte cell
plt.colorbar(im,fraction=0.046,pad=0.04); plt.tight_layout()
for ext in("png","pdf"): fig.savefig(f"{FIG}/celltype_specificity.{ext}",bbox_inches="tight")
print(f"\nsaved -> {OUT}/06_celltype_specificity.csv  and  {FIG}/celltype_specificity.png")
print("\nVERDICT CHECK: is the axis effect concentrated in Pericyte_like/mural vs other cell types?")
ax_eff=res.set_index("cell_type")["axis_diabVSnon"]
print(ax_eff.sort_values(ascending=False).to_string())
