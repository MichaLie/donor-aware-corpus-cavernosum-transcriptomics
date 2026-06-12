#!/usr/bin/env python3
"""
Donor QC table + sensitivity analyses to show the SH3BP5 diabetic-vs-non-diabetic-ED
pericyte effect is not an obvious artifact (cell number, depth, mito, stress, annotation).
"""
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad, itertools
from pathlib import Path
ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2])); OUT=f"{ROOT}/results/tables/repro"
PB=f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes"
x=ad.read_h5ad(f"{ROOT}/data/processed/GSE206528/GSE206528_scanpy_preprocessed.h5ad")
ann=pd.read_csv(f"{ROOT}/results/tables/GSE206528_cell_mural_subtype_annotations.tsv.gz",sep="\t").rename(columns={"Unnamed: 0":"cell"}).set_index("cell")
x.obs["mural"]=ann["mural_subtype"].reindex(x.obs_names)
x.obs["donor"]=ann["donor_id"].reindex(x.obs_names)
x.obs["grp"]=ann["ed_subgroup"].reindex(x.obs_names)
peri=x[x.obs.mural=="Pericyte_like"].copy()
sc.tl.score_genes(peri,["NR4A1","NR4A2","NR4A3","MT1A","RGS16"],score_name="stress5",random_state=0)
sc.tl.score_genes(peri,["RGS5","PDGFRB","NOTCH3"],score_name="peri_marker",random_state=0)
def e(obj,g):
    v=obj[:,g].X; return np.asarray(v.todense()).ravel() if hasattr(v,"todense") else np.asarray(v).ravel()
TYPE={"DMED1":"T1D long-standing","DMED2":"T1D long-standing","ED1":"non-diabetic","ED2":"non-diabetic","ED3":"non-diabetic","N1":"none","N2":"none","N3":"none"}
SRC={"DMED1":"prosthesis biopsy","DMED2":"prosthesis biopsy","ED1":"prosthesis biopsy","ED2":"prosthesis biopsy","ED3":"prosthesis biopsy","N1":"tumor-margin","N2":"tumor-margin","N3":"tumor-margin"}
od=peri.obs
qc=[]
for d,g in od.groupby("donor",observed=True):
    sub=peri[peri.obs.donor==d]
    qc.append({"donor":d,"group":od[od.donor==d]["grp"].iloc[0],"diabetes":TYPE.get(d,""),"tissue":SRC.get(d,""),
               "pericyte_cells":sub.n_obs,
               "median_genes":int(np.median(sub.obs.n_genes_by_counts)),
               "median_UMIs":int(np.median(sub.obs.total_counts)),
               "mito_pct":round(float((sub.obs.mt_counts/sub.obs.total_counts*100).mean()),2),
               "stress_score":round(float(sub.obs.stress5.mean()),3),
               "peri_marker_score":round(float(sub.obs.peri_marker.mean()),3),
               "SH3BP5_mean":round(float(e(sub,"SH3BP5").mean()),3)})
qc=pd.DataFrame(qc).set_index("donor").reindex(["N1","N2","N3","ED1","ED2","ED3","DMED1","DMED2"])
pd.set_option("display.width",200)
print("=== Donor-level QC table (Pericyte_like) ===")
print(qc.to_string())
qc.to_csv(f"{OUT}/07_donor_QC_table.csv")

# ---- correlate SH3BP5 with QC nuisance vars across donors (is the signal an artifact?) ----
print("\n=== Is SH3BP5 confounded by a QC variable? (Spearman across 8 donors) ===")
for v in ["pericyte_cells","median_genes","median_UMIs","mito_pct"]:
    print(f"  SH3BP5_mean vs {v:16s}: rho={qc['SH3BP5_mean'].corr(qc[v],method='spearman'):+.2f}")
print("  (diabetic donors DMED1/2 should top SH3BP5 regardless of these)")

# ---- equal-cell downsampling: recompute SH3BP5 diabetic-vs-nonED, multiple seeds ----
print("\n=== Equal-cell downsampling (per-donor pericytes -> N each), SH3BP5 diabetic vs non-diabetic ED ===")
ed_donors=["DMED1","DMED2","ED1","ED2","ED3"]; dia=["DMED1","DMED2"]
counts=od[od.donor.isin(ed_donors)].groupby("donor",observed=True).size()
N=int(counts.drop("ED2").min()) if "ED2" in counts.index else int(counts.min())   # ED2=16 is an outlier-low profile
print(f"  downsample N={N}/donor (ED2 has only {counts.get('ED2','?')} cells; included only if >=N)")
use=[d for d in ed_donors if counts.get(d,0)>=N]
for seed in range(5):
    rng=np.random.default_rng(seed); vals={}
    for d in use:
        idx=np.where(peri.obs.donor.values==d)[0]; pick=rng.choice(idx,N,replace=False)
        vals[d]=float(e(peri[pick],"SH3BP5").mean())
    diff=np.mean([vals[d] for d in dia])-np.mean([vals[d] for d in use if d not in dia])
    print(f"  seed{seed}: SH3BP5 diabetic-mean − nonED-mean = {diff:+.3f}  (per-donor: "+", ".join(f'{d}={vals[d]:.2f}' for d in use)+")")

# ---- stricter downsampling: ALL five ED donors at the minimum cell count (incl. ED2=16), many seeds ----
Nmin=int(counts.min())   # = ED2's 16 pericytes; include ALL five ED donors
print(f"\n=== Stricter equal-cell downsampling: N={Nmin}/donor including ALL {len(ed_donors)} ED donors (100 seeds) ===")
rows16=[]
for seed in range(100):
    rng=np.random.default_rng(seed); vals={}
    for d in ed_donors:
        idx=np.where(peri.obs.donor.values==d)[0]; pick=rng.choice(idx,Nmin,replace=False)
        vals[d]=float(e(peri[pick],"SH3BP5").mean())
    diff=np.mean([vals[d] for d in dia])-np.mean([vals[d] for d in ed_donors if d not in dia])
    rows16.append({"seed":seed,"SH3BP5_diab_minus_nonED":round(diff,4),**{f"SH3BP5_{d}":round(vals[d],4) for d in ed_donors}})
ds=pd.DataFrame(rows16); ds.to_csv(f"{OUT}/07_downsampling_16cells_alldonors.csv",index=False)
d16=ds["SH3BP5_diab_minus_nonED"]
print(f"  SH3BP5 diabetic−nonED across 100 seeds: min {d16.min():+.3f} / mean {d16.mean():+.3f} / max {d16.max():+.3f} "
      f"({(d16>0).mean()*100:.0f}% of seeds positive)")
print(f"  saved -> {OUT}/07_downsampling_16cells_alldonors.csv")

# ---- alternative annotation: broad mural vs pericyte_like ----
print("\n=== Alternative pericyte/mural definitions: SH3BP5 diabetic vs non-diabetic ED (donor-mean diff) ===")
for label,mask in [("Pericyte_like (primary)", x.obs.mural=="Pericyte_like"),
                   ("Broad mural (all subtypes)", x.obs.mural.isin(["Pericyte_like","Vessel_SMC_like","CC_SMC_like"]))]:
    sub=x[mask & x.obs.donor.isin(ed_donors)]
    dm=pd.Series(e(sub,"SH3BP5"),index=sub.obs.donor.values).groupby(level=0).mean()
    diff=dm.reindex(dia).mean()-dm.reindex([d for d in ed_donors if d not in dia]).mean()
    print(f"  {label:30s}: SH3BP5 diabetic−nonED = {diff:+.3f}")
print(f"\nsaved -> {OUT}/07_donor_QC_table.csv")
