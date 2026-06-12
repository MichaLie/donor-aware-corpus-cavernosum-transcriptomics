#!/usr/bin/env python3
"""
Reproduce the stress decomposition within GSE206528 using the archived pipeline.

Question: is the diabetes-associated axis (SH3BP5 up / ERAP2 down) independent of the
immediate-early/stress program (NR4A1/2/3, MT1A, RGS16), or just a stress readout?

Steps:
 1. Load preprocessed (log-norm) GSE206528, join mural annotations, keep Pericyte_like cells.
 2. Per-cell scores (scanpy score_genes): stress5, and per-cell log-norm expr of SH3BP5/ERAP2/NR4A3.
 3. Stress-correct each per cell: residual of (gene/score) ~ stress5 across all pericytes (OLS).
 4. Aggregate residuals to per-donor means; run the exact 2-vs-3 permutation (diabetic vs non-diab ED)
    on the STRESS-CORRECTED scores; compare to uncorrected.
 Expect: SH3BP5/ERAP2 axis stays separated (p=0.10) after correction; NR4A3 (stress) loses separation.
"""
import itertools, os
import os
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc

ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
H5=f"{ROOT}/data/processed/GSE206528/GSE206528_scanpy_preprocessed.h5ad"
ANN=f"{ROOT}/results/tables/GSE206528_cell_mural_subtype_annotations.tsv.gz"
OUT=f"{ROOT}/results/tables/repro"; os.makedirs(OUT, exist_ok=True)
STRESS5=["NR4A1","NR4A2","NR4A3","MT1A","RGS16"]

x=ad.read_h5ad(H5)
ann=pd.read_csv(ANN, sep="\t").rename(columns={"Unnamed: 0":"cell"}).set_index("cell")
x=x[x.obs_names.isin(ann.index)].copy()
x.obs["mural_subtype"]=ann.loc[x.obs_names,"mural_subtype"].values
x.obs["donor"]=ann.loc[x.obs_names,"donor_id"].values
x.obs["grp"]=ann.loc[x.obs_names,"ed_subgroup"].values
peri=x[x.obs.mural_subtype=="Pericyte_like"].copy()
print(f"Pericyte_like cells: {peri.n_obs}")
print("per-donor pericyte cells:", peri.obs.groupby('donor', observed=True).size().to_dict())

sc.tl.score_genes(peri, [g for g in STRESS5 if g in peri.var_names], score_name="stress5", ctrl_size=50, random_state=0)
def expr(g):
    import numpy as np
    return np.asarray(peri[:, g].X.todense()).ravel() if hasattr(peri[:,g].X,"todense") else np.asarray(peri[:,g].X).ravel()
df=pd.DataFrame({"donor":peri.obs.donor.values,"grp":peri.obs["grp"].values,
                 "stress5":peri.obs.stress5.values,
                 "SH3BP5":expr("SH3BP5"),"ERAP2":expr("ERAP2"),"NR4A3":expr("NR4A3")})
df["diab_specific"]=df.SH3BP5 - df.ERAP2

def residualize(y, s):                    # OLS residual of y ~ stress5 (remove stress component)
    A=np.c_[np.ones_like(s), s]; beta,_,_,_=np.linalg.lstsq(A,y,rcond=None); return y-A@beta
for col in ["SH3BP5","ERAP2","NR4A3","diab_specific"]:
    df[col+"_corr"]=residualize(df[col].values, df["stress5"].values)

donmean=df.groupby("donor").mean(numeric_only=True)
donmean["grp"]=df.groupby("donor")["grp"].first()
ed=donmean[donmean["grp"].isin(["diabetic_ED","non_diabetic_ED"])].copy()
diabetic=list(ed.index[ed["grp"]=="diabetic_ED"])

def exact_p(series):                      # one-sided: true diabetic pair most extreme (g2-g3)
    s=series.to_dict(); donors=list(series.index); stats=[]
    for pair in itertools.combinations(donors,2):
        g2=np.mean([s[c] for c in pair]); g3=np.mean([s[c] for c in donors if c not in pair])
        stats.append((g2-g3,set(pair)))
    stats.sort(key=lambda t:-t[0]); true=set(diabetic)
    rank=[i for i,(_,p) in enumerate(stats) if p==true][0]+1
    return rank,len(stats),rank/len(stats)

print("\nPer-donor mean stress5 by group (justifies the concern):")
print(donmean.groupby("grp")["stress5"].mean().to_string())
print("\nExact 2-vs-3 (diabetic vs non-diabetic ED), z-scored across 5 ED donors, UNcorrected vs STRESS-corrected:")
rows=[]
for label,raw,corr in [("SH3BP5","SH3BP5","SH3BP5_corr"),("ERAP2(down)","ERAP2","ERAP2_corr"),
                       ("NR4A3(stress)","NR4A3","NR4A3_corr"),("diab-specific SH3BP5-ERAP2","diab_specific","diab_specific_corr")]:
    zr=(ed[raw]-ed[raw].mean())/ed[raw].std()
    zc=(ed[corr]-ed[corr].mean())/ed[corr].std()
    sign = -1 if "down" in label else 1     # ERAP2 oriented down=diabetic
    r0=exact_p(sign*zr); r1=exact_p(sign*zc)
    rows.append({"score":label,"p_uncorrected":round(r0[2],3),"rank_unc":f"{r0[0]}/{r0[1]}",
                 "p_stress_corrected":round(r1[2],3),"rank_corr":f"{r1[0]}/{r1[1]}"})
    print(f"  {label:28s} uncorrected p={r0[2]:.3f} ({r0[0]}/{r0[1]})   stress-corrected p={r1[2]:.3f} ({r1[0]}/{r1[1]})")
pd.DataFrame(rows).to_csv(f"{OUT}/03_stress_corrected_exact_GSE206528.csv", index=False)
print(f"\nsaved -> {OUT}/03_stress_corrected_exact_GSE206528.csv")

# ---- per-cell effect-size attenuation: diabetic vs non-diabetic ED, raw vs stress-corrected ----
# (descriptive effect size on thousands of pericytes; donor-level exact test above is the inferential part)
cells = df[df["grp"].isin(["diabetic_ED","non_diabetic_ED"])].copy()
isd = (cells["grp"]=="diabetic_ED").values
def cohend(v):
    a,b=v[isd],v[~isd]
    sp=np.sqrt(((a.var(ddof=1)*(len(a)-1))+(b.var(ddof=1)*(len(b)-1)))/(len(a)+len(b)-2))
    return (a.mean()-b.mean())/sp if sp>0 else np.nan
print("\nPer-cell effect size (Cohen's d, diabetic vs non-diabetic ED pericytes), RAW vs STRESS-CORRECTED:")
att=[]
for g in ["SH3BP5","ERAP2","NR4A3","diab_specific"]:
    d_raw=cohend(cells[g].values); d_corr=cohend(cells[g+"_corr"].values)
    ret = d_corr/d_raw if abs(d_raw)>1e-9 else np.nan
    att.append({"gene":g,"d_raw":round(d_raw,3),"d_stress_corrected":round(d_corr,3),"pct_retained":round(100*ret,0)})
    print(f"  {g:14s} d_raw={d_raw:+.2f}  d_corr={d_corr:+.2f}  retained={100*ret:4.0f}%")
pd.DataFrame(att).to_csv(f"{OUT}/03b_stress_attenuation_effectsize.csv", index=False)
print(f"saved -> {OUT}/03b_stress_attenuation_effectsize.csv")
