#!/usr/bin/env python3
"""
Exact 2-vs-3 label permutation for DIABETES-SPECIFICITY of the pericyte signature,
within GSE206528 — the n=2-robust significance test (mirrors the archived
exact 2-vs-3 directional-score permutation, but for diabetic vs NON-diabetic ED).

Among the 5 organic-ED pericyte donors (2 diabetic: DMED1,DMED2; 3 non-diabetic: ED1,ED2,ED3),
ask whether the TRUE diabetic pair has the most extreme directional score among all C(5,2)=10
possible 2-vs-3 labelings. Exact one-sided p = rank/10.

Scores (per donor, z-scored log2CPM across the 5 ED donors):
  diabetes-associated  = z(SH3BP5) - z(ERAP2)            (both oriented: diabetic=high)
  shared-ED control  = mean z(IL6,ADAMTS4,PLA2G2A,CYR61,PI15)   (should NOT separate -> high p)
  full 12-gene       = directional score reference
Also single genes SH3BP5 (up) and ERAP2 (down).
"""
import itertools, os
import os
from pathlib import Path
import numpy as np, pandas as pd

ROOT = (os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
CNT = f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
META = f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_metadata.tsv"
OUT = f"{ROOT}/results/tables/repro"
os.makedirs(OUT, exist_ok=True)

UP = ["IL6","ADAMTS4","PLA2G2A","CYR61","PI15","NR4A1","NR4A2","NR4A3","MT1A","RGS16","SH3BP5"]  # up in diabetic ED
DOWN = ["ERAP2"]                                                                                 # down
SHARED = ["IL6","ADAMTS4","PLA2G2A","CYR61","PI15"]

meta = pd.read_csv(META, sep="\t")
peri = meta[meta.mural_subtype=="Pericyte_like"].set_index("pseudobulk_id")
counts = pd.read_csv(CNT, sep="\t", index_col=0)
cols = [c for c in counts.columns if c in peri.index]
cpm = counts[cols].div(counts[cols].sum(0), axis=1) * 1e6
logcpm = np.log2(cpm + 1)

donor = peri.loc[cols, "donor_id"]
sub   = peri.loc[cols, "ed_subgroup"]
ed_donors = [c for c in cols if sub[c] in ("diabetic_ED","non_diabetic_ED")]   # the 5 ED donors
diabetic = [c for c in ed_donors if sub[c]=="diabetic_ED"]                     # DMED1, DMED2

# z-score each gene across the 5 ED donors
sub_lc = logcpm[ed_donors]
z = sub_lc.sub(sub_lc.mean(1), axis=0).div(sub_lc.std(1).replace(0,np.nan), axis=0)

def score_diab_specific(c):  return z.loc["SH3BP5", c] - z.loc["ERAP2", c]
def score_shared(c):         return z.loc[SHARED, c].mean()
def score_full12(c):         return z.loc[UP, c].sum() - z.loc[DOWN, c].sum()
def score_gene(g, c, sign):  return sign*z.loc[g, c]

def exact_p(scorer):
    """one-sided exact p that the TRUE diabetic pair has the highest (group2_mean - group3_mean)."""
    s = {c: scorer(c) for c in ed_donors}
    stats = []
    for pair in itertools.combinations(ed_donors, 2):
        g2 = np.mean([s[c] for c in pair])
        g3 = np.mean([s[c] for c in ed_donors if c not in pair])
        stats.append((g2 - g3, set(pair)))
    stats.sort(key=lambda t: -t[0])                      # most-separated first
    true = set(diabetic)
    rank = [i for i,(_,p) in enumerate(stats) if p==true][0] + 1
    obs  = [v for v,p in stats if p==true][0]
    return rank, len(stats), rank/len(stats), obs

print(f"5 ED pericyte donors: {[donor[c] for c in ed_donors]}  (diabetic = {[donor[c] for c in diabetic]})\n")
rows=[]
for name, scorer in [
    ("diabetes-associated (SH3BP5 - ERAP2)", score_diab_specific),
    ("SH3BP5 (up)",            lambda c: score_gene('SH3BP5', c, +1)),
    ("ERAP2 (down)",           lambda c: score_gene('ERAP2',  c, -1)),
    ("shared-ED control (IL6/ADAMTS4/PLA2G2A/CYR61/PI15)", score_shared),
    ("full 12-gene directional", score_full12),
]:
    rank, n, p, obs = exact_p(scorer)
    rows.append({"score":name, "true_pair_rank":f"{rank}/{n}", "exact_one_sided_p":round(p,3), "stat":round(obs,3)})
    print(f"{name:52s}  rank {rank}/{n}   exact p={p:.3f}   stat={obs:+.2f}")

pd.DataFrame(rows).to_csv(f"{OUT}/02_exact_permutation_diabetes_specificity.csv", index=False)
print(f"\nsaved -> {OUT}/02_exact_permutation_diabetes_specificity.csv")
print("\nPer-donor diabetes-associated score (SH3BP5 - ERAP2):")
for c in ed_donors:
    print(f"  {donor[c]:6s} ({sub[c]:16s})  {score_diab_specific(c):+.2f}")
