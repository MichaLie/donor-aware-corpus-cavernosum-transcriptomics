#!/usr/bin/env python3
"""
GSE131882 nephropathy cross-tissue: HONEST donor-level reanalysis of the kidney-pericyte
"non-stress core" signal. A previously circulated headline p=3.8e-6 was cell-level
(pseudoreplication); this script keeps donor as the inferential unit.
Here: donor = unit. 3 diabetic vs 3 control donors, pseudobulk over the refined mural/pericyte cluster.

Inputs (built by 05b/05c in Docker from re-fetched GEO RDS):
  data/processed/GSE131882/GSE131882_mural_pseudobulk_counts_bySymbol_refined.csv   (gene x 6 donors)
  data/processed/GSE131882/GSE131882_mural_percell_activation_lognorm_refined.csv   (per-cell, for the
        DESCRIPTIVE cell-level effect/pseudorep contrast ONLY -- not inferential)

Non-stress core = pericyte activation signature MINUS the stress/immediate-early genes.
  activation (shared-ED) genes : IL6, ADAMTS4, PLA2G2A, CYR61/CCN1, PI15  (+ SH3BP5, ERAP2 diabetes axis)
  stress/IEG genes (excluded)  : NR4A1, NR4A2, NR4A3, MT1A, RGS16
Of these, present & non-trivially measurable in GSE131882 kidney mural cells: CCN1(=CYR61), ADAMTS4,
  SH3BP5, ERAP2 (PI15/PLA2G2A/IL6 ~absent). So the testable non-stress core here = {CCN1, ADAMTS4,
  SH3BP5, ERAP2}; RGS16 is a stress gene reported separately, NOT inside the non-stress core.

Donor-level tests:
  (1) pseudobulk DESeq2 ~condition (diabetes vs control) on all genes -> LFC/p for the core genes.
  (2) per-donor non-stress-core SCORE = mean of per-gene z-scored logCPM (across 6 donors);
      exact 3-vs-3 label permutation: is the true diabetic trio the most extreme? all-above-all-control?
  (3) per-gene "all-diabetic-above-all-control?" separation check.
Contrast explicitly with the cell-level Mann-Whitney p (which we compute, label INVALID = pseudoreplication).
"""
import itertools, os
import os
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
PROC=f"{ROOT}/data/processed/GSE131882"
OUT=f"{ROOT}/results/tables/repro"; os.makedirs(OUT, exist_ok=True)

DONORS=["control_s1","control_s2","control_s3","diabetes_s1","diabetes_s2","diabetes_s3"]
COND={d:("diabetes" if d.startswith("diabetes") else "control") for d in DONORS}
DIAB=[d for d in DONORS if COND[d]=="diabetes"]
CTRL=[d for d in DONORS if COND[d]=="control"]

# activation / non-stress-core vs stress sets (full reframe vocabulary)
ACTIVATION = ["IL6","ADAMTS4","PLA2G2A","CCN1","PI15","SH3BP5","ERAP2"]   # CCN1 == CYR61
STRESS     = ["NR4A1","NR4A2","NR4A3","MT1A","RGS16"]

pb = pd.read_csv(f"{PROC}/GSE131882_mural_pseudobulk_counts_bySymbol_refined.csv", index_col=0)
pb = pb[DONORS]
print("pseudobulk shape (genes x donors):", pb.shape)
print("library size (total counts) per donor:\n", pb.sum(0).to_string())

# which activation genes are actually measurable (nonzero in >=2 donors)
measurable = [g for g in ACTIVATION if g in pb.index and (pb.loc[g]>0).sum()>=2]
print("\nActivation genes present & measurable (nonzero in >=2 donors):", measurable)
print("Activation genes absent/too-sparse:", [g for g in ACTIVATION if g not in measurable])
core = measurable  # the testable NON-STRESS core in this dataset

# ---------------- (1) pseudobulk DESeq2 diabetes vs control ----------------
meta = pd.DataFrame({"condition": pd.Categorical([COND[d] for d in DONORS],
                                                 categories=["control","diabetes"])}, index=DONORS)
X = pb.T.astype(int)
keep = (X.sum(0) >= 10) | X.columns.isin(core+STRESS)
X = X.loc[:, keep]
dds = DeseqDataSet(counts=X, metadata=meta, design="~condition", refit_cooks=False, quiet=True)
dds.deseq2()
st = DeseqStats(dds, contrast=["condition","diabetes","control"], quiet=True); st.summary()
res = st.results_df
report_genes = [g for g in (core+STRESS+["RGS5","PDGFRB","NOTCH3"]) if g in res.index]
de = res.loc[report_genes, ["log2FoldChange","pvalue","padj","baseMean"]].copy()
print("\n===== (1) Pseudobulk DESeq2  diabetes vs control (donor-level, n=3 vs 3) =====")
pd.set_option("display.float_format", lambda v:f"{v:.3g}", "display.width", 160)
print(de.to_string())

# ---------------- per-donor logCPM ----------------
lib = pb.sum(0); logcpm = np.log2(pb.div(lib, axis=1)*1e6 + 1)

# ---------------- (2) per-donor non-stress-core SCORE + exact 3-vs-3 permutation ----------------
def zrows(genes):
    z = logcpm.loc[genes].copy()
    z = z.sub(z.mean(1), axis=0).div(z.std(1).replace(0,np.nan), axis=0).dropna(how="all")
    return z
zc = zrows(core)
core_score = zc.mean(0)                       # per-donor mean z of non-stress core
print("\n===== (2) Per-donor NON-STRESS CORE score (mean z of", core, ") =====")
tab = pd.DataFrame({"condition":[COND[d] for d in DONORS], "core_score":core_score.loc[DONORS]})
print(tab.to_string())

def exact_perm(score, diab):
    """one-sided: how extreme is the TRUE diabetic trio (mean_diab - mean_rest) among all C(6,3)=20 splits?"""
    d=score.to_dict(); donors=list(score.index); stats_=[]
    for trio in itertools.combinations(donors,3):
        m_in=np.mean([d[c] for c in trio]); m_out=np.mean([d[c] for c in donors if c not in trio])
        stats_.append((m_in-m_out, set(trio)))
    stats_.sort(key=lambda t:-t[0]); true=set(diab)
    rank=[i for i,(_,s) in enumerate(stats_) if s==true][0]+1
    # two-sided exact p for a 3-vs-3 design = rank-based; report one-sided (rank/20) which is the natural test
    return rank, len(stats_), rank/len(stats_)
r,n,p = exact_perm(core_score, DIAB)
print(f"\nExact 3-vs-3 permutation on the non-stress-core score: true diabetic trio rank {r}/{n}  (one-sided p={p:.3f})")
allabove = core_score.loc[DIAB].min() > core_score.loc[CTRL].max()
print(f"all-diabetic-above-all-control on core score? {allabove}  "
      f"(min diab {core_score.loc[DIAB].min():.2f} vs max ctrl {core_score.loc[CTRL].max():.2f})")

# per-gene contribution: which genes drive it (mean diab logCPM - mean ctrl logCPM), and separation
print("\n===== per-gene drivers (donor-level) =====")
drv=[]
for g in core + ["RGS16"]:   # include RGS16 as a stress-gene comparator
    if g not in logcpm.index: continue
    md=logcpm.loc[g,DIAB].mean(); mc=logcpm.loc[g,CTRL].mean()
    sep = logcpm.loc[g,DIAB].min() > logcpm.loc[g,CTRL].max()
    # leave-one-out: does dropping the single highest diabetic donor kill the difference?
    drv.append({"gene":g, "in_core": g in core,
                "mean_logCPM_diab":round(md,3),"mean_logCPM_ctrl":round(mc,3),
                "delta_diab_minus_ctrl":round(md-mc,3),
                "all_diab_above_all_ctrl":bool(sep),
                "diab_vals":";".join(f"{logcpm.loc[g,d]:.2f}" for d in DIAB),
                "ctrl_vals":";".join(f"{logcpm.loc[g,d]:.2f}" for d in CTRL)})
drv=pd.DataFrame(drv); print(drv.to_string(index=False))

# leave-one-diabetic-donor-out sensitivity on the core score difference
print("\n===== leave-one-out (drop one diabetic donor): does the core-score gap survive? =====")
loo=[]
base_gap = core_score.loc[DIAB].mean()-core_score.loc[CTRL].mean()
for drop in DIAB:
    keepd=[d for d in DIAB if d!=drop]
    gap = core_score.loc[keepd].mean()-core_score.loc[CTRL].mean()
    loo.append({"dropped":drop,"core_gap":round(gap,3)})
print(f"full core gap (diab-ctrl): {base_gap:.3f}")
print(pd.DataFrame(loo).to_string(index=False))

# ---------------- (3) CELL-LEVEL contrast (invalid inference; descriptive only) ----------------
print("\n===== (3) CELL-LEVEL contrast — INVALID inference (pseudoreplication), shown only to expose it =====")
pc = pd.read_csv(f"{PROC}/GSE131882_mural_percell_activation_lognorm_refined.csv")
isd = (pc["cond"]=="diabetes").values
cell_rows=[]
for g in core + ["RGS16"]:
    if g not in pc.columns: continue
    a=pc.loc[isd,g].values; b=pc.loc[~isd,g].values
    U,pmw = stats.mannwhitneyu(a,b,alternative="two-sided")
    # Cohen's d (descriptive effect size, allowed)
    sp=np.sqrt(((a.var(ddof=1)*(len(a)-1))+(b.var(ddof=1)*(len(b)-1)))/(len(a)+len(b)-2))
    d=(a.mean()-b.mean())/sp if sp>0 else np.nan
    cell_rows.append({"gene":g,"n_diab_cells":len(a),"n_ctrl_cells":len(b),
                      "cell_MannWhitney_p_INVALID":f"{pmw:.2e}","cohens_d_descr":round(d,3)})
cell=pd.DataFrame(cell_rows); print(cell.to_string(index=False))
# combined non-stress-core per-cell score, cell-level MWU, to mirror the non-donor-level summary
zcols=[g for g in core if g in pc.columns]
zc_cell = (pc[zcols]-pc[zcols].mean())/pc[zcols].std(ddof=0)
pc["core_cellscore"]=zc_cell.mean(1)
Uc,pc_mwu = stats.mannwhitneyu(pc.loc[isd,"core_cellscore"], pc.loc[~isd,"core_cellscore"], alternative="two-sided")
print(f"\nNON-STRESS-CORE combined per-cell score, CELL-LEVEL Mann-Whitney p = {pc_mwu:.2e}  "
      f"(n={isd.sum()} diabetic vs {(~isd).sum()} control CELLS)  <-- this is the pseudoreplicated kind of number")
print("   -> Honest DONOR-level result for the SAME score: exact 3-vs-3 p =",
      f"{p:.3f} (rank {r}/{n}), all-diab-above-all-ctrl={allabove}.")

# ---------------- assemble output CSV ----------------
summary_rows=[]
for g in core + ["RGS16"]:
    if g not in logcpm.index: continue
    row={"gene":g, "set":("non_stress_core" if g in core else "stress(excluded from core)"),
         "deseq2_log2FC_diab_vs_ctrl": round(float(res.loc[g,"log2FoldChange"]),3) if g in res.index else np.nan,
         "deseq2_pvalue": float(res.loc[g,"pvalue"]) if g in res.index else np.nan,
         "deseq2_padj":   float(res.loc[g,"padj"]) if g in res.index else np.nan,
         "mean_logCPM_diab": round(logcpm.loc[g,DIAB].mean(),3),
         "mean_logCPM_ctrl": round(logcpm.loc[g,CTRL].mean(),3),
         "delta_logCPM": round(logcpm.loc[g,DIAB].mean()-logcpm.loc[g,CTRL].mean(),3),
         "all_diab_above_all_ctrl": bool(logcpm.loc[g,DIAB].min()>logcpm.loc[g,CTRL].max())}
    cr = cell[cell.gene==g]
    row["cell_MannWhitney_p_INVALID"] = cr["cell_MannWhitney_p_INVALID"].values[0] if len(cr) else ""
    row["cohens_d_cell_descr"] = cr["cohens_d_descr"].values[0] if len(cr) else np.nan
    summary_rows.append(row)
out = pd.DataFrame(summary_rows)
# append the headline donor-level vs cell-level rows as metadata
meta_rows = pd.DataFrame([
 {"gene":"__NONSTRESS_CORE_SCORE__","set":"summary",
  "deseq2_log2FC_diab_vs_ctrl":np.nan,"deseq2_pvalue":np.nan,"deseq2_padj":np.nan,
  "mean_logCPM_diab":round(core_score.loc[DIAB].mean(),3),
  "mean_logCPM_ctrl":round(core_score.loc[CTRL].mean(),3),
  "delta_logCPM":round(core_score.loc[DIAB].mean()-core_score.loc[CTRL].mean(),3),
  "all_diab_above_all_ctrl":bool(allabove),
  "cell_MannWhitney_p_INVALID":f"{pc_mwu:.2e}",
  "cohens_d_cell_descr":np.nan},
 {"gene":"__DONOR_exact_3v3_p__","set":"summary","deseq2_log2FC_diab_vs_ctrl":np.nan,
  "deseq2_pvalue":round(p,3),"deseq2_padj":np.nan,"mean_logCPM_diab":np.nan,"mean_logCPM_ctrl":np.nan,
  "delta_logCPM":np.nan,"all_diab_above_all_ctrl":bool(allabove),
  "cell_MannWhitney_p_INVALID":f"rank {r}/{n}","cohens_d_cell_descr":np.nan},
])
# also bring in the loose-gate (cluster-9, 621-cell) cell-level core score for transparency
try:
    pcl = pd.read_csv(f"{PROC}/GSE131882_pericyte_percell_activation_lognorm.csv")
    loose_core = [g for g in ACTIVATION if g in pcl.columns]
    zl = (pcl[loose_core]-pcl[loose_core].mean())/pcl[loose_core].std(ddof=0)
    pcl["score"]=zl.mean(1); isdl=(pcl["cond"]=="diabetes").values
    _,p_loose2 = stats.mannwhitneyu(pcl.loc[isdl,"score"],pcl.loc[~isdl,"score"],alternative="two-sided")
    _,p_loose1 = stats.mannwhitneyu(pcl.loc[isdl,"score"],pcl.loc[~isdl,"score"],alternative="greater")
    loose_note = (pd.DataFrame([
      {"gene":"__LOOSE_gate_cellscore_MWU__","set":"summary","deseq2_log2FC_diab_vs_ctrl":np.nan,
       "deseq2_pvalue":np.nan,"deseq2_padj":np.nan,"mean_logCPM_diab":np.nan,"mean_logCPM_ctrl":np.nan,
       "delta_logCPM":np.nan,"all_diab_above_all_ctrl":np.nan,
       "cell_MannWhitney_p_INVALID":f"2sided={p_loose2:.2e};1sided={p_loose1:.2e};n={isdl.sum()}v{(~isdl).sum()}",
       "cohens_d_cell_descr":np.nan}]))
    meta_rows = pd.concat([meta_rows, loose_note], ignore_index=True)
except Exception as e:
    print("loose-gate note skipped:", e)
cell_level_note = pd.DataFrame([
  {"gene":"__previous_cell_level_p__","set":"summary","deseq2_log2FC_diab_vs_ctrl":np.nan,
   "deseq2_pvalue":np.nan,"deseq2_padj":np.nan,"mean_logCPM_diab":np.nan,"mean_logCPM_ctrl":np.nan,
   "delta_logCPM":np.nan,"all_diab_above_all_ctrl":np.nan,
   "cell_MannWhitney_p_INVALID":"3.8e-6 (cell-level, pseudoreplication; NOT reproduced by us)",
   "cohens_d_cell_descr":np.nan}])
meta_rows = pd.concat([meta_rows, cell_level_note], ignore_index=True)
out = pd.concat([out, meta_rows], ignore_index=True)
out.to_csv(f"{OUT}/05_nephropathy.csv", index=False)
print(f"\nsaved -> {OUT}/05_nephropathy.csv")
print("\nDONE")
