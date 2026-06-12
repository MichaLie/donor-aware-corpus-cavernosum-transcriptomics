#!/usr/bin/env python3
"""
Reproduce (independently, in pydeseq2) the donor-level pericyte-like pseudobulk
differential expression underlying the gene-tier reanalysis, using the archived
GSE206528 pericyte pseudobulk counts and annotation.

Goal: establish, within GSE206528 (no cross-dataset confound), which of the 12
signature genes are
  (a) shared-ED      : up in both diabetic and non-diabetic ED vs normal,
                       but NOT different between diabetic and non-diabetic ED
  (b) diabetes-associated: different in diabetic vs NON-diabetic ED.

Three contrasts on the 8 Pericyte_like donor pseudobulks (pydeseq2, ~condition):
  diab_ED vs normal      (compare to reference edgeR output as a method cross-check)
  nondiab_ED vs normal   (shared-ED signal)
  diab_ED vs nondiab_ED  (diabetic-versus-non-diabetic-ED contrast test)

Engine note: the primary analysis used edgeR-QL; pydeseq2 on the same
counts is an independent re-derivation and a DE-engine robustness check.
"""
import os
from pathlib import Path
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = (os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
CNT = f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
META = f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_metadata.tsv"
EDGER = f"{ROOT}/results/tables/GSE206528_mural_subtype_edger/all_mural_subtype_edgeR_results.tsv.gz"
OUTDIR = f"{ROOT}/results/tables/repro"
os.makedirs(OUTDIR, exist_ok=True)

SIG12 = ["IL6","ADAMTS4","PLA2G2A","CYR61","PI15","NR4A1","NR4A2","NR4A3","MT1A","RGS16","SH3BP5","ERAP2"]
STRESS5 = ["NR4A1","NR4A2","NR4A3","MT1A","RGS16"]   # immediate-early / stress subset

# ---- load pericyte_like pseudobulk (samples x genes) ----
counts = pd.read_csv(CNT, sep="\t", index_col=0)
meta_all = pd.read_csv(META, sep="\t")
peri = meta_all[meta_all.mural_subtype == "Pericyte_like"].copy()
group = {"diabetic_ED":"diab_ED", "non_diabetic_ED":"nondiab_ED", "normal":"normal"}
peri["condition"] = peri["ed_subgroup"].map(group)
peri = peri.set_index("pseudobulk_id")
samples = [s for s in counts.columns if s in peri.index]            # 8 pericyte_like donors
cnt = counts[samples].T.astype(int)                                 # samples x genes
md = peri.loc[samples, ["donor_id","condition","n_cells"]].copy()
print("Pericyte_like donor pseudobulks:")
print(md.to_string()); print()

# ---- filter genes (keep >=10 total counts; force-keep the 12) ----
keep = (cnt.sum(0) >= 10) | cnt.columns.isin(SIG12)
cnt = cnt.loc[:, keep]
print(f"genes after filter: {cnt.shape[1]} (of {counts.shape[0]})")

md["condition"] = pd.Categorical(md["condition"], categories=["normal","nondiab_ED","diab_ED"])
dds = DeseqDataSet(counts=cnt, metadata=md, design="~condition", refit_cooks=False, quiet=True)
dds.deseq2()

def contrast(a, b):
    st = DeseqStats(dds, contrast=["condition", a, b], quiet=True)
    st.summary()
    r = st.results_df.loc[SIG12, ["log2FoldChange","pvalue","padj"]].copy()
    r.columns = [f"lfc_{a}_vs_{b}", f"p_{a}_vs_{b}", f"padj_{a}_vs_{b}"]
    return r

c_dn  = contrast("diab_ED","normal")        # diabetic vs normal
c_nn  = contrast("nondiab_ED","normal")     # shared-ED vs normal
c_dne = contrast("diab_ED","nondiab_ED")    # diabetic-versus-non-diabetic-ED contrast (NEW)
res = pd.concat([c_dn, c_nn, c_dne], axis=1).loc[SIG12]

# ---- reference edgeR cross-check (diabetic_ED_vs_normal) ----
ed = pd.read_csv(EDGER, sep="\t")
edp = ed[(ed.mural_subtype=="Pericyte_like") & (ed.comparison=="diabetic_ED_vs_normal")]
edp = edp.set_index("gene_symbol").reindex(SIG12)[["logFC","FDR"]]
edp.columns = ["edgeR_logFC_DvN","edgeR_FDR_DvN"]

# ---- classify ----
def cls(g):
    r = res.loc[g]
    up_nondiab = (r["lfc_nondiab_ED_vs_normal"] > 0.5) and (r["p_nondiab_ED_vs_normal"] < 0.05)
    diff_dvsne = (r["p_diab_ED_vs_nondiab_ED"] < 0.05)
    if diff_dvsne:
        return "diabetes-associated"
    if up_nondiab:
        return "shared-ED"
    return "ambiguous"
res["class"] = [cls(g) for g in res.index]
out = res.join(edp)

pd.set_option("display.width", 240, "display.max_columns", 40, "display.float_format", lambda x: f"{x:.3g}")
print("\n===== 12 signature genes: 3 contrasts (pydeseq2) + edgeR cross-check =====")
show = out[["lfc_diab_ED_vs_normal","padj_diab_ED_vs_normal","edgeR_logFC_DvN","edgeR_FDR_DvN",
            "lfc_nondiab_ED_vs_normal","padj_nondiab_ED_vs_normal",
            "lfc_diab_ED_vs_nondiab_ED","p_diab_ED_vs_nondiab_ED","class"]]
print(show.to_string())
out.to_csv(f"{OUTDIR}/01_pericyte_3contrast_DE_GSE206528.csv", float_format="%.12g")
print(f"\nsaved -> {OUTDIR}/01_pericyte_3contrast_DE_GSE206528.csv")

# correlation of pydeseq2 vs edgeR LFC (diabetic vs normal) = engine-robustness check
m = out[["lfc_diab_ED_vs_normal","edgeR_logFC_DvN"]].dropna()
print(f"\npydeseq2 vs edgeR LFC (diab vs normal), Pearson r = {m.iloc[:,0].corr(m.iloc[:,1]):.3f} (n={len(m)})")
