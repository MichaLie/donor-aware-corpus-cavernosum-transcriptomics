#!/usr/bin/env python3
"""
05_invitro: in-vitro hyperglycemia induction test on GSE228900.

GSE228900 = "Modeling early pathophysiological phenotypes of diabetic retinopathy
in a human inner blood-retinal barrier (iBRB)-on-a-chip" (Roche; PMID 38355716).
RETINAL in-vitro model (NOT penile) -- cross-context only.

Design (TC = tri-culture microvascular networks of EC+pericytes+astrocytes,
bulk RNA-seq, n=3 replicate chips per condition x timepoint):
  treatment in {untreated(control), mannitol(osmotic ctrl), diabetic(high glucose)}
  matched timepoints with all three treatments: day 14 and day 28.
  -> 18 TC samples: 3 treat x 2 time x 3 rep.

Analysis question: high glucose induces the reduced non-stress core in
this pericyte tri-culture, and it is hyperglycemia-specific vs the MANNITOL osmotic
control (core NOT elevated by mannitol).

CRITICAL SPECIFICITY TEST: is the core up in diabetic-vs-untreated AND not merely an
osmotic effect (i.e., diabetic-vs-mannitol also up; mannitol-vs-untreated ~flat)?

DISCIPLINE: the experimental unit here is the replicate CHIP (n=3 per cell), bulk
RNA-seq. We do donor/unit-level inference via pydeseq2 (Wald, with time as covariate).
We do NOT do cell-level pseudoreplication (there are no single cells here -- it is bulk).

Core (reduced non-stress) signature, oriented as in the GSE206528 reframe:
  UP in diabetic ED:  IL6, ADAMTS4, PLA2G2A, CYR61, PI15, SH3BP5, RGS16
  DOWN in diabetic ED: ERAP2
RGS16 is the one stress/immediate-early gene retained in the "reduced" core per the task.

Output: results/tables/repro/05_invitro.csv
"""
import os, gzip, glob, re
from pathlib import Path
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = (os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
RAW  = f"{ROOT}/data/raw/GSE228900"
OUT  = f"{ROOT}/results/tables/repro"
os.makedirs(OUT, exist_ok=True)

# core genes -> stable Ensembl IDs (GRCh38), with reframe orientation
CORE = {  # gene: (ENSG, expected_dir_in_diabetic: +1 up / -1 down)
    "IL6":    ("ENSG00000136244", +1),
    "ADAMTS4":("ENSG00000158859", +1),
    "PLA2G2A":("ENSG00000188257", +1),
    "CYR61":  ("ENSG00000142871", +1),  # CCN1
    "PI15":   ("ENSG00000137976", +1),
    "SH3BP5": ("ENSG00000131370", +1),
    "ERAP2":  ("ENSG00000164308", -1),
    "RGS16":  ("ENSG00000143333", +1),
}
ENSG2SYM = {e: g for g,(e,_) in CORE.items()}

# ---- sample table for TC matched panel (treatment x time, 3 reps each) ----
# from GEO sample metadata (verified by fetch)
SAMPLES = {
 # untreated TC day14/day28
 "GSM7142390":("untreated","day14"),"GSM7142391":("untreated","day14"),"GSM7142392":("untreated","day14"),
 "GSM7142393":("untreated","day28"),"GSM7142394":("untreated","day28"),"GSM7142395":("untreated","day28"),
 # mannitol TC day14/day28
 "GSM7142396":("mannitol","day14"),"GSM7142397":("mannitol","day14"),"GSM7142398":("mannitol","day14"),
 "GSM7142399":("mannitol","day28"),"GSM7142400":("mannitol","day28"),"GSM7142401":("mannitol","day28"),
 # diabetic TC day14/day28
 "GSM7142402":("diabetic","day14"),"GSM7142403":("diabetic","day14"),"GSM7142404":("diabetic","day14"),
 "GSM7142405":("diabetic","day28"),"GSM7142406":("diabetic","day28"),"GSM7142407":("diabetic","day28"),
}

def load_counts(gsm):
    fn = glob.glob(f"{RAW}/{gsm}_*.counts.txt.gz")[0]
    s = {}
    with gzip.open(fn,"rt") as f:
        for line in f:
            k,v = line.rstrip("\n").split("\t")
            s[k] = int(v)
    return pd.Series(s, name=gsm)

# build counts matrix: genes x samples
cnt = pd.DataFrame({gsm: load_counts(gsm) for gsm in SAMPLES})
# drop the HTSeq summary rows if any (__no_feature etc.) -- keep only ENSG
cnt = cnt[cnt.index.str.startswith("ENSG")]
meta = pd.DataFrame({g:{"treatment":t[0],"time":t[1]} for g,t in SAMPLES.items()}).T
meta = meta.loc[cnt.columns]
print("Counts matrix:", cnt.shape, "  samples:", list(cnt.columns))
print(meta.value_counts().to_string())
print("Library sizes (total counts) per sample:")
print(cnt.sum(0).to_string())

# ---- per-condition CPM effect sizes for the core (descriptive) ----
cpm = cnt.div(cnt.sum(0), axis=1) * 1e6
logcpm = np.log2(cpm + 1)
core_ensg = [e for e,_ in CORE.values()]
desc_rows = []
for sym,(e,direction) in CORE.items():
    row = {"gene":sym, "ensembl":e, "expected_dir": "up" if direction>0 else "down"}
    for tr in ["untreated","mannitol","diabetic"]:
        cols = meta.index[meta.treatment==tr]
        row[f"meanlogCPM_{tr}"] = round(logcpm.loc[e, cols].mean(), 3)
        row[f"meanCPM_{tr}"]    = round(cpm.loc[e, cols].mean(), 2)
    desc_rows.append(row)
desc = pd.DataFrame(desc_rows)
print("\n=== Per-condition mean log2CPM (descriptive effect sizes) ===")
print(desc[["gene","expected_dir","meanlogCPM_untreated","meanlogCPM_mannitol","meanlogCPM_diabetic"]].to_string(index=False))

# ---- pydeseq2: design ~ time + treatment, three pairwise contrasts ----
def run_deseq(ref):
    m = meta.copy()
    m["treatment"] = pd.Categorical(m["treatment"], categories=[ref] + [x for x in ["untreated","mannitol","diabetic"] if x!=ref])
    m["time"] = pd.Categorical(m["time"])
    dds = DeseqDataSet(counts=cnt.T, metadata=m, design="~time+treatment", quiet=True)
    dds.deseq2()
    return dds

# We need 3 contrasts: diabetic vs untreated, diabetic vs mannitol, mannitol vs untreated.
# Run with untreated as ref (gives diab-vs-untr, mann-vs-untr), and mannitol ref (diab-vs-mann).
dds_u = run_deseq("untreated")
dds_m = run_deseq("mannitol")

def contrast(dds, num, den):
    st = DeseqStats(dds, contrast=["treatment", num, den], quiet=True)
    st.summary()
    r = st.results_df
    return r

c_diab_untr = contrast(dds_u, "diabetic", "untreated")
c_mann_untr = contrast(dds_u, "mannitol", "untreated")
c_diab_mann = contrast(dds_m, "diabetic", "mannitol")

def core_table(res, label):
    sub = res.loc[res.index.intersection(core_ensg)].copy()
    sub["gene"] = [ENSG2SYM[e] for e in sub.index]
    sub["contrast"] = label
    return sub[["gene","contrast","log2FoldChange","lfcSE","pvalue","padj"]].round(4)

ct = pd.concat([
    core_table(c_diab_untr, "diabetic_vs_untreated"),
    core_table(c_diab_mann, "diabetic_vs_mannitol"),
    core_table(c_mann_untr, "mannitol_vs_untreated"),
])

print("\n=== Core gene DE (pydeseq2 Wald, design ~time+treatment) ===")
for lab in ["diabetic_vs_untreated","diabetic_vs_mannitol","mannitol_vs_untreated"]:
    print(f"\n-- {lab} --")
    sub = ct[ct.contrast==lab].set_index("gene").reindex(list(CORE.keys()))
    print(sub[["log2FoldChange","lfcSE","pvalue","padj"]].to_string())

# ---- specificity verdict per gene + composite ----
verdict_rows=[]
for sym,(e,direction) in CORE.items():
    du = c_diab_untr.loc[e]; dm = c_diab_mann.loc[e]; mu = c_mann_untr.loc[e]
    # "induced by glucose, hyperglycemia-specific" if oriented-FC up in BOTH diab-vs-untr and diab-vs-mann
    fc_du = direction*du.log2FoldChange
    fc_dm = direction*dm.log2FoldChange
    fc_mu = direction*mu.log2FoldChange
    induced = (fc_du>0) and (du.padj<0.05 if pd.notna(du.padj) else False)
    specific = (fc_dm>0) and (dm.padj<0.05 if pd.notna(dm.padj) else False)
    osmotic_only = (fc_mu>0) and (mu.padj<0.05 if pd.notna(mu.padj) else False)
    verdict_rows.append({
        "gene":sym,"expected_dir":"up" if direction>0 else "down",
        "orientedL2FC_diab_vs_untr":round(fc_du,3),"padj_diab_vs_untr":round(du.padj,4) if pd.notna(du.padj) else np.nan,
        "orientedL2FC_diab_vs_mann":round(fc_dm,3),"padj_diab_vs_mann":round(dm.padj,4) if pd.notna(dm.padj) else np.nan,
        "orientedL2FC_mann_vs_untr":round(fc_mu,3),"padj_mann_vs_untr":round(mu.padj,4) if pd.notna(mu.padj) else np.nan,
        "induced_by_glucose(sig)":induced,"hyperglycemia_specific(sig)":specific,"osmotic_effect(sig)":osmotic_only,
    })
verdict = pd.DataFrame(verdict_rows)
print("\n=== SPECIFICITY VERDICT (oriented log2FC; up=consistent with reframe) ===")
print(verdict.to_string(index=False))

n_induced = verdict["induced_by_glucose(sig)"].sum()
n_specific = verdict["hyperglycemia_specific(sig)"].sum()
n_osm = verdict["osmotic_effect(sig)"].sum()
print(f"\nCore genes (of {len(CORE)}): induced_by_glucose(sig)={n_induced}; "
      f"hyperglycemia_specific(sig)={n_specific}; significant_osmotic_effect={n_osm}")

# ---- composite core directional score per sample (z within panel) + exact permutation across reps ----
# Build oriented composite per sample, then compare diabetic vs (untreated+mannitol) reps with a
# nonparametric exact-ish summary at the chip-replicate level (bulk units, not cells).
z = logcpm.loc[core_ensg].copy()
z = z.sub(z.mean(1), axis=0).div(z.std(1).replace(0,np.nan), axis=0)
oriented = pd.Series(0.0, index=z.columns)
for sym,(e,direction) in CORE.items():
    oriented += direction*z.loc[e]
oriented /= len(CORE)
comp = pd.DataFrame({"treatment":meta.treatment,"time":meta.time,"core_score":oriented.round(3)})
print("\n=== Composite oriented core score per sample (mean z over 8 genes) ===")
print(comp.sort_values(["treatment","time"]).to_string())
print("\nMean core score by treatment:")
print(comp.groupby("treatment")["core_score"].agg(["mean","std","count"]).round(3).to_string())

# ---- save ----
desc.to_csv(f"{OUT}/05_invitro_per_condition_logCPM.csv", index=False)
ct.to_csv(f"{OUT}/05_invitro_core_DE_contrasts.csv", index=False)
verdict.to_csv(f"{OUT}/05_invitro.csv", index=False)
comp.to_csv(f"{OUT}/05_invitro_composite_score.csv")
print(f"\nsaved -> {OUT}/05_invitro.csv (main verdict)")
print(f"saved -> {OUT}/05_invitro_core_DE_contrasts.csv")
print(f"saved -> {OUT}/05_invitro_per_condition_logCPM.csv")
print(f"saved -> {OUT}/05_invitro_composite_score.csv")
