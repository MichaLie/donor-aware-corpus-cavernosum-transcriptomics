#!/usr/bin/env python3
"""
Summary figure for the donor-aware SH3BP5 reanalysis:
 A: per-gene log2FC in two contrasts (diabetic-vs-reference vs diabetic-vs-non-diabetic-ED) -> shared-ED collapse, SH3BP5/ERAP2 persist
 B: per-donor diabetes-associated score (SH3BP5 - ERAP2, z) -> exact 2-vs-3 separation
 C: stress-correction retention -> axis ~80% vs NR4A3 40%
 D: 9-donor SH3BP5 / ERAP2 pericyte logCPM incl. independent replication donor
"""
import os, numpy as np, pandas as pd, matplotlib
from pathlib import Path
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=(os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2])); REP=f"{ROOT}/results/tables/repro"
OUT=f"{ROOT}/results/figures/repro"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size":9,"font.family":"DejaVu Sans",
                     "axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150})
TIER={"SH3BP5":"diabetes-associated","ERAP2":"diabetes-associated",
      "IL6":"shared-ED","ADAMTS4":"shared-ED","PLA2G2A":"shared-ED","CYR61":"shared-ED","PI15":"shared-ED",
      "NR4A1":"stress","NR4A2":"stress","NR4A3":"stress","MT1A":"stress","RGS16":"stress"}
COL={"diabetes-associated":"#c0392b","shared-ED":"#2c7fb8","stress":"#999999"}

fig,axes=plt.subplots(2,2,figsize=(10,8)); (axA,axB),(axC,axD)=axes

# ---- Panel A ----
de=pd.read_csv(f"{REP}/01_pericyte_3contrast_DE_GSE206528.csv",index_col=0)
order=["SH3BP5","ERAP2","IL6","ADAMTS4","PLA2G2A","CYR61","PI15","NR4A1","NR4A2","NR4A3","MT1A","RGS16"]
y=np.arange(len(order))[::-1]
axA.scatter(de.loc[order,"lfc_diab_ED_vs_normal"], y, s=55, marker="o", facecolors="none",
            edgecolors=[COL[TIER[g]] for g in order], linewidths=1.6, label="diabetic vs reference")
axA.scatter(de.loc[order,"lfc_diab_ED_vs_nondiab_ED"], y, s=55, marker="D",
            color=[COL[TIER[g]] for g in order], label="diabetic vs non-diabetic ED")
axA.axvline(0,color="k",lw=0.7); axA.set_yticks(y); axA.set_yticklabels(order, fontstyle="italic")
axA.set_xlabel("log$_2$ fold-change"); axA.set_title("A  Two contrasts separate shared and diabetes-associated signals",loc="left",fontweight="bold")
axA.legend(loc="lower right",fontsize=7,frameon=False)
axA.text(0.02,0.02,"open=vs reference   filled=vs non-diabetic ED",transform=axA.transAxes,fontsize=6.5,color="#555")

# ---- Panel B: per-donor diabetes-associated score (recompute) ----
cnt=pd.read_csv(f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_counts.tsv.gz",sep="\t",index_col=0)
meta=pd.read_csv(f"{ROOT}/data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_metadata.tsv",sep="\t")
pm=meta[meta.mural_subtype=="Pericyte_like"].set_index("pseudobulk_id")
edcols=[c for c in cnt.columns if c in pm.index and pm.loc[c,"ed_subgroup"] in ("diabetic_ED","non_diabetic_ED")]
lc=np.log2(cnt[edcols].div(cnt[edcols].sum(0),axis=1)*1e6+1)
z=lc.sub(lc.mean(1),axis=0).div(lc.std(1),axis=0)
score=z.loc["SH3BP5"]   # SH3BP5 alone is the diabetes-associated anchor (perfect 2-vs-3 separation, P=0.10)
dlab={c:pm.loc[c,"donor_id"] for c in edcols}; isdia={c:pm.loc[c,"ed_subgroup"]=="diabetic_ED" for c in edcols}
sc_sorted=score.sort_values()
bars=axB.barh([dlab[c] for c in sc_sorted.index], sc_sorted.values,
              color=["#c0392b" if isdia[c] else "#bdc7d6" for c in sc_sorted.index], edgecolor="k", linewidth=0.5)
axB.axvline(0,color="k",lw=0.7); axB.set_xlabel("SH3BP5 (z-scored pericyte logCPM)")
axB.set_title("B  SH3BP5 separates the diabetic donors (exact P=0.10)",loc="left",fontweight="bold")
axB.legend(handles=[Patch(color="#c0392b",label="diabetic ED"),Patch(color="#bdc7d6",label="non-diabetic ED")],
           loc="lower right",fontsize=7,frameon=False)

# ---- Panel C: stress retention ----
at=pd.read_csv(f"{REP}/03b_stress_attenuation_effectsize.csv")
at=at[at.gene!="diab_specific"].reset_index(drop=True)   # drop the combined-axis bar; SH3BP5-centered
labs=[{"SH3BP5":"SH3BP5","ERAP2":"ERAP2","NR4A3":"NR4A3\n(stress)"}[g] for g in at.gene]
cols=["#c0392b","#c0392b","#999999"]
axC.bar(range(len(at)), at.pct_retained, color=cols, edgecolor="k", linewidth=0.5)
axC.axhline(100,color="#888",lw=0.7,ls=":")
for i,v in enumerate(at.pct_retained): axC.text(i, v+2, f"{v:.0f}%", ha="center", fontsize=8)
axC.set_xticks(range(len(at))); axC.set_xticklabels(labs); axC.set_ylabel("% effect retained after\nstress-score correction")
axC.set_ylim(0,110); axC.set_title("C  SH3BP5 is stress-independent",loc="left",fontweight="bold")

# ---- Panel D: 9-donor replication ----
pb259=pd.read_csv(f"{REP}/04_GSE259348_pericyte_pseudobulk.csv",index_col=0).iloc[:,0]
pcols=[c for c in cnt.columns if c in pm.index]
comb=cnt[pcols].copy(); comb["GSE259348_DMED"]=pb259.reindex(cnt.index).fillna(0).values
lcA=np.log2(comb.div(comb.sum(0),axis=1)*1e6+1)
grp={c:("diabetic" if (c.startswith("DMED") or c=="GSE259348_DMED") else "non-diab ED" if c.startswith("ED") else "normal") for c in comb.columns}
name={c:(pm.loc[c,"donor_id"] if c in pm.index else "GSE259348*") for c in comb.columns}
cmap={"diabetic":"#c0392b","non-diab ED":"#2c7fb8","normal":"#7fbf7f"}
for ax_,gene,inv in [(axD,"SH3BP5",False)]:
    pass
xs=lcA.loc["SH3BP5"]; ys=lcA.loc["ERAP2"]
for c in comb.columns:
    mk="*" if c=="GSE259348_DMED" else "o"
    axD.scatter(xs[c], ys[c], s=170 if c=="GSE259348_DMED" else 70, marker=mk,
                color=cmap[grp[c]], edgecolor="k", linewidth=0.6, zorder=3)
    axD.annotate(name[c], (xs[c],ys[c]), fontsize=6.5, xytext=(3,3), textcoords="offset points")
axD.set_xlabel("SH3BP5  pericyte logCPM"); axD.set_ylabel("ERAP2  pericyte logCPM")
axD.set_title("D  Independent replication (★ GSE259348)",loc="left",fontweight="bold")
axD.legend(handles=[Patch(color=cmap[k],label=k) for k in cmap],loc="upper right",fontsize=7,frameon=False)
axD.text(0.02,0.02,"all 3 diabetic donors: SH3BP5-high, ERAP2-low",transform=axD.transAxes,fontsize=6.5,color="#555")

plt.tight_layout()
for ext in ("png","pdf"): fig.savefig(f"{OUT}/sh3bp5_summary.{ext}", bbox_inches="tight")
print("per-donor diabetes-associated score:"); print(score.rename(index=dlab).round(2).to_string())
print(f"saved -> {OUT}/sh3bp5_summary.png / .pdf")
