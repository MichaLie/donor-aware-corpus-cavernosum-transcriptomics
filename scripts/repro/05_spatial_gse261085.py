#!/usr/bin/env python3
"""
05_spatial -- Spatial perivascular localization of the non-stress pericyte core in
GSE261085 (normal human corpus cavernosum spatial transcriptome, single NORMAL donor).

DATA (re-fetched from GEO FTP 2026-06-09):
  GSE261085 / GSM8134145 "Human, CC, Normal, 1" (HCCNC1).
  Platform: BMKMANU S1000 array-based spatial (NOT 10x Visium); PMID 39299236.
  GSM8134145_HCCNC1.expression.tar.gz contains 8 BIN levels (L1=finest 2.29M spots ...
  L13=coarsest 4,898 spots) of the SAME tissue cross-section (all cover x:[0,1000],
  y:[0,1005]; total UMIs ~4.08e7 conserved across levels -> they are bin aggregations,
  NOT separate regions). We analyze a Visium-like bin level (L7: 18,095 spots,
  median ~2,262 UMI / ~1,176 genes per spot) and replicate at L13 / L6 for robustness.

QUESTION: does the non-stress pericyte/perivascular core localize to perivascular /
endothelial-rich regions? We test both a simple spot-level correlation and a
deconvolution/enrichment approach, and report whether the conclusion is estimator-dependent.

NON-STRESS CORE (reframe orientation, up in diabetic ED pericytes; ERAP2 down):
  UP:   IL6, ADAMTS4, PLA2G2A, CYR61, PI15, SH3BP5
  DOWN: ERAP2
We score the core as a directional signature (down genes negated) AND as an unsigned
"pericyte-activation up-core" for transparency.

PERIVASCULAR reference = mural (pericyte/SMC) + endothelial markers:
  MURAL: RGS5, PDGFRB, NOTCH3, KCNJ8, HIGD1B, ACTA2, MYH11, TAGLN, CNN1, DES
  EC:    PECAM1, VWF, CLDN5, CDH5, FLT1, EGFL7

DISCIPLINE / PSEUDOREPLICATION GUARD:
  This is ONE normal donor. A spot is NOT a donor. We therefore treat spot-level
  correlations / enrichment as DESCRIPTIVE spatial structure within one sample, and we
  do NOT claim donor-level / disease inference from any spot p-value. No disease
  contrast is possible (no diabetic spatial sample here). The deliverable is purely
  whether the two estimators AGREE on perivascular localization within this one normal
  sample (the estimator-switch question).

Outputs: results/tables/repro/05_spatial.csv
"""
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path
from scipy.io import mmread
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
from scipy.optimize import nnls

ROOT = (os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
RAW  = f"{ROOT}/data/raw/GSE261085"
OUT  = f"{ROOT}/results/tables/repro"; os.makedirs(OUT, exist_ok=True)
sc.settings.verbosity = 1
RNG = np.random.default_rng(0)

CORE_UP   = ["IL6","ADAMTS4","PLA2G2A","CYR61","PI15","SH3BP5"]
CORE_DOWN = ["ERAP2"]
MURAL = ["RGS5","PDGFRB","NOTCH3","KCNJ8","HIGD1B","ACTA2","MYH11","TAGLN","CNN1","DES"]
EC    = ["PECAM1","VWF","CLDN5","CDH5","FLT1","EGFL7"]
PERIVASC = MURAL + EC

def load_level(level):
    d = f"{RAW}/{level}_heAuto"
    M = mmread(f"{d}/matrix.mtx.gz").tocsr()        # genes x spots
    feat = pd.read_csv(f"{d}/features.tsv.gz", sep="\t", header=None)
    bc   = pd.read_csv(f"{d}/barcodes.tsv.gz", header=None)[0].values
    pos  = pd.read_csv(f"{d}/barcodes_pos.tsv.gz", sep="\t", header=None,
                       names=["barcode","x","y"]).set_index("barcode")
    a = ad.AnnData(M.T.tocsr(), obs=pd.DataFrame(index=bc),
                   var=pd.DataFrame(index=feat[1].values))   # spots x genes
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    pos = pos.reindex(a.obs_names)
    a.obs["x"] = pos["x"].values; a.obs["y"] = pos["y"].values
    return a

def qc_norm(a):
    sc.pp.calculate_qc_metrics(a, inplace=True, percent_top=None)
    # keep real tissue spots: at least 200 genes & 300 UMI
    a = a[(a.obs.n_genes_by_counts >= 200) & (a.obs.total_counts >= 300)].copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a

def score(a, genes, name):
    g = [x for x in genes if x in a.var_names]
    sc.tl.score_genes(a, g, score_name=name, ctrl_size=50, random_state=0)
    return g

# ---- build a full-cell-type reference (for NNLS deconvolution) from GSE206528 ----
REF_CELLTYPES = ["Fibroblast","Smooth_muscle","Endothelial","Pericyte","Macrophage",
                 "T_cell","B_or_plasma_cell","Mast_cell","Schwann_or_neural","Epithelial"]
PERIVASC_CT   = ["Pericyte","Smooth_muscle","Endothelial"]   # perivascular cell types

def build_reference():
    """Mean log-norm profile per cell type from GSE206528 (cells -> argmax marker score)."""
    a = sc.read_h5ad("./data/processed/"
                     "GSE206528/GSE206528_scanpy_preprocessed.h5ad")
    scols = [f"score_{c}" for c in REF_CELLTYPES]
    lab = a.obs[scols].values.argmax(1)
    a.obs["ct"] = [REF_CELLTYPES[i] for i in lab]
    X = a.X  # log-norm
    prof = {}
    for c in REF_CELLTYPES:
        m = (a.obs["ct"].values == c)
        prof[c] = np.asarray(X[m].mean(0)).ravel()
    ref = pd.DataFrame(prof, index=a.var_names)   # genes x celltypes
    return ref

def signed_core_score(a):
    # directional: mean(up log-norm) - mean(down log-norm), z within sample
    up  = [g for g in CORE_UP   if g in a.var_names]
    dn  = [g for g in CORE_DOWN if g in a.var_names]
    Xu = np.asarray(a[:, up].X.todense())
    Xd = np.asarray(a[:, dn].X.todense()) if dn else np.zeros((a.n_obs,1))
    # z-score each gene across spots, then signed mean
    def zc(X):
        m = X.mean(0); s = X.std(0); s[s==0]=1
        return (X-m)/s
    su = zc(Xu).mean(1)
    sd = zc(Xd).mean(1) if dn else np.zeros(a.n_obs)
    return su - sd, up, dn

def analyze(level, REF):
    print(f"\n========== {level} ==========")
    a = load_level(level)
    print(f"loaded spots x genes = {a.shape}; coord x:[{a.obs.x.min():.1f},{a.obs.x.max():.1f}] "
          f"y:[{a.obs.y.min():.1f},{a.obs.y.max():.1f}]")
    a = qc_norm(a)
    print(f"after QC (>=200 genes & >=300 UMI): {a.n_obs} spots; "
          f"median UMI={np.median(a.obs.total_counts):.0f}, median genes={np.median(a.obs.n_genes_by_counts):.0f}")

    # marker scores
    score(a, MURAL, "sc_mural"); score(a, EC, "sc_ec"); score(a, PERIVASC, "sc_perivasc")
    score(a, CORE_UP, "sc_core_up")     # unsigned up-core
    a.obs["sc_core_signed"], up, dn = signed_core_score(a)

    rows = []
    n = a.n_obs
    # ---------- (a) SIMPLE spot-level correlations ----------
    for cscore in ["sc_core_signed","sc_core_up"]:
        for mscore,mlab in [("sc_perivasc","perivascular(mural+EC)"),
                            ("sc_mural","mural"),("sc_ec","endothelial")]:
            cv = a.obs[cscore].astype(np.float64).values
            mv = a.obs[mscore].astype(np.float64).values
            rho,prho = spearmanr(cv, mv)
            r,pr     = pearsonr(cv, mv)
            rows.append(dict(level=level, n_spots=n, approach="(a) simple spot correlation",
                             core=cscore, reference=mlab,
                             spearman_rho=round(rho,4), spearman_p=f"{prho:.2e}",
                             pearson_r=round(r,4), pearson_p=f"{pr:.2e}",
                             stat=round(rho,4),
                             conclusion=("perivascular-enriched" if rho>0.1 else
                                         ("anti/none" if rho<=0.1 else "weak"))))

    # ---------- (b) ENRICHMENT / deconvolution-style ----------
    # Define "perivascular-rich" spots = top quartile of perivascular marker score
    # (a coarse, transparent deconvolution surrogate: high mural+EC content).
    q = a.obs["sc_perivasc"].quantile(0.75)
    periv_mask = (a.obs["sc_perivasc"] >= q).values
    for cscore in ["sc_core_signed","sc_core_up"]:
        hi = a.obs[cscore].values[periv_mask]
        lo = a.obs[cscore].values[~periv_mask]
        U,pmw = mannwhitneyu(hi, lo, alternative="greater")
        # Cliff's delta effect size
        cliff = 2*U/(len(hi)*len(lo)) - 1
        rows.append(dict(level=level, n_spots=n,
                         approach="(b) enrichment: core in perivasc-rich(top25%) vs rest",
                         core=cscore, reference="perivascular(mural+EC) top quartile",
                         spearman_rho=np.nan, spearman_p="",
                         pearson_r=np.nan, pearson_p=f"{pmw:.2e}",
                         stat=round(float(cliff),4),
                         mean_hi=round(float(hi.mean()),4), mean_lo=round(float(lo.mean()),4),
                         conclusion=("perivascular-enriched" if cliff>0.1 else
                                     ("anti/none" if cliff<=0.1 else "weak"))))

    # ---------- (b2) gene-set deconvolution via cluster identity ----------
    # Leiden cluster the spots; identify the cluster MOST perivascular; test core there.
    av = a.copy()
    sc.pp.highly_variable_genes(av, n_top_genes=2000)
    av = av[:, av.var.highly_variable].copy()
    sc.pp.scale(av, max_value=10); sc.tl.pca(av, n_comps=30, random_state=0)
    sc.pp.neighbors(av, n_neighbors=15, random_state=0)
    sc.tl.leiden(av, resolution=0.5, random_state=0, flavor="igraph",
                 n_iterations=2, directed=False)
    a.obs["leiden"] = av.obs["leiden"].values
    cl = a.obs.groupby("leiden", observed=True)[["sc_perivasc","sc_core_signed","sc_core_up"]].mean()
    periv_cl = cl["sc_perivasc"].idxmax()
    cl_sorted0 = cl["sc_perivasc"].sort_values(ascending=False)
    margin0 = round(float(cl_sorted0.iloc[0] - cl_sorted0.iloc[1]), 4)
    in_mask = (a.obs["leiden"]==periv_cl).values
    for cscore in ["sc_core_signed","sc_core_up"]:
        hi = a.obs[cscore].values[in_mask]; lo = a.obs[cscore].values[~in_mask]
        U,pmw = mannwhitneyu(hi, lo, alternative="greater")
        cliff = 2*U/(len(hi)*len(lo)) - 1
        rows.append(dict(level=level, n_spots=n,
                         approach=f"(b2) deconv-by-cluster: core in most-perivasc Leiden cl{periv_cl} vs rest",
                         core=cscore, reference=f"perivasc Leiden cluster {periv_cl} (n={int(in_mask.sum())})",
                         spearman_rho=np.nan, spearman_p="",
                         pearson_r=np.nan, pearson_p=f"{pmw:.2e}",
                         stat=round(float(cliff),4),
                         mean_hi=round(float(hi.mean()),4), mean_lo=round(float(lo.mean()),4),
                         perivasc_cluster_margin=margin0,
                         conclusion=("perivascular-enriched (FRAGILE: tiny cluster margin)" if cliff>0.1 else
                                     ("anti/none" if cliff<=0.1 else "weak"))))
    print(f"\n[{level}] (b2) perivasc-cluster margin (top minus runner-up perivasc score) = {margin0:.3f}")

    # ---------- (c) NNLS DECONVOLUTION (the Tangram-analog) ----------
    # Per spot, fit non-negative cell-type proportions against a GSE206528 reference,
    # then test whether the core tracks the PERIVASCULAR proportion (Pericyte+SMC+EC).
    if REF is not None:
        genes = [g for g in REF.index if g in a.var_names]
        # use markery/variable-ish genes: union of marker sets + ref HVGs to keep it light
        Rsub = REF.loc[genes]
        # restrict to genes with some variance across cell types
        keepg = Rsub.index[(Rsub.std(1) > 0.05)]
        Rsub = Rsub.loc[keepg]
        Xs = np.asarray(a[:, list(Rsub.index)].X.todense(), dtype=np.float64)  # spots x genes (log-norm)
        Rm = Rsub.values.astype(np.float64)   # genes x celltypes
        props = np.zeros((Xs.shape[0], Rm.shape[1]))
        for i in range(Xs.shape[0]):
            sol,_ = nnls(Rm, Xs[i])
            s = sol.sum()
            props[i] = sol/s if s>0 else sol
        prop = pd.DataFrame(props, columns=Rsub.columns, index=a.obs_names)
        a.obs["prop_perivasc"] = prop[PERIVASC_CT].sum(1).values
        for cscore in ["sc_core_signed","sc_core_up"]:
            cv = a.obs[cscore].astype(np.float64).values
            pv = a.obs["prop_perivasc"].astype(np.float64).values
            rho,prho = spearmanr(cv, pv)
            r,pr     = pearsonr(cv, pv)
            rows.append(dict(level=level, n_spots=n,
                             approach="(c) NNLS deconv (Tangram-analog): core vs perivasc cell-type proportion",
                             core=cscore, reference="Pericyte+SMC+EC NNLS proportion",
                             spearman_rho=round(rho,4), spearman_p=f"{prho:.2e}",
                             pearson_r=round(r,4), pearson_p=f"{pr:.2e}",
                             stat=round(rho,4),
                             conclusion=("perivascular-enriched" if rho>0.1 else "anti/none")))
        print(f"[{level}] (c) NNLS mean perivasc proportion = {a.obs['prop_perivasc'].mean():.3f}; "
              f"signed-core~perivasc-prop spearman = "
              f"{spearmanr(a.obs['sc_core_signed'].astype(float), a.obs['prop_perivasc'].astype(float))[0]:.3f}")

    print("\nLeiden cluster means (perivasc / core):")
    print(cl.round(3).to_string())
    return rows, a

if __name__ == "__main__":
    print("building GSE206528 cell-type reference for NNLS deconvolution ...")
    try:
        REF = build_reference()
        print("reference cell types:", REF.columns.tolist(), "genes:", REF.shape[0])
    except Exception as e:
        print("reference build failed, NNLS deconv skipped:", e); REF = None
    all_rows = []
    for lvl in ["L7","L13","L6"]:
        rows, _ = analyze(lvl, REF)
        all_rows += rows
    df = pd.DataFrame(all_rows)
    df.to_csv(f"{OUT}/05_spatial.csv", index=False)
    pd.set_option("display.width", 200, "display.max_columns", 30)
    print("\n\n================ SUMMARY (saved -> 05_spatial.csv) ================")
    show = df[["level","approach","core","reference","stat","spearman_rho","pearson_p","conclusion"]]
    print(show.to_string(index=False))
