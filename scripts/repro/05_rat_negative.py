#!/usr/bin/env python3
"""
P05 — rat STZ-DMED corpus cavernosum (GSE259299, BD Rhapsody RSEC) negative control.

Goal (honest negative): the human pericyte "shared-ED activation" signature does NOT
conserve / inverts in diabetic-rat corpus-cavernosum pericytes, AND the rat design is
depth/platform-confounded (DMED samples much shallower; BD Rhapsody not 10x). We test the
DONOR-level direction of the core in DMED (4) vs control (2) pericyte pseudobulks. We NEVER
report a cell-level p as inference (that is a pseudoreplication failure mode);
cell-level numbers are descriptive only. Inference = exact 4-vs-2 label permutation on the
donor-level signature score.

Data are re-fetched RSEC MolsPerCell CSVs from GEO (GSE259299_RAW.tar), reference
Rattus_norvegicus.mRatBN7.2.106 (per the BD pipeline header in each file). We read them
directly here; no external intermediate is used.

Orthology note (verified against the actual rat gene panel in the CSV header):
  Cyr61 -> Ccn1 (modern symbol; present)
  Sh3bp5 -> ABSENT in panel (only paralog Sh3bp5l) -> diabetes-associated axis NOT testable in rat
  Erap2  -> ABSENT in panel (only Erap1; rat lacks a 1:1 Erap2 ortholog) -> NOT testable
  Higd1b -> ABSENT (only Higd1a/c) -> dropped from pericyte marker set
"""
import os, glob, itertools, numpy as np, pandas as pd, scanpy as sc, anndata as ad
from pathlib import Path

ROOT = (os.environ.get("PROJECT_ROOT") or str(Path(__file__).resolve().parents[2]))
RAW  = f"{ROOT}/data/raw/GSE259299"
OUT  = f"{ROOT}/results/tables/repro"
os.makedirs(OUT, exist_ok=True)

# sample -> condition (from GEO series: N control, H STZ-diabetic ED)
SAMPLES = {
    "GSM8113509_N1": "control", "GSM8113510_N2": "control",
    "GSM8113511_H1": "DMED",    "GSM8113512_H2": "DMED",
    "GSM8113513_H3": "DMED",    "GSM8113514_H4": "DMED",
}

# ---- rat orthologs of the human core (human -> rat symbol present in panel) ----
# shared-ED activation core + stress gene Rgs16 (all map); diabetes axis (Sh3bp5/Erap2) absent.
CORE_MAP = {"IL6":"Il6","ADAMTS4":"Adamts4","PLA2G2A":"Pla2g2a",
            "CYR61":"Ccn1","PI15":"Pi15","RGS16":"Rgs16"}
ABSENT   = {"SH3BP5":"absent (only paralog Sh3bp5l)","ERAP2":"absent (only Erap1; no 1:1 rat ortholog)"}
SHARED_CORE = ["Il6","Adamts4","Pla2g2a","Ccn1","Pi15"]   # the 5 shared-ED activation genes

PERI_MARK = ["Rgs5","Pdgfrb","Notch3","Kcnj8"]   # Higd1b absent in panel
SMC_MARK  = ["Acta2","Myh11","Tagln","Cnn1","Des"]


def read_rsec(path):
    """Read one BD Rhapsody RSEC MolsPerCell CSV (## header lines) -> AnnData (cells x genes)."""
    # count comment lines to skip
    skip = 0
    import gzip
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                skip += 1
            else:
                break
    df = pd.read_csv(path, skiprows=skip, index_col=0)   # first col = Cell_Index
    df.index = df.index.astype(str)
    a = ad.AnnData(X=df.values.astype(np.float32),
                   obs=pd.DataFrame(index=df.index),
                   var=pd.DataFrame(index=df.columns))
    return a


def annotate_mural(a):
    """Light QC + leiden; score pericyte vs SMC; return pericyte-cell mask + descriptive stats."""
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    sc.pp.filter_cells(a, min_genes=100)
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.raw = a
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    av = a[:, a.var.highly_variable].copy()
    sc.pp.scale(av, max_value=10)
    sc.tl.pca(av, n_comps=30, random_state=0)
    sc.pp.neighbors(av, n_neighbors=15, random_state=0)
    sc.tl.leiden(av, resolution=1.0, random_state=0, flavor="igraph",
                 n_iterations=2, directed=False)
    a.obs["leiden"] = av.obs["leiden"].values
    pm = [g for g in PERI_MARK if g in a.var_names]
    sm = [g for g in SMC_MARK if g in a.var_names]
    sc.tl.score_genes(a, pm, score_name="sc_Pericyte")
    sc.tl.score_genes(a, sm, score_name="sc_SMC")
    # Mural compartment of corpus cavernosum is overwhelmingly SMC; pericytes are the Rgs5-high
    # / Acta2-relatively-LOW subset. Robust 2-step call (works in every sample):
    #   1) mural clusters = leiden clusters with mean Rgs5 > 0.5 (clearly expressed pericyte marker)
    #   2) pericyte cluster(s) = mural clusters whose mean Acta2 is below the mural-cluster median
    #      (the pericyte-leaning, less-contractile mural cells), requiring Rgs5 > SMC-cluster Rgs5.
    obs = a.obs.copy()
    for g in ["Rgs5", "Pdgfrb", "Acta2", "Myh11", "Tagln"]:
        obs[g] = np.asarray(a[:, g].X).ravel() if g in a.var_names else 0.0
    cl = obs.groupby("leiden", observed=True).agg(
        sc_Pericyte=("sc_Pericyte", "mean"), sc_SMC=("sc_SMC", "mean"),
        Rgs5=("Rgs5", "mean"), Pdgfrb=("Pdgfrb", "mean"), Acta2=("Acta2", "mean"),
        Myh11=("Myh11", "mean"), Tagln=("Tagln", "mean"))
    cl["n"] = obs.groupby("leiden", observed=True).size()
    mural = cl[cl["Rgs5"] > 0.5]
    peri_cl = []
    if len(mural):
        # pericyte = mural cluster(s) with the lowest Acta2 (most pericyte-like vs contractile SMC)
        amed = mural["Acta2"].median()
        cand = mural[mural["Acta2"] <= amed]
        # if all mural clusters tie, take the single lowest-Acta2 cluster
        if len(cand) == len(mural):
            cand = mural.nsmallest(1, "Acta2")
        peri_cl = list(cand.index)
    a.obs["is_peri"] = a.obs["leiden"].isin(peri_cl)
    a.uns["cl_table"] = cl
    return a, peri_cl, pm, sm


def main():
    rng = np.random.default_rng(0)
    rows_depth = []   # per-sample depth/QC for confound table
    pb = {}           # sample -> pericyte pseudobulk (gene -> summed counts)
    peri_desc = {}    # sample -> dict of pericyte marker means + n_peri

    for samp, cond in SAMPLES.items():
        path = glob.glob(f"{RAW}/{samp}_RSEC_MolsPerCell.csv.gz")[0]
        a = read_rsec(path)
        n_cells_raw = a.n_obs
        # raw depth metrics (before filtering) for the confound quantification
        counts = np.asarray(a.X)
        umis = counts.sum(1)
        genes = (counts > 0).sum(1)
        depth_row = dict(sample=samp, condition=cond, n_cells_raw=int(n_cells_raw),
                         median_umi=float(np.median(umis)),
                         median_genes=float(np.median(genes)))
        a, peri_cl, pm, sm = annotate_mural(a)
        peri = a[a.obs.is_peri]
        n_peri = peri.n_obs
        depth_row["n_cells_qc"] = int(a.n_obs)
        depth_row["n_pericyte"] = int(n_peri)
        depth_row["peri_clusters"] = ";".join(peri_cl)
        rows_depth.append(depth_row)
        # descriptive pericyte marker check (log-norm means) — NOT inference
        if n_peri > 0:
            peri_desc[samp] = {g: float(np.asarray(peri[:, g].X).mean()) for g in pm + sm if g in peri.var_names}
        # pericyte pseudobulk = sum of RAW counts over pericyte cells
        cmat = np.asarray(peri.layers["counts"]).sum(0).ravel()
        pb[samp] = pd.Series(cmat, index=peri.var_names)
        print(f"[{samp} {cond}] raw cells={n_cells_raw} qc={a.n_obs} peri={n_peri} "
              f"medUMI={depth_row['median_umi']:.0f} medGenes={depth_row['median_genes']:.0f} "
              f"periClusters={peri_cl}")

    depth = pd.DataFrame(rows_depth)
    print("\n=== depth / platform confound ===")
    print(depth.to_string(index=False))

    # ---- build pericyte pseudobulk matrix (genes x 6 samples), logCPM ----
    allg = sorted(set().union(*[set(s.index) for s in pb.values()]))
    M = pd.DataFrame({s: pb[s].reindex(allg).fillna(0) for s in SAMPLES}, index=allg)
    lib = M.sum(0)
    logcpm = np.log2(M.div(lib, axis=1) * 1e6 + 1)

    core_present = [g for g in SHARED_CORE if g in logcpm.index] + (["Rgs16"] if "Rgs16" in logcpm.index else [])
    print(f"\ncore genes present in rat pericyte pseudobulk: {core_present}")
    print(f"core genes ABSENT (no testable rat ortholog): {ABSENT}")

    # Only samples with a real pericyte pool (>=20 cells) enter the donor-level test;
    # samples with 0 pericytes are NOT imputed to 0 (that would fabricate signal).
    MIN_PERI = 20
    n_peri_by = depth.set_index("sample")["n_pericyte"]
    used = [s for s in SAMPLES if n_peri_by[s] >= MIN_PERI]
    ctrl = [s for s in used if SAMPLES[s] == "control"]
    dmed = [s for s in used if SAMPLES[s] == "DMED"]
    print(f"\nsamples with >= {MIN_PERI} pericytes used in donor test: "
          f"control={ctrl} DMED={dmed}")

    # ---- per-gene donor-level direction: mean logCPM DMED vs control ----
    g_rows = []
    for g in core_present:
        v = logcpm.loc[g]
        m_ctrl = v[ctrl].mean(); m_dmed = v[dmed].mean()
        g_rows.append(dict(human="".join([k for k, r in CORE_MAP.items() if r == g] or [g.upper()]),
                           rat=g, mean_logCPM_control=round(m_ctrl, 3),
                           mean_logCPM_DMED=round(m_dmed, 3),
                           delta_DMEDminusCtrl=round(m_dmed - m_ctrl, 3),
                           direction_up_in_DMED=bool(m_dmed > m_ctrl)))
    gtab = pd.DataFrame(g_rows)
    print("\n=== per-gene donor-level direction (rat pericyte pseudobulk, logCPM) ===")
    print(gtab.to_string(index=False))

    # human signature is UP in ED (shared-ED activation). Conservation => UP in DMED vs control.
    # signature score = mean logCPM of the 5 shared-ED core genes, per sample (z within gene to balance scales).
    core5 = [g for g in SHARED_CORE if g in logcpm.index]
    L = logcpm.loc[core5, used]   # z-score ONLY across the used samples
    z = L.sub(L.mean(1), axis=0).div(L.std(1).replace(0, np.nan), axis=0).fillna(0)
    score = z.mean(0)   # per-sample signature score (used samples only)
    s_ctrl = score[ctrl].mean(); s_dmed = score[dmed].mean()
    obs_diff = s_dmed - s_ctrl   # >0 = conserved (up in DMED); <0 = inverts
    print("\n=== per-sample shared-ED signature score (z of 5 core genes, used samples) ===")
    print(score.round(3).to_string())
    print(f"score: control mean={s_ctrl:.3f}  DMED mean={s_dmed:.3f}  diff(DMED-ctrl)={obs_diff:.3f}")

    # ---- exact label permutation over the USED samples (donor-level inference; NOT cell-level) ----
    samps = list(used)
    n_dmed = len(dmed)
    diffs = []
    for combo in itertools.combinations(samps, n_dmed):   # all assignments of n_dmed "DMED" labels
        d = score[list(combo)].mean() - score[[s for s in samps if s not in combo]].mean()
        diffs.append(d)
    diffs = np.array(diffs)
    # one-sided p for CONSERVATION (obs_diff as extreme & positive): fraction of perms >= obs
    p_cons = float(np.mean(diffs >= obs_diff))
    # two-sided
    p_two = float(np.mean(np.abs(diffs) >= abs(obs_diff)))
    n_perm = len(diffs)
    print(f"\nexact label permutation ({len(dmed)}-vs-{len(ctrl)} over {len(used)} used samples, "
          f"{n_perm} assignments):")
    print(f"  one-sided p (signature conserved/UP in DMED): {p_cons:.3f}")
    print(f"  two-sided p: {p_two:.3f}")

    # fraction of core genes that go the "right" (conserved=up) way
    n_up = int(gtab.loc[gtab.rat.isin(core5), "direction_up_in_DMED"].sum())
    print(f"  shared-ED core genes up in DMED: {n_up}/{len(core5)} "
          f"({'inverts/not conserved' if n_up < len(core5)/2 else 'partial'} )")

    # ---- depth confound headline numbers ----
    md_ctrl = depth.loc[depth.condition=="control","median_genes"]
    md_dmed = depth.loc[depth.condition=="DMED","median_genes"]
    print(f"\ndepth confound: median genes/cell  control={md_ctrl.mean():.0f}  "
          f"DMED={md_dmed.mean():.0f}  (ratio {md_ctrl.mean()/md_dmed.mean():.2f}x)")

    # ---- write output table ----
    # one tidy CSV with: depth block + per-gene direction block + verdict block
    depth_out = depth.assign(block="depth_confound")
    gene_out = gtab.assign(block="gene_direction")
    verdict = pd.DataFrame([dict(block="verdict",
        signature="shared_ED_core_5gene",
        samples_used=f"DMED={len(dmed)}({';'.join(dmed)}),control={len(ctrl)}({';'.join(ctrl)})",
        score_control=round(s_ctrl,3), score_DMED=round(s_dmed,3),
        score_diff_DMEDminusCtrl=round(obs_diff,3), perm_p_one_sided_conserved=round(p_cons,3),
        perm_p_two_sided=round(p_two,3), n_perm=n_perm, core_up_in_DMED=f"{n_up}/{len(core5)}",
        sh3bp5_erap2="absent_in_rat_panel_not_testable",
        median_genes_control_all6=round(md_ctrl.mean(),0), median_genes_DMED_all6=round(md_dmed.mean(),0),
        depth_ratio_ctrl_over_dmed=round(md_ctrl.mean()/md_dmed.mean(),2))])
    out_path = f"{OUT}/05_rat_negative.csv"
    with open(out_path, "w") as fh:
        fh.write("# GSE259299 rat STZ-DMED corpus cavernosum pericyte pseudobulk — negative control\n")
        fh.write("# block=depth_confound : per-sample raw/QC depth & pericyte counts\n")
        depth_out.to_csv(fh, index=False)
        fh.write("\n# block=gene_direction : donor-level direction of core orthologs (logCPM, DMED vs control)\n")
        gene_out.to_csv(fh, index=False)
        fh.write("\n# block=verdict : signature score + exact 4-vs-2 donor-level permutation (NOT cell-level)\n")
        verdict.to_csv(fh, index=False)
    # also save the per-sample logCPM of core genes for transparency
    logcpm.loc[core_present].T.to_csv(f"{OUT}/05_rat_pericyte_core_logCPM.csv")
    print(f"\nsaved -> {out_path}")
    print(f"saved -> {OUT}/05_rat_pericyte_core_logCPM.csv")

    # descriptive pericyte marker means (sanity that we annotated pericytes)
    print("\n=== descriptive pericyte marker means (log-norm, per sample) — sanity, NOT inference ===")
    dd = pd.DataFrame(peri_desc).T
    print(dd.round(2).to_string())


if __name__ == "__main__":
    main()
