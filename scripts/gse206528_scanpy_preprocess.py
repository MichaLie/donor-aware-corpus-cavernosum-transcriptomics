#!/usr/bin/env python3
"""Initial Scanpy preprocessing and exploratory clustering for GSE206528.

This step is for visualization and major-cell-type annotation. It does not
perform disease differential expression.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import scanpy as sc
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
IN_FILE = ROOT / "data" / "processed" / "GSE206528" / "GSE206528_raw_counts_combined.h5ad"
OUT_FILE = ROOT / "data" / "processed" / "GSE206528" / "GSE206528_scanpy_preprocessed.h5ad"
FIG_DIR = ROOT / "results" / "figures" / "GSE206528"
TABLE_DIR = ROOT / "results" / "tables"
MARKERS = ROOT / "metadata" / "GSE206528" / "canonical_marker_sets.tsv"


def load_marker_sets() -> dict[str, list[str]]:
    markers = pd.read_csv(MARKERS, sep="\t")
    return {
        row.cell_type: [gene.strip() for gene in row.markers.split(",") if gene.strip()]
        for row in markers.itertuples(index=False)
    }


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = FIG_DIR
    sc.settings.verbosity = 2

    adata = sc.read_h5ad(IN_FILE)
    print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    adata.layers["counts"] = adata.X.copy()

    sc.pp.filter_genes(adata, min_cells=10)
    print(f"After gene filtering: {adata.n_obs} cells x {adata.n_vars} genes")

    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=3000,
        batch_key="donor_id",
    )

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    marker_sets = load_marker_sets()
    for cell_type, genes in marker_sets.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if present:
            sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{cell_type}")

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.uns["pca"] = adata_hvg.uns["pca"]
    adata.varm["PCs"] = pd.DataFrame(index=adata.var_names)

    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, min_dist=0.3)

    for resolution in [0.3, 0.5, 0.8, 1.2]:
        key = f"leiden_r{str(resolution).replace('.', '_')}"
        sc.tl.leiden(adata, resolution=resolution, key_added=key, flavor="igraph")

    cluster_key = "leiden_r0_5"
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        n_genes=100,
        pts=True,
    )
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(TABLE_DIR / "GSE206528_leiden_r0_5_marker_genes.tsv", sep="\t", index=False)

    composition = (
        adata.obs.groupby(["donor_id", "primary_group", "ed_subgroup", cluster_key], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    composition.to_csv(TABLE_DIR / "GSE206528_leiden_r0_5_cluster_counts_by_donor.tsv", sep="\t", index=False)

    score_cols = [col for col in adata.obs.columns if col.startswith("score_")]
    cluster_scores = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    cluster_scores.to_csv(TABLE_DIR / "GSE206528_leiden_r0_5_marker_scores_by_cluster.tsv", sep="\t")

    sc.pl.umap(
        adata,
        color=["donor_id", "ed_subgroup", cluster_key],
        wspace=0.35,
        save="_GSE206528_donor_group_cluster.png",
        show=False,
    )
    marker_panel = [
        "COL1A1",
        "ACTA2",
        "PECAM1",
        "RGS5",
        "C1QA",
        "CD3D",
        "TPSAB1",
        "S100B",
        "EPCAM",
    ]
    marker_panel = [gene for gene in marker_panel if gene in adata.var_names]
    sc.pl.dotplot(
        adata,
        var_names=marker_panel,
        groupby=cluster_key,
        standard_scale="var",
        save="_GSE206528_leiden_r0_5_canonical_markers.png",
        show=False,
    )

    adata.write_h5ad(OUT_FILE, compression="gzip")
    print(f"Wrote {OUT_FILE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
