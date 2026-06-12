#!/usr/bin/env python3
"""Build a combined GSE206528 AnnData object from per-donor h5ad files."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "processed" / "GSE206528" / "per_sample_h5ad_manifest.tsv"
CELL_QC = ROOT / "results" / "tables" / "GSE206528_cell_qc_metrics.tsv.gz"
OUT_DIR = ROOT / "data" / "processed" / "GSE206528"
OUT_FILE = OUT_DIR / "GSE206528_raw_counts_combined.h5ad"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(MANIFEST, sep="\t")
    qc = pd.read_csv(CELL_QC, sep="\t")
    qc["qc_key"] = qc["geo_accession"].astype(str) + "|" + qc["cell_id"].astype(str)
    qc = qc.set_index("qc_key")

    adatas = []
    for row in manifest.itertuples(index=False):
        path = ROOT / row.h5ad_file
        sample = ad.read_h5ad(path)

        sample.obs["original_cell_id"] = sample.obs_names.astype(str)
        sample.obs_names = (
            sample.obs["geo_accession"].astype(str)
            + "_"
            + sample.obs["original_cell_id"].astype(str).str.replace(" ", "_", regex=False)
        )
        sample.obs_names_make_unique()

        sample.obs["qc_key"] = (
            sample.obs["geo_accession"].astype(str)
            + "|"
            + sample.obs["original_cell_id"].astype(str)
        )
        qc_cols = ["total_counts", "n_genes_by_counts", "mt_counts", "pct_mt"]
        sample.obs = sample.obs.join(qc[qc_cols], on="qc_key")
        sample.obs = sample.obs.drop(columns=["qc_key"])

        adatas.append(sample)
        print(f"Loaded {row.sample_title}: {sample.n_obs} cells x {sample.n_vars} genes")

    combined = ad.concat(
        adatas,
        axis=0,
        join="inner",
        merge="same",
        uns_merge="unique",
        label="sample_title_from_concat",
        keys=[a.obs["sample_title"].iloc[0] for a in adatas],
        index_unique=None,
    )
    combined.obs = combined.obs.drop(columns=["sample_title_from_concat"])
    combined.obs_names_make_unique()
    combined.var["gene_symbol"] = combined.var_names.astype(str)
    combined.uns["dataset_id"] = "GSE206528"
    combined.uns["description"] = "Combined raw-count AnnData object from donor-level GEO matrices."

    combined.write_h5ad(OUT_FILE, compression="gzip")
    print(f"Wrote {OUT_FILE.relative_to(ROOT)}")
    print(f"Combined shape: {combined.n_obs} cells x {combined.n_vars} genes")


if __name__ == "__main__":
    main()
