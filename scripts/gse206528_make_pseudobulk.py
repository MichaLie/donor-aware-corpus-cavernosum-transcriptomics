#!/usr/bin/env python3
"""Create donor-by-cell-type pseudobulk count matrices for GSE206528."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
RAW_H5AD = ROOT / "data" / "processed" / "GSE206528" / "GSE206528_raw_counts_combined.h5ad"
PREPROCESSED_H5AD = ROOT / "data" / "processed" / "GSE206528" / "GSE206528_scanpy_preprocessed.h5ad"
ANNOTATION = ROOT / "metadata" / "GSE206528" / "major_celltype_leiden_r0_5_annotations.tsv"
OUT_DIR = ROOT / "data" / "processed" / "GSE206528" / "pseudobulk"
TABLE_DIR = ROOT / "results" / "tables"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    annotations = pd.read_csv(ANNOTATION, sep="\t")
    annotations["leiden_r0_5"] = annotations["leiden_r0_5"].astype(str)
    cluster_to_type = annotations.set_index("leiden_r0_5")["major_cell_type"].to_dict()

    pre = ad.read_h5ad(PREPROCESSED_H5AD, backed="r")
    obs_annotation = pre.obs[["leiden_r0_5"]].copy()
    pre.file.close()
    obs_annotation["major_cell_type"] = (
        obs_annotation["leiden_r0_5"].astype(str).map(cluster_to_type)
    )
    if obs_annotation["major_cell_type"].isna().any():
        missing = obs_annotation.loc[
            obs_annotation["major_cell_type"].isna(), "leiden_r0_5"
        ].unique()
        raise ValueError(f"Missing annotations for clusters: {missing}")

    raw = ad.read_h5ad(RAW_H5AD)
    raw.obs = raw.obs.join(obs_annotation, how="left")
    if raw.obs["major_cell_type"].isna().any():
        raise ValueError("Some raw cells did not receive major cell-type annotation")

    group_cols = ["donor_id", "primary_group", "ed_subgroup", "diabetes_status", "major_cell_type"]
    pb_meta = raw.obs[group_cols].drop_duplicates().sort_values(group_cols).reset_index(drop=True)
    pb_meta["pseudobulk_id"] = (
        pb_meta["donor_id"].astype(str)
        + "__"
        + pb_meta["major_cell_type"].astype(str)
    )

    group_lookup = {
        tuple(row[col] for col in group_cols): i for i, row in pb_meta.iterrows()
    }
    group_index = raw.obs[group_cols].apply(lambda row: group_lookup[tuple(row)], axis=1).to_numpy()

    membership = sparse.csr_matrix(
        (
            np.ones(raw.n_obs, dtype=np.int8),
            (group_index, np.arange(raw.n_obs)),
        ),
        shape=(len(pb_meta), raw.n_obs),
    )
    pb = membership @ raw.X
    pb = pb.tocsr()

    cell_counts = (
        raw.obs.groupby(["donor_id", "major_cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    pb_meta = pb_meta.merge(cell_counts, on=["donor_id", "major_cell_type"], how="left")

    gene_symbols = raw.var_names.astype(str)
    count_df = pd.DataFrame.sparse.from_spmatrix(
        pb.T,
        index=gene_symbols,
        columns=pb_meta["pseudobulk_id"],
    )
    count_df.index.name = "gene_symbol"

    counts_out = OUT_DIR / "GSE206528_pseudobulk_counts_by_donor_celltype.tsv.gz"
    meta_out = OUT_DIR / "GSE206528_pseudobulk_sample_metadata.tsv"
    obs_out = TABLE_DIR / "GSE206528_major_celltype_counts_by_donor.tsv"

    count_df.to_csv(counts_out, sep="\t", compression="gzip")
    pb_meta.to_csv(meta_out, sep="\t", index=False)
    cell_counts.to_csv(obs_out, sep="\t", index=False)

    raw.obs[["donor_id", "ed_subgroup", "leiden_r0_5", "major_cell_type"]].to_csv(
        TABLE_DIR / "GSE206528_cell_major_celltype_annotations.tsv.gz",
        sep="\t",
        index=True,
    )

    print(f"Wrote {counts_out.relative_to(ROOT)}")
    print(f"Wrote {meta_out.relative_to(ROOT)}")
    print(f"Wrote {obs_out.relative_to(ROOT)}")
    print(f"Pseudobulk matrix: {count_df.shape[0]} genes x {count_df.shape[1]} donor-cell-type profiles")


if __name__ == "__main__":
    main()
