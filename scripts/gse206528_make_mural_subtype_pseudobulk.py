#!/usr/bin/env python3
"""Create GSE206528 donor pseudobulk profiles for refined mural subtypes."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
RAW_H5AD = ROOT / "data/processed/GSE206528/GSE206528_raw_counts_combined.h5ad"
PREPROCESSED_H5AD = ROOT / "data/processed/GSE206528/GSE206528_scanpy_preprocessed.h5ad"
OUT_DIR = ROOT / "data/processed/GSE206528/pseudobulk_mural_subtypes"
TABLE_DIR = ROOT / "results/tables"
META_DIR = ROOT / "metadata/GSE206528"


MURAL_SUBTYPE_MAP = {
    "6": "Pericyte_like",
    "8": "Vessel_SMC_like",
    "12": "Vessel_SMC_like",
    "13": "CC_SMC_like",
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    META_DIR.mkdir(parents=True, exist_ok=True)

    subtype_table = pd.DataFrame(
        [
            {
                "leiden_r0_5": cluster,
                "mural_subtype": subtype,
                "annotation_basis": "author_marker_crosscheck",
            }
            for cluster, subtype in MURAL_SUBTYPE_MAP.items()
        ]
    ).sort_values("leiden_r0_5")
    subtype_table.to_csv(META_DIR / "mural_subtype_leiden_r0_5_annotations.tsv", sep="\t", index=False)

    pre = ad.read_h5ad(PREPROCESSED_H5AD, backed="r")
    obs_annotation = pre.obs[["leiden_r0_5"]].copy()
    pre.file.close()
    obs_annotation["mural_subtype"] = obs_annotation["leiden_r0_5"].astype(str).map(MURAL_SUBTYPE_MAP)
    obs_annotation = obs_annotation.dropna(subset=["mural_subtype"])

    raw = ad.read_h5ad(RAW_H5AD)
    raw = raw[obs_annotation.index].copy()
    raw.obs = raw.obs.join(obs_annotation[["leiden_r0_5", "mural_subtype"]], how="left")
    if raw.obs["mural_subtype"].isna().any():
        raise ValueError("Some selected mural cells did not receive a refined subtype")

    group_cols = ["donor_id", "primary_group", "ed_subgroup", "diabetes_status", "mural_subtype"]
    pb_meta = raw.obs[group_cols].drop_duplicates().sort_values(group_cols).reset_index(drop=True)
    pb_meta["pseudobulk_id"] = (
        pb_meta["donor_id"].astype(str) + "__" + pb_meta["mural_subtype"].astype(str)
    )

    group_lookup = {tuple(row[col] for col in group_cols): i for i, row in pb_meta.iterrows()}
    group_index = raw.obs[group_cols].apply(lambda row: group_lookup[tuple(row)], axis=1).to_numpy()

    membership = sparse.csr_matrix(
        (
            np.ones(raw.n_obs, dtype=np.int8),
            (group_index, np.arange(raw.n_obs)),
        ),
        shape=(len(pb_meta), raw.n_obs),
    )
    pb = (membership @ raw.X).tocsr()

    cell_counts = (
        raw.obs.groupby(["donor_id", "mural_subtype"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    pb_meta = pb_meta.merge(cell_counts, on=["donor_id", "mural_subtype"], how="left")

    count_df = pd.DataFrame.sparse.from_spmatrix(
        pb.T,
        index=raw.var_names.astype(str),
        columns=pb_meta["pseudobulk_id"],
    )
    count_df.index.name = "gene_symbol"

    count_df.to_csv(OUT_DIR / "GSE206528_mural_subtype_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")
    pb_meta.to_csv(OUT_DIR / "GSE206528_mural_subtype_pseudobulk_metadata.tsv", sep="\t", index=False)
    cell_counts.to_csv(TABLE_DIR / "GSE206528_mural_subtype_counts_by_donor.tsv", sep="\t", index=False)
    raw.obs[["donor_id", "ed_subgroup", "leiden_r0_5", "mural_subtype"]].to_csv(
        TABLE_DIR / "GSE206528_cell_mural_subtype_annotations.tsv.gz",
        sep="\t",
        index=True,
    )

    print(f"Refined mural cells: {raw.n_obs}")
    print(f"Pseudobulk matrix: {count_df.shape[0]} genes x {count_df.shape[1]} donor-subtype profiles")


if __name__ == "__main__":
    main()
