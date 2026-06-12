#!/usr/bin/env python3
"""Convert GSE206528 GEO dense CSV matrices to per-donor sparse h5ad files."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw" / "GSE206528"
SAMPLE_TABLE = ROOT / "metadata" / "GSE206528" / "samples.tsv"
OUT_DIR = ROOT / "data" / "processed" / "GSE206528" / "per_sample_h5ad"
MANIFEST_OUT = ROOT / "data" / "processed" / "GSE206528" / "per_sample_h5ad_manifest.tsv"


def convert_sample(sample: pd.Series, chunksize: int = 250) -> dict:
    matrix_path = RAW_DIR / sample["matrix_file"]
    out_path = OUT_DIR / f"{sample['geo_accession']}_{sample['sample_title']}.h5ad"

    if not matrix_path.exists():
        raise FileNotFoundError(matrix_path)

    x_chunks = []
    gene_names = []
    cell_ids = None

    for chunk in pd.read_csv(matrix_path, index_col=0, chunksize=chunksize):
        if cell_ids is None:
            cell_ids = chunk.columns.astype(str).to_list()
        gene_names.extend(chunk.index.astype(str).to_list())
        x_chunks.append(sparse.csr_matrix(chunk.to_numpy(dtype="int32").T))

    if cell_ids is None:
        raise ValueError(f"No cells found in {matrix_path}")

    x = sparse.hstack(x_chunks, format="csr", dtype="int32")

    obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
    for column in [
        "geo_accession",
        "sample_title",
        "donor_id",
        "primary_group",
        "ed_subgroup",
        "diabetes_status",
        "tissue",
        "organism",
    ]:
        obs[column] = sample[column]

    var = pd.DataFrame(index=pd.Index(gene_names, name="gene_symbol"))
    var["gene_symbol"] = var.index

    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.uns["dataset_id"] = "GSE206528"
    adata.uns["source_matrix_file"] = sample["matrix_file"]
    adata.write_h5ad(out_path, compression="gzip")

    return {
        "geo_accession": sample["geo_accession"],
        "sample_title": sample["sample_title"],
        "donor_id": sample["donor_id"],
        "primary_group": sample["primary_group"],
        "ed_subgroup": sample["ed_subgroup"],
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "nnz": int(adata.X.nnz),
        "h5ad_file": str(out_path.relative_to(ROOT)),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(SAMPLE_TABLE, sep="\t")

    manifest_rows = []
    for _, sample in samples.iterrows():
        row = convert_sample(sample)
        manifest_rows.append(row)
        print(
            f"{row['sample_title']}: wrote {row['h5ad_file']} "
            f"({row['n_cells']} cells, {row['n_genes']} genes, {row['nnz']} nonzero counts)"
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(MANIFEST_OUT, sep="\t", index=False)
    print(f"Wrote {MANIFEST_OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
