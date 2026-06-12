#!/usr/bin/env python3
"""Stream QC metrics for GSE206528 dense GEO count matrices.

The GEO supplementary files are gzipped CSV matrices with genes in rows and
single cells in columns. Loading all donors as one dense matrix is unnecessary
and memory-heavy, so this script computes basic per-cell QC in chunks.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw" / "GSE206528"
SAMPLE_TABLE = ROOT / "metadata" / "GSE206528" / "samples.tsv"
OUT_DIR = ROOT / "results" / "tables"
CELL_QC_OUT = OUT_DIR / "GSE206528_cell_qc_metrics.tsv.gz"
SAMPLE_QC_OUT = OUT_DIR / "GSE206528_sample_qc_summary.tsv"


def _normalise_gene_index(index: pd.Index) -> pd.Series:
    return pd.Series(index.astype(str), index=index).str.upper()


def stream_sample_qc(sample: pd.Series, chunksize: int = 250) -> tuple[pd.DataFrame, dict]:
    matrix_path = RAW_DIR / sample["matrix_file"]
    if not matrix_path.exists():
        raise FileNotFoundError(matrix_path)

    total_counts = None
    n_genes_by_counts = None
    mt_counts = None
    cell_ids = None
    n_gene_rows = 0

    for chunk in pd.read_csv(matrix_path, index_col=0, chunksize=chunksize):
        if cell_ids is None:
            cell_ids = chunk.columns.astype(str).to_numpy()
            n_cells = len(cell_ids)
            total_counts = np.zeros(n_cells, dtype=np.int64)
            n_genes_by_counts = np.zeros(n_cells, dtype=np.int32)
            mt_counts = np.zeros(n_cells, dtype=np.int64)

        values = chunk.to_numpy(dtype=np.int32, copy=False)
        total_counts += values.sum(axis=0, dtype=np.int64)
        n_genes_by_counts += (values > 0).sum(axis=0, dtype=np.int32)

        gene_index = _normalise_gene_index(chunk.index)
        mt_mask = gene_index.str.startswith("MT-").to_numpy()
        if mt_mask.any():
            mt_counts += values[mt_mask, :].sum(axis=0, dtype=np.int64)

        n_gene_rows += chunk.shape[0]

    if cell_ids is None:
        raise ValueError(f"No cells found in {matrix_path}")

    pct_mt = np.divide(
        mt_counts,
        total_counts,
        out=np.zeros_like(mt_counts, dtype=float),
        where=total_counts > 0,
    ) * 100.0

    cell_metrics = pd.DataFrame(
        {
            "geo_accession": sample["geo_accession"],
            "sample_title": sample["sample_title"],
            "donor_id": sample["donor_id"],
            "primary_group": sample["primary_group"],
            "ed_subgroup": sample["ed_subgroup"],
            "diabetes_status": sample["diabetes_status"],
            "cell_id": cell_ids,
            "total_counts": total_counts,
            "n_genes_by_counts": n_genes_by_counts,
            "mt_counts": mt_counts,
            "pct_mt": pct_mt,
        }
    )

    summary = {
        "geo_accession": sample["geo_accession"],
        "sample_title": sample["sample_title"],
        "donor_id": sample["donor_id"],
        "primary_group": sample["primary_group"],
        "ed_subgroup": sample["ed_subgroup"],
        "diabetes_status": sample["diabetes_status"],
        "matrix_file": sample["matrix_file"],
        "n_gene_rows": n_gene_rows,
        "n_cells": len(cell_ids),
        "median_total_counts": float(np.median(total_counts)),
        "median_n_genes_by_counts": float(np.median(n_genes_by_counts)),
        "median_pct_mt": float(np.median(pct_mt)),
        "cells_pct_mt_gt_20": int((pct_mt > 20).sum()),
        "cells_total_counts_lt_500": int((total_counts < 500).sum()),
        "cells_n_genes_lt_200": int((n_genes_by_counts < 200).sum()),
    }

    return cell_metrics, summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(SAMPLE_TABLE, sep="\t")

    all_cell_metrics = []
    summaries = []
    for _, sample in samples.iterrows():
        cell_metrics, summary = stream_sample_qc(sample)
        all_cell_metrics.append(cell_metrics)
        summaries.append(summary)
        print(
            f"{sample['sample_title']}: "
            f"{summary['n_cells']} cells, "
            f"median genes {summary['median_n_genes_by_counts']:.0f}, "
            f"median UMIs {summary['median_total_counts']:.0f}, "
            f"median mt {summary['median_pct_mt']:.2f}%"
        )

    pd.concat(all_cell_metrics, ignore_index=True).to_csv(
        CELL_QC_OUT, sep="\t", index=False
    )
    pd.DataFrame(summaries).to_csv(SAMPLE_QC_OUT, sep="\t", index=False)
    print(f"Wrote {CELL_QC_OUT}")
    print(f"Wrote {SAMPLE_QC_OUT}")


if __name__ == "__main__":
    main()
