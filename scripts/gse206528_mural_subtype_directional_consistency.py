#!/usr/bin/env python3
"""Donor-direction checks for refined GSE206528 mural subtype edgeR results."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PSEUDOBULK_DIR = ROOT / "data/processed/GSE206528/pseudobulk_mural_subtypes"
RESULT_DIR = ROOT / "results/tables/GSE206528_mural_subtype_edger"

COUNTS_FILE = PSEUDOBULK_DIR / "GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
META_FILE = PSEUDOBULK_DIR / "GSE206528_mural_subtype_pseudobulk_metadata.tsv"
RESULTS_FILE = RESULT_DIR / "all_mural_subtype_edgeR_results.tsv.gz"

COMPARISON_GROUPS = {
    "organic_ED_vs_normal": ("primary_group", "normal", "organic_ED"),
    "non_diabetic_ED_vs_normal": ("ed_subgroup", "normal", "non_diabetic_ED"),
    "diabetic_ED_vs_normal": ("ed_subgroup", "normal", "diabetic_ED"),
}


def main() -> None:
    counts = pd.read_csv(COUNTS_FILE, sep="\t", index_col=0)
    meta = pd.read_csv(META_FILE, sep="\t").set_index("pseudobulk_id")
    results = pd.read_csv(RESULTS_FILE, sep="\t")

    counts = counts.loc[:, meta.index]
    library_sizes = counts.sum(axis=0)
    log_cpm = np.log2((counts + 0.5).div(library_sizes + 1, axis=1) * 1_000_000)

    rows: list[dict[str, object]] = []
    candidates = results.loc[results["FDR"] < 0.10].copy()

    for row in candidates.itertuples(index=False):
        group_column, control_label, disease_label = COMPARISON_GROUPS[row.comparison]
        sample_meta = meta.loc[
            (meta["mural_subtype"] == row.mural_subtype)
            & (meta["n_cells"] >= 50)
            & (meta[group_column].isin([control_label, disease_label]))
        ].copy()

        if row.gene_symbol not in log_cpm.index or sample_meta.empty:
            continue

        values = log_cpm.loc[row.gene_symbol, sample_meta.index]
        control_values = values.loc[sample_meta[group_column] == control_label]
        disease_values = values.loc[sample_meta[group_column] == disease_label]
        if control_values.empty or disease_values.empty:
            continue

        direction = "up_in_ED" if row.logFC > 0 else "down_in_ED"
        control_median = float(control_values.median())

        if row.logFC > 0:
            consistent = disease_values > control_median
            beyond_control_range = disease_values > control_values.max()
        else:
            consistent = disease_values < control_median
            beyond_control_range = disease_values < control_values.min()

        rows.append(
            {
                "mural_subtype": row.mural_subtype,
                "comparison": row.comparison,
                "gene_symbol": row.gene_symbol,
                "direction": direction,
                "logFC": row.logFC,
                "logCPM": row.logCPM,
                "PValue": row.PValue,
                "FDR": row.FDR,
                "n_control_donors": int(control_values.shape[0]),
                "n_disease_donors": int(disease_values.shape[0]),
                "control_mean_logCPM": float(control_values.mean()),
                "disease_mean_logCPM": float(disease_values.mean()),
                "disease_minus_control_mean_logCPM": float(disease_values.mean() - control_values.mean()),
                "disease_donors_consistent_vs_control_median": int(consistent.sum()),
                "disease_fraction_consistent_vs_control_median": float(consistent.mean()),
                "disease_donors_beyond_control_range": int(beyond_control_range.sum()),
                "control_donor_values_logCPM": ";".join(
                    f"{sample}:{value:.3f}" for sample, value in control_values.items()
                ),
                "disease_donor_values_logCPM": ";".join(
                    f"{sample}:{value:.3f}" for sample, value in disease_values.items()
                ),
            }
        )

    consistency = pd.DataFrame(rows).sort_values(
        ["mural_subtype", "comparison", "FDR", "PValue", "gene_symbol"]
    )
    consistency.to_csv(
        RESULT_DIR / "FDR_0_10_genes_with_directional_consistency.tsv",
        sep="\t",
        index=False,
    )

    robust = consistency.loc[
        (consistency["FDR"] < 0.05)
        & (consistency["disease_fraction_consistent_vs_control_median"] >= 0.8)
    ].copy()
    robust.to_csv(
        RESULT_DIR / "FDR_0_05_directionally_consistent_candidates.tsv",
        sep="\t",
        index=False,
    )

    overview = (
        consistency.assign(
            is_fdr_0_05=lambda df: df["FDR"] < 0.05,
            is_directionally_consistent=lambda df: df["disease_fraction_consistent_vs_control_median"] >= 0.8,
            is_fdr_0_05_directionally_consistent=lambda df: df["is_fdr_0_05"]
            & df["is_directionally_consistent"],
        )
        .groupby(["mural_subtype", "comparison"], as_index=False)
        .agg(
            n_FDR_0_10=("gene_symbol", "size"),
            n_FDR_0_05=("is_fdr_0_05", "sum"),
            n_FDR_0_05_directionally_consistent=("is_fdr_0_05_directionally_consistent", "sum"),
        )
    )
    overview.to_csv(
        RESULT_DIR / "directional_consistency_overview.tsv",
        sep="\t",
        index=False,
    )
    print("Wrote refined mural directional consistency summaries to:", RESULT_DIR)


if __name__ == "__main__":
    main()
