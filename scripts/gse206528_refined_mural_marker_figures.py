#!/usr/bin/env python3
"""Create donor-level and marker-support figure outputs.

Outputs:
- donor-level logCPM table and plots for the 12 diabetic ED pericyte-like genes;
- refined mural subtype marker-support table and plot for the pericyte-like/perivascular label.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]

PSEUDOBULK_DIR = ROOT / "data/processed/GSE206528/pseudobulk_mural_subtypes"
PSEUDOBULK_COUNTS = PSEUDOBULK_DIR / "GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
PSEUDOBULK_META = PSEUDOBULK_DIR / "GSE206528_mural_subtype_pseudobulk_metadata.tsv"
EDGER_RESULTS = ROOT / "results/tables/GSE206528_mural_subtype_edger/all_mural_subtype_edgeR_results.tsv.gz"
PREPROCESSED_H5AD = ROOT / "data/processed/GSE206528/GSE206528_scanpy_preprocessed.h5ad"
AUTHOR_CROSSWALK = ROOT / "results/tables/GSE206528_author_marker_crosswalk_top3.tsv"

FIG_DIR = ROOT / "results/figures/GSE206528"
TABLE_DIR = ROOT / "results/tables/GSE206528_refined_mural_marker_support"

MURAL_SUBTYPE_MAP = {
    "6": "Pericyte_like",
    "8": "Vessel_SMC_like",
    "12": "Vessel_SMC_like",
    "13": "CC_SMC_like",
}

CANDIDATE_GENES = [
    "ADAMTS4",
    "NR4A3",
    "IL6",
    "PI15",
    "PLA2G2A",
    "SH3BP5",
    "MT1A",
    "NR4A2",
    "ERAP2",
    "NR4A1",
    "CYR61",
    "RGS16",
]

MAIN_TEXT_REPRESENTATIVE_GENES = ["IL6", "ADAMTS4", "PLA2G2A", "CYR61", "PI15"]

MARKER_PANEL = [
    "PDGFRB",
    "RGS5",
    "CSPG4",
    "MCAM",
    "NOTCH3",
    "ABCC9",
    "ACTA2",
    "MYH11",
    "TAGLN",
    "CNN1",
    "DES",
    "MYLK",
    "PDGFRA",
    "DCN",
    "COL1A1",
    "PECAM1",
    "VWF",
    "PTPRC",
]

SUBTYPE_ORDER = ["Pericyte_like", "Vessel_SMC_like", "CC_SMC_like"]
SUBTYPE_LABELS = {
    "Pericyte_like": "Pericyte-like\nperivascular",
    "Vessel_SMC_like": "Vessel SMC-like",
    "CC_SMC_like": "CC SMC-like",
}
GROUP_LABELS = {
    "normal": "Normal",
    "diabetic_ED": "Diabetic ED",
}
GROUP_COLORS = {
    "normal": "#0072B2",
    "diabetic_ED": "#D55E00",
}


def save_figure(fig: plt.Figure, stem: str) -> None:
    for ext in ["png", "pdf"]:
        fig.savefig(FIG_DIR / f"{stem}.{ext}", dpi=450 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def load_pseudobulk_logcpm() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(PSEUDOBULK_COUNTS, sep="\t", index_col=0)
    meta = pd.read_csv(PSEUDOBULK_META, sep="\t").set_index("pseudobulk_id")
    counts = counts.loc[:, meta.index]
    library_sizes = counts.sum(axis=0)
    log_cpm = np.log2((counts + 0.5).div(library_sizes + 1, axis=1) * 1_000_000)
    return log_cpm, meta


def build_donor_candidate_table() -> pd.DataFrame:
    log_cpm, meta = load_pseudobulk_logcpm()
    edger = pd.read_csv(EDGER_RESULTS, sep="\t")
    stats = edger.loc[
        (edger["mural_subtype"] == "Pericyte_like")
        & (edger["comparison"] == "diabetic_ED_vs_normal")
        & (edger["gene_symbol"].isin(CANDIDATE_GENES)),
        ["gene_symbol", "logFC", "logCPM", "PValue", "FDR"],
    ].copy()

    selected_meta = meta.loc[
        (meta["mural_subtype"] == "Pericyte_like")
        & (meta["n_cells"] >= 50)
        & (meta["ed_subgroup"].isin(["normal", "diabetic_ED"]))
    ].copy()
    selected_meta["display_group"] = selected_meta["ed_subgroup"].map(GROUP_LABELS)

    missing = [gene for gene in CANDIDATE_GENES if gene not in log_cpm.index]
    if missing:
        raise ValueError(f"Candidate genes missing from pseudobulk matrix: {missing}")

    rows: list[dict[str, object]] = []
    for gene in CANDIDATE_GENES:
        gene_stats = stats.loc[stats["gene_symbol"] == gene]
        if gene_stats.empty:
            raise ValueError(f"Candidate gene missing from edgeR result table: {gene}")
        stat_row = gene_stats.iloc[0].to_dict()
        for pseudobulk_id, meta_row in selected_meta.iterrows():
            rows.append(
                {
                    "gene_symbol": gene,
                    "pseudobulk_id": pseudobulk_id,
                    "donor_id": meta_row["donor_id"],
                    "ed_subgroup": meta_row["ed_subgroup"],
                    "display_group": meta_row["display_group"],
                    "n_cells": int(meta_row["n_cells"]),
                    "log2_cpm": float(log_cpm.loc[gene, pseudobulk_id]),
                    "logFC_diabetic_ED_vs_normal": float(stat_row["logFC"]),
                    "edgeR_logCPM": float(stat_row["logCPM"]),
                    "PValue": float(stat_row["PValue"]),
                    "FDR": float(stat_row["FDR"]),
                }
            )

    table = pd.DataFrame(rows)
    table["gene_symbol"] = pd.Categorical(table["gene_symbol"], CANDIDATE_GENES, ordered=True)
    table["ed_subgroup"] = pd.Categorical(table["ed_subgroup"], ["normal", "diabetic_ED"], ordered=True)
    table = table.sort_values(["gene_symbol", "ed_subgroup", "donor_id"])
    table.to_csv(
        TABLE_DIR / "diabetic_ed_pericyte_like_12_gene_donor_logcpm.tsv",
        sep="\t",
        index=False,
    )
    return table


def plot_candidate_donor_expression(table: pd.DataFrame, genes: list[str], stem: str, ncols: int) -> None:
    gene_stats = (
        table[["gene_symbol", "logFC_diabetic_ED_vs_normal", "FDR"]]
        .drop_duplicates()
        .set_index("gene_symbol")
    )
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.15, nrows * 2.1), sharex=True)
    axes_flat = np.asarray(axes).reshape(-1)
    x_positions = {"normal": 0, "diabetic_ED": 1}
    donor_offsets = {
        "N1": -0.12,
        "N2": 0.0,
        "N3": 0.12,
        "DMED1": -0.08,
        "DMED2": 0.08,
    }

    for ax, gene in zip(axes_flat, genes):
        sub = table.loc[table["gene_symbol"].astype(str) == gene].copy()
        for _, row in sub.iterrows():
            x = x_positions[row["ed_subgroup"]] + donor_offsets.get(row["donor_id"], 0.0)
            ax.scatter(
                x,
                row["log2_cpm"],
                s=34,
                color=GROUP_COLORS[row["ed_subgroup"]],
                edgecolor="black",
                linewidth=0.35,
                zorder=3,
            )
            ax.text(
                x,
                row["log2_cpm"] + 0.08,
                row["donor_id"],
                ha="center",
                va="bottom",
                fontsize=5.8,
            )

        for subgroup, xpos in x_positions.items():
            values = sub.loc[sub["ed_subgroup"] == subgroup, "log2_cpm"].to_numpy()
            if len(values):
                ax.hlines(
                    np.median(values),
                    xpos - 0.24,
                    xpos + 0.24,
                    color="black",
                    linewidth=1.0,
                    zorder=2,
                )

        normal_median = sub.loc[sub["ed_subgroup"] == "normal", "log2_cpm"].median()
        ax.axhline(normal_median, color="#999999", linestyle=":", linewidth=0.8, zorder=1)
        stat = gene_stats.loc[gene]
        direction = "up" if stat["logFC_diabetic_ED_vs_normal"] > 0 else "down"
        ax.set_title(
            f"{gene} ({direction}, FDR={stat['FDR']:.3g})",
            fontsize=8,
            pad=5,
        )
        ax.set_xlim(-0.45, 1.45)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Normal", "Diabetic\nED"], fontsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.5, alpha=0.8)

    for ax in axes_flat[len(genes) :]:
        ax.axis("off")

    fig.supylabel("Donor pseudobulk expression (log2 CPM)", x=0.025, fontsize=9)
    fig.suptitle(
        "Diabetic ED pericyte-like candidates shown at donor level",
        y=1.005,
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout(w_pad=1.0, h_pad=1.2)
    save_figure(fig, stem)


def extract_dense_matrix(adata: ad.AnnData, gene_indices: list[int]) -> np.ndarray:
    x = adata[:, gene_indices].X
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def extract_dense_layer(adata: ad.AnnData, layer: str, gene_indices: list[int]) -> np.ndarray:
    x = adata[:, gene_indices].layers[layer]
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def build_marker_support_table() -> pd.DataFrame:
    adata = ad.read_h5ad(PREPROCESSED_H5AD)
    adata.obs["mural_subtype"] = adata.obs["leiden_r0_5"].astype(str).map(MURAL_SUBTYPE_MAP)
    adata = adata[adata.obs["mural_subtype"].notna()].copy()

    present_markers = [gene for gene in MARKER_PANEL if gene in adata.var_names]
    missing_markers = sorted(set(MARKER_PANEL) - set(present_markers))
    if missing_markers:
        pd.Series(missing_markers, name="missing_marker_gene").to_csv(
            TABLE_DIR / "refined_mural_marker_panel_missing_genes.tsv",
            sep="\t",
            index=False,
        )

    gene_indices = [adata.var_names.get_loc(gene) for gene in present_markers]
    expr = extract_dense_matrix(adata, gene_indices)
    counts = extract_dense_layer(adata, "counts", gene_indices)

    rows: list[dict[str, object]] = []
    for subtype in SUBTYPE_ORDER:
        mask = (adata.obs["mural_subtype"].to_numpy() == subtype)
        subtype_expr = expr[mask]
        subtype_counts = counts[mask]
        for j, gene in enumerate(present_markers):
            rows.append(
                {
                    "mural_subtype": subtype,
                    "gene_symbol": gene,
                    "mean_log1p_normalized_expression": float(np.mean(subtype_expr[:, j])),
                    "pct_cells_detected": float(np.mean(subtype_counts[:, j] > 0) * 100.0),
                    "n_cells": int(mask.sum()),
                }
            )

    marker_table = pd.DataFrame(rows)
    marker_table["mural_subtype"] = pd.Categorical(marker_table["mural_subtype"], SUBTYPE_ORDER, ordered=True)
    marker_table["gene_symbol"] = pd.Categorical(marker_table["gene_symbol"], present_markers, ordered=True)
    marker_table = marker_table.sort_values(["mural_subtype", "gene_symbol"])
    marker_table.to_csv(
        TABLE_DIR / "refined_mural_subtype_marker_support.tsv",
        sep="\t",
        index=False,
    )

    author = pd.read_csv(AUTHOR_CROSSWALK, sep="\t")
    author_support = author.loc[
        (author["comparison_level"] == "author_11_subclusters")
        & (author["leiden_r0_5"].astype(str).isin(MURAL_SUBTYPE_MAP.keys()))
    ].copy()
    author_support["mural_subtype"] = author_support["leiden_r0_5"].astype(str).map(MURAL_SUBTYPE_MAP)
    author_support = author_support.sort_values(
        ["leiden_r0_5", "n_overlap_top100"],
        ascending=[True, False],
    )
    author_support.to_csv(
        TABLE_DIR / "refined_mural_author_marker_crosscheck_top3.tsv",
        sep="\t",
        index=False,
    )
    return marker_table


def plot_marker_support(marker_table: pd.DataFrame) -> None:
    plot_table = marker_table.copy()
    genes = [gene for gene in MARKER_PANEL if gene in plot_table["gene_symbol"].astype(str).unique()]
    plot_table["gene_symbol"] = pd.Categorical(plot_table["gene_symbol"].astype(str), genes, ordered=True)
    plot_table["mural_subtype"] = pd.Categorical(
        plot_table["mural_subtype"].astype(str),
        SUBTYPE_ORDER,
        ordered=True,
    )

    size_min, size_max = 22, 265
    expr_values = plot_table["mean_log1p_normalized_expression"].to_numpy()
    expr_min, expr_max = float(np.nanmin(expr_values)), float(np.nanmax(expr_values))

    fig, ax = plt.subplots(figsize=(8.6, 2.9))
    for _, row in plot_table.iterrows():
        x = genes.index(str(row["gene_symbol"]))
        y = SUBTYPE_ORDER.index(str(row["mural_subtype"]))
        size = size_min + (row["pct_cells_detected"] / 100.0) * (size_max - size_min)
        color_value = 0.5 if expr_max == expr_min else (
            (row["mean_log1p_normalized_expression"] - expr_min) / (expr_max - expr_min)
        )
        ax.scatter(
            x,
            y,
            s=size,
            c=[plt.cm.viridis(color_value)],
            edgecolor="black",
            linewidth=0.3,
        )

    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(SUBTYPE_ORDER)))
    ax.set_yticklabels([SUBTYPE_LABELS[subtype] for subtype in SUBTYPE_ORDER], fontsize=8)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(SUBTYPE_ORDER) - 0.45, -0.55)
    ax.set_title(
        "Marker support for refined mural subtype labels",
        fontsize=10,
        fontweight="bold",
        pad=8,
    )
    ax.set_xlabel("Canonical marker genes", fontsize=8)
    ax.set_ylabel("")
    ax.grid(color="#E5E5E5", linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    cax = fig.add_axes([0.81, 0.33, 0.016, 0.42])
    norm = matplotlib.colors.Normalize(vmin=expr_min, vmax=expr_max)
    sm = matplotlib.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Mean log1p expression", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    legend_sizes = [10, 50, 90]
    handles = [
        ax.scatter([], [], s=size_min + (pct / 100.0) * (size_max - size_min), color="#BDBDBD", edgecolor="black")
        for pct in legend_sizes
    ]
    ax.legend(
        handles,
        [f"{pct}% detected" for pct in legend_sizes],
        title="Detected cells",
        frameon=False,
        bbox_to_anchor=(1.18, 0.92),
        loc="upper left",
        fontsize=6.5,
        title_fontsize=7,
        borderaxespad=0.0,
    )
    fig.subplots_adjust(left=0.17, right=0.78, bottom=0.33, top=0.82)
    save_figure(fig, "GSE206528_refined_mural_subtype_marker_support_dotplot")


def plot_author_marker_crosscheck() -> None:
    author_support = pd.read_csv(TABLE_DIR / "refined_mural_author_marker_crosscheck_top3.tsv", sep="\t")
    author_support["cluster_label"] = (
        "Leiden "
        + author_support["leiden_r0_5"].astype(str)
        + "\n"
        + author_support["mural_subtype"].map(SUBTYPE_LABELS).fillna(author_support["mural_subtype"])
    )
    cluster_order = ["6", "8", "12", "13"]
    cluster_labels = (
        author_support.drop_duplicates("leiden_r0_5")
        .assign(leiden_r0_5=lambda df: df["leiden_r0_5"].astype(str))
        .set_index("leiden_r0_5")
        .loc[cluster_order, "cluster_label"]
        .tolist()
    )
    author_clusters = ["PC", "Vessel SMC", "CC SMC"]
    colors = {
        "PC": "#009E73",
        "Vessel SMC": "#0072B2",
        "CC SMC": "#CC79A7",
    }

    fig, ax = plt.subplots(figsize=(5.8, 2.9))
    x = np.arange(len(cluster_order))
    width = 0.24
    for i, author_cluster in enumerate(author_clusters):
        values = []
        for cluster in cluster_order:
            sub = author_support.loc[
                (author_support["leiden_r0_5"].astype(str) == cluster)
                & (author_support["author_cluster"] == author_cluster)
            ]
            values.append(float(sub["n_overlap_top100"].iloc[0]) if not sub.empty else 0.0)
        ax.bar(
            x + (i - 1) * width,
            values,
            width=width,
            label=author_cluster,
            color=colors[author_cluster],
            edgecolor="black",
            linewidth=0.35,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(cluster_labels, fontsize=7)
    ax.set_ylabel("Overlap with author top-100 markers", fontsize=8)
    ax.set_title(
        "Author-marker cross-check supports refined mural labels",
        fontsize=10,
        fontweight="bold",
        pad=8,
    )
    ax.legend(frameon=False, fontsize=7, ncols=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.5)
    fig.tight_layout()
    save_figure(fig, "GSE206528_refined_mural_author_marker_crosscheck_top3")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    donor_table = build_donor_candidate_table()
    plot_candidate_donor_expression(
        donor_table,
        CANDIDATE_GENES,
        "GSE206528_diabetic_ED_pericyte_like_12_gene_donor_logcpm",
        ncols=4,
    )
    plot_candidate_donor_expression(
        donor_table,
        MAIN_TEXT_REPRESENTATIVE_GENES,
        "GSE206528_diabetic_ED_pericyte_like_representative_gene_donor_logcpm",
        ncols=5,
    )

    marker_table = build_marker_support_table()
    plot_marker_support(marker_table)
    plot_author_marker_crosscheck()
    print(f"Wrote marker-support figures to {FIG_DIR.relative_to(ROOT)}")
    print(f"Wrote marker-support tables to {TABLE_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
