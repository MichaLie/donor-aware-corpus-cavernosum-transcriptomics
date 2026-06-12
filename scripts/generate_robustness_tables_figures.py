#!/usr/bin/env python3
"""Generate robustness tables and figures.

This script prepares selected module gene sets, donor-level composition tables,
module-score sensitivity plots, and redesigned candidate/validation figures.
The limma cameraPR analyses are run separately by
scripts/r/gse206528_robustness_cameraPR.R after this script writes the gene-set
file, then this script can be rerun to build all figures from the R outputs.
"""

from __future__ import annotations

from pathlib import Path

import gseapy as gp
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme(style="white", context="paper", font_scale=0.9)
plt.rcParams.update(
    {
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.titlesize": 10,
    }
)


ROOT = Path(__file__).resolve().parents[1]
GENE_SET_OUT = ROOT / "metadata/gene_sets/robustness_module_gene_sets.tsv"
TABLE_DIR = ROOT / "results/tables/robustness_checks"
FIG_DIR = ROOT / "results/figures/robustness_checks"

PSEUDOBULK_COUNTS = ROOT / "data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
PSEUDOBULK_META = ROOT / "data/processed/GSE206528/pseudobulk_mural_subtypes/GSE206528_mural_subtype_pseudobulk_metadata.tsv"
EDGER_RESULTS = ROOT / "results/tables/GSE206528_mural_subtype_edger/all_mural_subtype_edgeR_results.tsv.gz"
DONOR_CANDIDATE_TABLE = ROOT / "results/tables/GSE206528_refined_mural_marker_support/diabetic_ed_pericyte_like_12_gene_donor_logcpm.tsv"

DONOR_ORDER = ["N1", "N2", "N3", "ED1", "ED2", "ED3", "DMED1", "DMED2"]
DISCOVERY_DONOR_ORDER = ["N1", "N2", "N3", "DMED1", "DMED2"]
GROUP_ORDER = ["normal", "non_diabetic_ED", "diabetic_ED"]
GROUP_LABELS = {
    "normal": "Normal",
    "non_diabetic_ED": "Non-diabetic ED",
    "diabetic_ED": "Diabetic ED",
}
GROUP_COLORS = {
    "normal": "#0072B2",
    "non_diabetic_ED": "#999999",
    "diabetic_ED": "#D55E00",
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
CANDIDATE_UP_GENES = [g for g in CANDIDATE_GENES if g != "ERAP2"]
PRIMARY_STORY_GENES = ["IL6", "ADAMTS4", "PLA2G2A", "CYR61", "PI15"]
CANDIDATE_DIRECTION = {g: 1 for g in CANDIDATE_UP_GENES} | {"ERAP2": -1}

SELECTED_HALLMARKS = {
    "TNF-alpha Signaling via NF-kB": ("hallmark_tnfa_nfkb", "Hallmark TNF-alpha/NF-kB signaling", "immune"),
    "Inflammatory Response": ("hallmark_inflammatory_response", "Hallmark inflammatory response", "immune"),
    "IL-6/JAK/STAT3 Signaling": ("hallmark_il6_jak_stat3", "Hallmark IL6/JAK/STAT3 signaling", "immune"),
    "Epithelial Mesenchymal Transition": ("hallmark_emt_ecm", "Hallmark epithelial-mesenchymal/ECM remodeling", "remodeling"),
    "Angiogenesis": ("hallmark_angiogenesis", "Hallmark angiogenesis", "vascular"),
    "Hypoxia": ("hallmark_hypoxia", "Hallmark hypoxia", "metabolism"),
    "Reactive Oxygen Species Pathway": ("hallmark_ros", "Hallmark reactive oxygen species pathway", "metabolism"),
    "Unfolded Protein Response": ("hallmark_unfolded_protein_response", "Hallmark unfolded protein response", "stress"),
    "Myogenesis": ("hallmark_myogenesis", "Hallmark myogenesis", "contractility"),
    "Complement": ("hallmark_complement", "Hallmark complement", "immune"),
}

MODULE_SCORE_ORDER = [
    "candidate_directional_signature_12",
    "candidate_inflammatory_remodeling_5",
    "immediate_early_stress",
    "hallmark_inflammatory_response",
    "hallmark_il6_jak_stat3",
    "hallmark_tnfa_nfkb",
    "ecm_fibrosis",
    "hallmark_emt_ecm",
    "hallmark_angiogenesis",
    "hallmark_hypoxia",
    "hallmark_ros",
    "smooth_muscle_contractility",
]


def ensure_dirs() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    GENE_SET_OUT.parent.mkdir(parents=True, exist_ok=True)


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=450, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def style_heatmap_colorbar(ax: plt.Axes, label: str) -> None:
    cbar = ax.collections[0].colorbar
    if cbar is not None:
        cbar.set_label(label, fontsize=7)
        cbar.ax.tick_params(labelsize=7)


def clean_gene_list(genes: list[str] | str) -> list[str]:
    if isinstance(genes, str):
        genes = genes.split(",")
    out = []
    for gene in genes:
        gene = str(gene).strip().upper()
        if gene and gene not in out:
            out.append(gene)
    return out


def prepare_gene_sets() -> pd.DataFrame:
    curated = pd.read_csv(ROOT / "metadata/gene_sets/curated_disease_programs.tsv", sep="\t")
    rows = curated.to_dict("records")

    rows.extend(
        [
            {
                "program_id": "candidate_signature_12",
                "program_label": "Diabetic ED pericyte-like 12-gene candidate signature",
                "program_family": "candidate_signature",
                "source_status": "analysis_result",
                "gene_symbols": ",".join(CANDIDATE_GENES),
                "rationale": "All FDR < 0.05 diabetic ED pericyte-like candidates, including ERAP2 as the single downregulated candidate.",
            },
            {
                "program_id": "candidate_up_signature_11",
                "program_label": "Diabetic ED pericyte-like 11-gene up-signature",
                "program_family": "candidate_signature",
                "source_status": "analysis_result",
                "gene_symbols": ",".join(CANDIDATE_UP_GENES),
                "rationale": "Upregulated diabetic ED pericyte-like candidates used for same-direction Peyronie's cameraPR testing.",
            },
            {
                "program_id": "candidate_inflammatory_remodeling_5",
                "program_label": "Inflammatory-remodeling candidate genes",
                "program_family": "candidate_signature",
                "source_status": "analysis_interpretation",
                "gene_symbols": ",".join(PRIMARY_STORY_GENES),
                "rationale": "Most interpretable candidate genes for the biological interpretation.",
            },
        ]
    )

    try:
        hallmark = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
        for hallmark_name, (program_id, label, family) in SELECTED_HALLMARKS.items():
            if hallmark_name not in hallmark:
                continue
            rows.append(
                {
                    "program_id": program_id,
                    "program_label": label,
                    "program_family": family,
                    "source_status": "MSigDB_Hallmark_2020_via_Enrichr",
                    "gene_symbols": ",".join(clean_gene_list(hallmark[hallmark_name])),
                    "rationale": f"Established MSigDB Hallmark gene set: {hallmark_name}.",
                }
            )
    except Exception as exc:
        print(f"WARNING: could not download MSigDB Hallmark library via gseapy: {exc}")

    gene_sets = pd.DataFrame(rows).drop_duplicates("program_id", keep="last")
    gene_sets["gene_symbols"] = gene_sets["gene_symbols"].map(lambda x: ",".join(clean_gene_list(x)))
    gene_sets["n_genes"] = gene_sets["gene_symbols"].str.split(",").map(len)
    gene_sets = gene_sets.sort_values(["program_family", "program_id"])
    gene_sets.to_csv(GENE_SET_OUT, sep="\t", index=False)
    return gene_sets


def load_logcpm() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(PSEUDOBULK_COUNTS, sep="\t", index_col=0)
    counts.index = counts.index.str.upper()
    meta = pd.read_csv(PSEUDOBULK_META, sep="\t").set_index("pseudobulk_id")
    counts = counts.loc[:, meta.index]
    library_sizes = counts.sum(axis=0)
    log_cpm = np.log2((counts + 0.5).div(library_sizes + 1, axis=1) * 1_000_000)
    return log_cpm, meta


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    means = df.mean(axis=1)
    stds = df.std(axis=1, ddof=0).replace(0, np.nan)
    return df.sub(means, axis=0).div(stds, axis=0).fillna(0)


def candidate_tier_table() -> pd.DataFrame:
    tier_map = {
        "IL6": ("Inflammatory-remodeling", "Inflammatory cytokine axis"),
        "ADAMTS4": ("Inflammatory-remodeling", "Matrix-remodeling/protease axis"),
        "PLA2G2A": ("Inflammatory-remodeling", "Inflammatory/lipid mediator-related signal"),
        "CYR61": ("Inflammatory-remodeling", "Tissue remodeling or angiogenic-stress context"),
        "PI15": ("Inflammatory-remodeling", "Pericyte-like candidate; broader organic ED signal"),
        "NR4A1": ("Activation/stress context", "Immediate-early/stress-activation gene"),
        "NR4A2": ("Activation/stress context", "Immediate-early/stress-activation gene"),
        "NR4A3": ("Activation/stress context", "Immediate-early/stress-activation gene"),
        "MT1A": ("Activation/stress context", "Stress-response candidate"),
        "RGS16": ("Activation/stress context", "Signaling regulator; activation/stress context"),
        "SH3BP5": ("Exploratory/uncertain", "Candidate signal; limited mechanistic interpretation"),
        "ERAP2": ("Exploratory/uncertain", "Downregulated candidate; cautious interpretation"),
    }
    edger = pd.read_csv(EDGER_RESULTS, sep="\t")
    edger["gene_symbol"] = edger["gene_symbol"].str.upper()
    stats = edger.loc[
        (edger["mural_subtype"] == "Pericyte_like")
        & (edger["comparison"] == "diabetic_ED_vs_normal")
        & (edger["gene_symbol"].isin(CANDIDATE_GENES)),
        ["gene_symbol", "logFC", "FDR"],
    ].copy()
    rows = []
    for gene in CANDIDATE_GENES:
        row = stats.loc[stats["gene_symbol"] == gene].iloc[0]
        tier, rationale = tier_map[gene]
        rows.append(
            {
                "gene_symbol": gene,
                "tier": tier,
                "direction": "up_in_diabetic_ED" if row["logFC"] > 0 else "down_in_diabetic_ED",
                "logFC": row["logFC"],
                "FDR": row["FDR"],
                "rationale": rationale,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "diabetic_ed_candidate_gene_tiers.tsv", sep="\t", index=False)
    return out


def build_cell_composition_tables_and_figure() -> None:
    sample_meta = pd.read_csv(ROOT / "metadata/GSE206528/samples.tsv", sep="\t")
    donor_group = sample_meta.set_index("donor_id")["ed_subgroup"].to_dict()

    major = pd.read_csv(ROOT / "results/tables/GSE206528_major_celltype_counts_by_donor.tsv", sep="\t")
    total_cells = major.groupby("donor_id")["n_cells"].sum()
    major["total_annotated_cells"] = major["donor_id"].map(total_cells)
    major["proportion_of_total"] = major["n_cells"] / major["total_annotated_cells"]
    major["ed_subgroup"] = major["donor_id"].map(donor_group)
    major.to_csv(TABLE_DIR / "GSE206528_major_celltype_composition_by_donor.tsv", sep="\t", index=False)

    refined = pd.read_csv(ROOT / "results/tables/GSE206528_mural_subtype_counts_by_donor.tsv", sep="\t")
    refined_total = refined.groupby("donor_id")["n_cells"].sum()
    refined["total_annotated_cells"] = refined["donor_id"].map(total_cells)
    refined["total_refined_mural_smc_cells"] = refined["donor_id"].map(refined_total)
    refined["proportion_of_total_cells"] = refined["n_cells"] / refined["total_annotated_cells"]
    refined["proportion_within_refined_mural_smc"] = refined["n_cells"] / refined["total_refined_mural_smc_cells"]
    refined["ed_subgroup"] = refined["donor_id"].map(donor_group)
    refined.to_csv(TABLE_DIR / "GSE206528_refined_mural_smc_composition_by_donor.tsv", sep="\t", index=False)

    cell_order = ["Endothelial", "Fibroblast", "Mural", "Smooth_muscle", "Macrophage", "T_cell", "Schwann_neural"]
    cell_colors = {
        "Endothelial": "#0072B2",
        "Fibroblast": "#009E73",
        "Mural": "#D55E00",
        "Smooth_muscle": "#CC79A7",
        "Macrophage": "#E69F00",
        "T_cell": "#56B4E9",
        "Schwann_neural": "#999999",
    }
    subtype_order = ["Pericyte_like", "Vessel_SMC_like", "CC_SMC_like"]
    subtype_colors = {
        "Pericyte_like": "#D55E00",
        "Vessel_SMC_like": "#0072B2",
        "CC_SMC_like": "#CC79A7",
    }

    fig = plt.figure(figsize=(9.6, 6.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.0], hspace=0.60, wspace=0.36)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])

    bottom = np.zeros(len(DONOR_ORDER))
    for cell in cell_order:
        values = (
            major.loc[major["major_cell_type"] == cell]
            .set_index("donor_id")
            .reindex(DONOR_ORDER)["proportion_of_total"]
            .fillna(0)
            .to_numpy()
        )
        ax1.bar(DONOR_ORDER, values, bottom=bottom, label=cell.replace("_", " "), color=cell_colors[cell], width=0.74)
        bottom += values
    ax1.set_ylabel("Fraction of annotated cells")
    ax1.set_title("A. Major cell types", loc="left", fontweight="bold", pad=6)
    ax1.legend(ncol=2, frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.02))
    ax1.tick_params(axis="x", rotation=0)
    for tick, donor in zip(ax1.get_xticklabels(), DONOR_ORDER):
        tick.set_color(GROUP_COLORS[donor_group[donor]])

    bottom = np.zeros(len(DONOR_ORDER))
    for subtype in subtype_order:
        values = (
            refined.loc[refined["mural_subtype"] == subtype]
            .set_index("donor_id")
            .reindex(DONOR_ORDER)["proportion_within_refined_mural_smc"]
            .fillna(0)
            .to_numpy()
        )
        ax2.bar(DONOR_ORDER, values, bottom=bottom, label=subtype.replace("_", " "), color=subtype_colors[subtype], width=0.74)
        bottom += values
    ax2.set_ylabel("Fraction within refined mural/SMC")
    ax2.set_title("B. Refined mural/SMC", loc="left", fontweight="bold", pad=6)
    ax2.tick_params(axis="x", labelrotation=45)
    ax2.legend(frameon=False, fontsize=7)

    pericyte = refined.loc[refined["mural_subtype"] == "Pericyte_like"].copy()
    x_map = {"normal": 0, "non_diabetic_ED": 1, "diabetic_ED": 2}
    for subgroup in GROUP_ORDER:
        sub = pericyte.loc[pericyte["ed_subgroup"] == subgroup]
        x = np.repeat(x_map[subgroup], len(sub)) + np.linspace(-0.08, 0.08, max(len(sub), 1))[: len(sub)]
        ax3.scatter(
            x,
            sub["proportion_of_total_cells"],
            color=GROUP_COLORS[subgroup],
            edgecolor="black",
            linewidth=0.4,
            s=52,
            label=GROUP_LABELS[subgroup],
            zorder=3,
        )
        for xi, (_, row) in zip(x, sub.iterrows()):
            ax3.text(xi, row["proportion_of_total_cells"] + 0.0007, row["donor_id"], ha="center", va="bottom", fontsize=7)
    ax3.set_xticks([0, 1, 2], [GROUP_LABELS[g] for g in GROUP_ORDER], rotation=20, ha="right")
    ax3.set_ylabel("Pericyte-like fraction of total cells")
    ax3.set_title("C. Pericyte-like abundance", loc="left", fontweight="bold", pad=6)
    ax3.set_ylim(0.0045, pericyte["proportion_of_total_cells"].max() + 0.0045)
    ax3.grid(axis="y", color="#E5E7EB")
    sns.despine(fig=fig)
    save_figure(fig, "GSE206528_descriptive_cell_composition_by_donor")


def build_donor_module_scores(gene_sets: pd.DataFrame) -> pd.DataFrame:
    log_cpm, meta = load_logcpm()
    selected_meta = meta.loc[
        (meta["mural_subtype"] == "Pericyte_like")
        & (meta["n_cells"] >= 50)
        & (meta["ed_subgroup"].isin(["normal", "diabetic_ED"]))
    ].copy()
    selected_meta = selected_meta.loc[[f"{donor}__Pericyte_like" for donor in DISCOVERY_DONOR_ORDER]]
    selected = log_cpm.loc[:, selected_meta.index]
    z = zscore_rows(selected)

    gs_lookup = {
        row.program_id: clean_gene_list(row.gene_symbols)
        for row in gene_sets.itertuples(index=False)
    }
    rows: list[dict[str, object]] = []
    for module_id in MODULE_SCORE_ORDER:
        if module_id == "candidate_directional_signature_12":
            genes = [g for g in CANDIDATE_GENES if g in z.index]
            score = pd.Series(0.0, index=selected_meta.index)
            for gene in genes:
                score = score.add(z.loc[gene] * CANDIDATE_DIRECTION[gene], fill_value=0)
            score = score / max(1, len(genes))
            label = "Candidate directional signature"
            family = "candidate_signature"
        else:
            genes = [g for g in gs_lookup.get(module_id, []) if g in z.index]
            if not genes:
                continue
            score = z.loc[genes].mean(axis=0)
            label = gene_sets.loc[gene_sets["program_id"] == module_id, "program_label"].iloc[0]
            family = gene_sets.loc[gene_sets["program_id"] == module_id, "program_family"].iloc[0]
        for pseudobulk_id, value in score.items():
            meta_row = selected_meta.loc[pseudobulk_id]
            rows.append(
                {
                    "module_id": module_id,
                    "module_label": label,
                    "program_family": family,
                    "pseudobulk_id": pseudobulk_id,
                    "donor_id": meta_row["donor_id"],
                    "ed_subgroup": meta_row["ed_subgroup"],
                    "n_genes_scored": len(genes),
                    "genes_scored": ",".join(genes),
                    "mean_z_score": value,
                }
            )
    scores = pd.DataFrame(rows)
    scores.to_csv(TABLE_DIR / "GSE206528_pericyte_like_donor_module_scores.tsv", sep="\t", index=False)
    return scores


def plot_candidate_heatmap(tiers: pd.DataFrame) -> None:
    table = pd.read_csv(DONOR_CANDIDATE_TABLE, sep="\t")
    pivot = table.pivot_table(index="gene_symbol", columns="donor_id", values="log2_cpm", aggfunc="first")
    pivot = pivot.reindex(index=CANDIDATE_GENES, columns=DISCOVERY_DONOR_ORDER)
    z = zscore_rows(pivot)
    tier_order = ["Inflammatory-remodeling", "Activation/stress context", "Exploratory/uncertain"]
    sorted_genes = []
    for tier in tier_order:
        sorted_genes.extend(tiers.loc[tiers["tier"] == tier, "gene_symbol"].tolist())
    z = z.reindex(sorted_genes)

    fig, ax = plt.subplots(figsize=(5.4, 5.2))
    sns.heatmap(
        z,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-2,
        vmax=2,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"shrink": 0.76, "pad": 0.03},
    )
    style_heatmap_colorbar(ax, "Donor row z-score\n(log2 CPM)")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Diabetic ED pericyte-like candidates", fontsize=9, fontweight="bold", pad=8)
    for tick in ax.get_xticklabels():
        tick.set_rotation(0)
        tick.set_color(GROUP_COLORS["normal"] if tick.get_text().startswith("N") else GROUP_COLORS["diabetic_ED"])
        tick.set_fontweight("bold")
    ax.tick_params(axis="y", pad=6)
    for y in [5, 10]:
        ax.hlines(y, *ax.get_xlim(), color="black", linewidth=1.2)
    ax.text(-1.08, 2.5, "Primary\nstory", va="center", ha="right", fontsize=7, clip_on=False)
    ax.text(-1.08, 7.5, "Activation/\nstress", va="center", ha="right", fontsize=7, clip_on=False)
    ax.text(-1.08, 11, "Exploratory", va="center", ha="right", fontsize=7, clip_on=False)
    save_figure(fig, "candidate_donor_heatmap")


def plot_module_sensitivity(scores: pd.DataFrame) -> None:
    module_results_path = TABLE_DIR / "GSE206528_refined_mural_module_cameraPR_results.tsv"
    if not module_results_path.exists():
        print("Module cameraPR results not found yet; skipping module sensitivity plot.")
        return
    module_results = pd.read_csv(module_results_path, sep="\t")
    context = module_results.loc[
        (module_results["mural_subtype"] == "Pericyte_like")
        & (module_results["comparison"] == "diabetic_ED_vs_normal")
        & (module_results["program_id"].isin([m for m in MODULE_SCORE_ORDER if m != "candidate_directional_signature_12"]))
    ].copy()
    context["neg_log10_FDR"] = -np.log10(context["camera_FDR_within_context"].clip(lower=1e-300))
    context["signed_neg_log10_FDR"] = np.where(context["camera_direction"] == "Up", 1, -1) * context["neg_log10_FDR"]
    context["program_id"] = pd.Categorical(context["program_id"], [m for m in MODULE_SCORE_ORDER if m != "candidate_directional_signature_12"], ordered=True)
    context = context.sort_values("program_id")

    score_pivot = scores.pivot_table(index="module_id", columns="donor_id", values="mean_z_score", aggfunc="first")
    score_pivot = score_pivot.reindex(index=[m for m in MODULE_SCORE_ORDER if m in score_pivot.index], columns=DISCOVERY_DONOR_ORDER)
    label_lookup = scores.drop_duplicates("module_id").set_index("module_id")["module_label"].to_dict()
    short_score_labels = {
        "candidate_directional_signature_12": "Candidate signature",
        "candidate_inflammatory_remodeling_5": "Inflammatory-remodeling 5 genes",
        "immediate_early_stress": "Immediate-early/stress",
        "hallmark_inflammatory_response": "Inflammatory response",
        "hallmark_il6_jak_stat3": "IL6/JAK/STAT3",
        "hallmark_tnfa_nfkb": "TNF/NF-kB",
        "ecm_fibrosis": "ECM/fibrosis",
        "hallmark_emt_ecm": "EMT/ECM remodeling",
        "hallmark_angiogenesis": "Angiogenesis",
        "hallmark_hypoxia": "Hypoxia",
        "hallmark_ros": "ROS",
        "smooth_muscle_contractility": "Smooth-muscle contractility",
    }
    score_pivot.index = [short_score_labels.get(i, label_lookup.get(i, i)) for i in score_pivot.index]

    fig = plt.figure(figsize=(10.4, 6.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.50)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    sns.heatmap(
        score_pivot,
        ax=ax1,
        cmap="vlag",
        center=0,
        vmin=-1.8,
        vmax=1.8,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"shrink": 0.72, "pad": 0.03},
    )
    style_heatmap_colorbar(ax1, "Mean gene z-score")
    ax1.set_title("A. Donor module scores", loc="left", fontweight="bold", pad=8)
    ax1.set_xlabel("")
    ax1.set_ylabel("")
    for tick in ax1.get_xticklabels():
        tick.set_rotation(0)
        tick.set_color(GROUP_COLORS["normal"] if tick.get_text().startswith("N") else GROUP_COLORS["diabetic_ED"])
        tick.set_fontweight("bold")

    plot_context = context.dropna(subset=["signed_neg_log10_FDR"]).copy()
    module_label_map = {
        "candidate_inflammatory_remodeling_5": "Inflammatory-remodeling 5 genes",
        "immediate_early_stress": "Immediate-early/stress",
        "hallmark_inflammatory_response": "Inflammatory response",
        "hallmark_il6_jak_stat3": "IL6/JAK/STAT3",
        "hallmark_tnfa_nfkb": "TNF/NF-kB",
        "ecm_fibrosis": "ECM/fibrosis",
        "hallmark_emt_ecm": "EMT/ECM remodeling",
        "hallmark_angiogenesis": "Angiogenesis",
        "hallmark_hypoxia": "Hypoxia",
        "hallmark_ros": "ROS",
        "smooth_muscle_contractility": "SM contractility",
    }
    labels = plot_context["program_id"].astype(str).map(module_label_map).fillna(plot_context["program_label"])
    colors = np.where(plot_context["signed_neg_log10_FDR"] >= 0, "#D55E00", "#0072B2")
    ax2.barh(range(len(plot_context)), plot_context["signed_neg_log10_FDR"], color=colors)
    ax2.axvline(0, color="black", linewidth=0.8)
    ax2.axvline(-np.log10(0.05), color="#999999", linestyle="--", linewidth=0.8)
    ax2.axvline(np.log10(0.05) * -1, color="#999999", linestyle="--", linewidth=0.8)
    ax2.set_yticks(range(len(plot_context)), labels, fontsize=7)
    ax2.set_xlabel("Signed -log10 cameraPR FDR\n(positive = up in diabetic ED)")
    ax2.set_title("B. Ranked enrichment", loc="left", fontweight="bold", pad=8)
    ax2.grid(axis="x", color="#E5E7EB")
    sns.despine(fig=fig)
    save_figure(fig, "GSE206528_pericyte_like_module_stress_sensitivity")


def plot_peyronie_robustness() -> None:
    loo_path = TABLE_DIR / "peyronie_leave_one_gene_out_cameraPR.tsv"
    rank_path = TABLE_DIR / "peyronie_candidate_gene_rank_table.tsv"
    if not loo_path.exists() or not rank_path.exists():
        print("Peyronie robustness R outputs not found yet; skipping Peyronie robustness plot.")
        return
    loo = pd.read_csv(loo_path, sep="\t")
    ranks = pd.read_csv(rank_path, sep="\t")

    disease_comps = ["PD_plaque_vs_control", "PF_vs_NF"]
    disease = loo.loc[loo["peyronie_comparison"].isin(disease_comps)].copy()
    label_map = {
        "PD_plaque_vs_control": "PD plaque vs control\nGSE126005",
        "PF_vs_NF": "PD fibroblast vs normal\nGSE146500",
        "PFPC_vs_PF": "PC-NV PF vs PF\nGSE146500",
        "PFPC_vs_NF": "PC-NV PF vs normal\nGSE146500",
    }
    disease["label"] = disease["peyronie_comparison"].map(label_map)
    disease["neg_log10_FDR"] = -np.log10(disease["camera_FDR_within_comparison"].clip(lower=1e-300))

    rank_subset = ranks.loc[ranks["gene_symbol"].isin(CANDIDATE_UP_GENES)].copy()
    rank_subset["contrast_label"] = rank_subset["peyronie_comparison"].map(label_map)
    contrast_order = [label_map[x] for x in ["PD_plaque_vs_control", "PF_vs_NF", "PFPC_vs_PF", "PFPC_vs_NF"]]
    heat = rank_subset.pivot_table(index="gene_symbol", columns="contrast_label", values="signed_sqrt_F", aggfunc="first")
    heat = heat.reindex(index=CANDIDATE_UP_GENES, columns=contrast_order)

    fig = plt.figure(figsize=(10.2, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.25], wspace=0.38)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    sns.stripplot(
        data=disease,
        x="label",
        y="neg_log10_FDR",
        hue="omitted_gene",
        dodge=False,
        ax=ax1,
        palette="tab20",
        size=5,
        edgecolor="black",
        linewidth=0.25,
        legend=False,
    )
    full = disease.loc[disease["omitted_gene"] == "NONE_FULL_SIGNATURE"]
    ax1.scatter(full["label"], full["neg_log10_FDR"], marker="D", s=70, color="black", label="Full signature", zorder=4)
    ax1.axhline(-np.log10(0.05), color="#999999", linestyle="--", linewidth=0.9)
    ax1.set_ylabel("-log10 cameraPR FDR")
    ax1.set_xlabel("")
    ax1.set_title("A. Leave-one-gene-out cameraPR", loc="left", fontweight="bold", pad=8)
    ax1.tick_params(axis="x", labelrotation=20)
    ax1.grid(axis="y", color="#E5E7EB")

    sns.heatmap(
        heat,
        ax=ax2,
        cmap="vlag",
        center=0,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"shrink": 0.82, "pad": 0.05},
    )
    style_heatmap_colorbar(ax2, "Peyronie's signed sqrt(F)")
    ax2.set_title("B. Per-gene Peyronie's direction", loc="left", fontweight="bold", pad=8)
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    ax2.tick_params(axis="x", labelrotation=25)
    sns.despine(fig=fig)
    save_figure(fig, "peyronie_signature_robustness")


def main() -> None:
    ensure_dirs()
    gene_sets = prepare_gene_sets()
    tiers = candidate_tier_table()
    build_cell_composition_tables_and_figure()
    scores = build_donor_module_scores(gene_sets)
    plot_candidate_heatmap(tiers)
    plot_module_sensitivity(scores)
    plot_peyronie_robustness()
    print(f"Wrote tables to {TABLE_DIR}")
    print(f"Wrote figures to {FIG_DIR}")
    print(f"Wrote gene sets to {GENE_SET_OUT}")


if __name__ == "__main__":
    main()
