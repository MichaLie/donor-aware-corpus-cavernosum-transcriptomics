from __future__ import annotations

import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures"
OUT_TABLE = TABLE_DIR / "reproducibility_checks"
OUT_FIG = FIG_DIR / "reproducibility_checks"

CANDIDATE_GENES = [
    "IL6",
    "ADAMTS4",
    "PLA2G2A",
    "CYR61",
    "PI15",
    "NR4A1",
    "NR4A2",
    "NR4A3",
    "MT1A",
    "RGS16",
    "SH3BP5",
    "ERAP2",
]
PRIMARY_GENES = ["IL6", "ADAMTS4", "PLA2G2A", "CYR61", "PI15"]
DOWN_IN_DMED = {"ERAP2"}
MIN_CELLS_FORMAL = 50


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs)


def write_table(df: pd.DataFrame, name: str) -> None:
    OUT_TABLE.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TABLE / f"{name}.tsv", sep="\t", index=False)
    df.to_csv(OUT_TABLE / f"{name}.csv", index=False)


def safe_z(values: pd.Series, center: float, scale: float) -> pd.Series:
    if scale == 0 or np.isnan(scale):
        return values * 0
    return (values - center) / scale


def build_exact_permutation() -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = read_tsv(TABLE_DIR / "robustness_checks" / "GSE206528_pericyte_like_donor_module_scores.tsv")
    scores = scores[scores["module_id"].eq("candidate_directional_signature_12")].copy()
    scores = scores.sort_values("donor_id")

    donor_scores = dict(zip(scores["donor_id"], scores["mean_z_score"]))
    true_dmed = {"DMED1", "DMED2"}
    donors = list(donor_scores)
    rows = []
    for pair in itertools.combinations(donors, 2):
        pair = set(pair)
        rest = set(donors) - pair
        diff = float(np.mean([donor_scores[d] for d in pair]) - np.mean([donor_scores[d] for d in rest]))
        rows.append(
            {
                "putative_diabetic_pair": ",".join(sorted(pair)),
                "remaining_reference_triplet": ",".join(sorted(rest)),
                "mean_putative_diabetic_score": float(np.mean([donor_scores[d] for d in pair])),
                "mean_remaining_reference_score": float(np.mean([donor_scores[d] for d in rest])),
                "score_difference": diff,
                "is_observed_diabetic_assignment": pair == true_dmed,
            }
        )
    out = pd.DataFrame(rows).sort_values("score_difference", ascending=False).reset_index(drop=True)
    out["rank_descending"] = np.arange(1, len(out) + 1)
    observed = out[out["is_observed_diabetic_assignment"]].iloc[0]
    one_sided = (out["score_difference"] >= observed["score_difference"] - 1e-12).sum() / len(out)
    two_sided = (out["score_difference"].abs() >= abs(observed["score_difference"]) - 1e-12).sum() / len(out)
    summary = pd.DataFrame(
        [
            {
                "observed_diabetic_pair": "DMED1,DMED2",
                "n_total_2_vs_3_assignments": len(out),
                "observed_score_difference": observed["score_difference"],
                "observed_rank_descending": int(observed["rank_descending"]),
                "exact_one_sided_p_ge_observed": one_sided,
                "exact_two_sided_p_abs_ge_observed": two_sided,
                "interpretation": "Sanity check only: observed diabetic labels have the largest 12-gene directional-score separation among the ten possible 2-vs-3 assignments.",
            }
        ]
    )
    return out, summary


def build_descriptive_scores() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = read_tsv(
        ROOT / "data" / "processed" / "GSE206528" / "pseudobulk_mural_subtypes" / "GSE206528_mural_subtype_pseudobulk_counts.tsv.gz"
    )
    meta = read_tsv(
        ROOT / "data" / "processed" / "GSE206528" / "pseudobulk_mural_subtypes" / "GSE206528_mural_subtype_pseudobulk_metadata.tsv"
    )
    meta = meta[meta["mural_subtype"].eq("Pericyte_like")].copy()
    pericyte_ids = meta["pseudobulk_id"].tolist()

    gene_counts = counts[counts["gene_symbol"].isin(CANDIDATE_GENES)][["gene_symbol"] + pericyte_ids].copy()
    long = gene_counts.melt(id_vars="gene_symbol", var_name="pseudobulk_id", value_name="count")
    long = long.merge(meta, on="pseudobulk_id", how="left")
    lib_sizes = counts[pericyte_ids].sum(axis=0).rename("library_size").reset_index().rename(columns={"index": "pseudobulk_id"})
    long = long.merge(lib_sizes, on="pseudobulk_id", how="left")
    long["log2_cpm"] = np.log2(((long["count"] + 0.5) / (long["library_size"] + 1.0)) * 1_000_000)

    normal = long[long["ed_subgroup"].eq("normal")]
    stats = (
        normal.groupby("gene_symbol")["log2_cpm"]
        .agg(normal_reference_mean="mean", normal_reference_sd=lambda x: x.std(ddof=1))
        .reset_index()
    )
    long = long.merge(stats, on="gene_symbol", how="left")
    long["normal_reference_z"] = np.where(
        (long["normal_reference_sd"].eq(0)) | (long["normal_reference_sd"].isna()),
        0,
        (long["log2_cpm"] - long["normal_reference_mean"]) / long["normal_reference_sd"],
    )
    long["directional_z"] = np.where(long["gene_symbol"].isin(DOWN_IN_DMED), -long["normal_reference_z"], long["normal_reference_z"])
    long["passes_primary_min_cell_threshold"] = long["n_cells"] >= MIN_CELLS_FORMAL

    module_rows = []
    for genes, module_id, label in [
        (CANDIDATE_GENES, "candidate_directional_signature_12", "12-gene directional score"),
        (PRIMARY_GENES, "candidate_inflammatory_remodeling_5", "5-gene inflammatory-remodeling score"),
    ]:
        sub = long[long["gene_symbol"].isin(genes)]
        score = (
            sub.groupby(["donor_id", "pseudobulk_id", "ed_subgroup", "n_cells", "passes_primary_min_cell_threshold"], as_index=False)
            .agg(mean_directional_reference_z=("directional_z", "mean"), n_genes_scored=("gene_symbol", "nunique"))
        )
        score["module_id"] = module_id
        score["module_label"] = label
        module_rows.append(score)
    scores = pd.concat(module_rows, ignore_index=True)
    display = {"normal": "Tumor-margin reference", "non_diabetic_ED": "Non-diabetic organic ED", "diabetic_ED": "Diabetic ED"}
    scores["display_group"] = scores["ed_subgroup"].map(display).fillna(scores["ed_subgroup"])
    group_summary = (
        scores.groupby(["module_id", "module_label", "ed_subgroup", "display_group"], as_index=False)
        .agg(
            n_donor_profiles=("donor_id", "nunique"),
            n_profiles_passing_primary_min_cell_threshold=("passes_primary_min_cell_threshold", "sum"),
            mean_score=("mean_directional_reference_z", "mean"),
            median_score=("mean_directional_reference_z", "median"),
            min_score=("mean_directional_reference_z", "min"),
            max_score=("mean_directional_reference_z", "max"),
        )
    )
    return long, scores, group_summary


def build_all_tests_table() -> pd.DataFrame:
    major = read_tsv(TABLE_DIR / "GSE206528_pseudobulk_edger" / "summary_pseudobulk_edgeR_results.tsv")
    major = major.rename(columns={"cell_type": "compartment_or_system"})
    major["dataset"] = "GSE206528"
    major["analysis_level"] = "major cell type pseudobulk"
    major["donor_or_sample_profiles"] = major["n_samples"]
    major["multiple_testing_family"] = "BH FDR across genes within each major-cell-type x contrast edgeR result"

    mural = read_tsv(TABLE_DIR / "GSE206528_mural_subtype_edger" / "summary_mural_subtype_edgeR_results.tsv")
    mural = mural.rename(columns={"mural_subtype": "compartment_or_system"})
    mural["dataset"] = "GSE206528"
    mural["analysis_level"] = "refined mural/SMC subtype pseudobulk"
    mural["donor_or_sample_profiles"] = mural["n_samples"]
    mural["multiple_testing_family"] = "BH FDR across genes within each refined-mural-subtype x contrast edgeR result"

    pey1 = read_tsv(TABLE_DIR / "GSE126005" / "GSE126005_edgeR_summary.tsv")
    pey1["compartment_or_system"] = "bulk tunica albuginea/plaque"
    pey1["analysis_level"] = "bulk RNA-seq edgeR"
    pey1["group_summary"] = "control=6;Peyronie_plaque=6"
    pey1["donor_or_sample_profiles"] = pey1["n_samples"]
    pey1["multiple_testing_family"] = "BH FDR across genes within GSE126005 edgeR contrast"

    pey2 = read_tsv(TABLE_DIR / "GSE146500" / "GSE146500_edgeR_summary.tsv")
    pey2["compartment_or_system"] = "cultured tunica albuginea fibroblast"
    pey2["analysis_level"] = "bulk RNA-seq edgeR"
    group_map = {
        "PF_vs_NF": "normal_fibroblast=4;Peyronie_fibroblast=4",
        "PFPC_vs_PF": "Peyronie_fibroblast=4;PC-NV_treated_Peyronie_fibroblast=4",
        "PFPC_vs_NF": "normal_fibroblast=4;PC-NV_treated_Peyronie_fibroblast=4",
    }
    pey2["group_summary"] = pey2["comparison"].map(group_map)
    pey2["donor_or_sample_profiles"] = pey2["n_samples"]
    pey2["multiple_testing_family"] = "BH FDR across genes within each GSE146500 edgeR contrast"

    all_tests = pd.concat([major, mural, pey1, pey2], ignore_index=True, sort=False)
    keep = [
        "dataset",
        "analysis_level",
        "compartment_or_system",
        "comparison",
        "group_summary",
        "donor_or_sample_profiles",
        "genes_tested",
        "n_genes_tested",
        "n_FDR_0_05",
        "n_FDR_0_10",
        "multiple_testing_family",
    ]
    all_tests = all_tests[keep].copy()
    all_tests["genes_tested"] = all_tests["genes_tested"].fillna(all_tests["n_genes_tested"])
    all_tests = all_tests.drop(columns=["n_genes_tested"])
    all_tests["display_location"] = "Complete all-tests table"
    focus = (
        all_tests["dataset"].eq("GSE206528")
        & all_tests["analysis_level"].eq("refined mural/SMC subtype pseudobulk")
        & all_tests["compartment_or_system"].eq("Pericyte_like")
        & all_tests["comparison"].eq("diabetic_ED_vs_normal")
    )
    all_tests.loc[focus, "display_location"] = "Primary focus; full row retained in complete output table"
    all_tests["selection_note"] = np.where(
        focus,
        "Strongest refined diabetic ED signal; interpreted as exploratory because n=2 diabetic ED donors and tumor-margin references.",
        "Reported for transparency as tested negative, secondary, or cross-context support.",
    )
    return all_tests.sort_values(["dataset", "analysis_level", "compartment_or_system", "comparison"]).reset_index(drop=True)


def build_lodo_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    diag = read_tsv(TABLE_DIR / "GSE206528_mural_subtype_edger_lodo" / "lodo_diagnostics.tsv")
    samples = read_tsv(ROOT / "metadata" / "GSE206528" / "samples.tsv")[["donor_id", "ed_subgroup"]]
    diag = diag.merge(samples.rename(columns={"ed_subgroup": "omitted_donor_group"}), left_on="leave_out_donor", right_on="donor_id", how="left")
    diag = diag.drop(columns=["donor_id"])
    diag["recommended_wording"] = np.where(
        (diag["comparison"].eq("diabetic_ED_vs_normal")) & (diag["omitted_donor_group"].eq("normal")),
        "normal-reference-donor leave-one-out sensitivity",
        np.where(
            (diag["comparison"].eq("diabetic_ED_vs_normal")) & (diag["omitted_donor_group"].eq("diabetic_ED")),
            "not runnable: omitting a diabetic ED donor leaves one diabetic ED donor",
            "leave-one-donor sensitivity/context row",
        ),
    )

    detail = read_tsv(TABLE_DIR / "GSE206528_mural_subtype_edger_lodo" / "candidate_lodo_edgeR_results.tsv.gz")
    detail = detail[
        detail["mural_subtype"].eq("Pericyte_like")
        & detail["comparison"].eq("diabetic_ED_vs_normal")
        & detail["gene_symbol"].isin(CANDIDATE_GENES)
        & detail["leave_out_donor"].isin(["N1", "N2", "N3"])
    ].copy()
    summary = (
        detail.groupby(["gene_symbol", "full_logFC", "full_FDR"], as_index=False)
        .agg(
            n_runnable_normal_reference_omissions=("leave_out_donor", "nunique"),
            min_lodo_logFC=("lodo_logFC", "min"),
            max_lodo_logFC=("lodo_logFC", "max"),
            max_lodo_FDR=("lodo_FDR", "max"),
            same_direction_all_runnable=("same_direction", "all"),
            nominal_p_lt_0_05_all_runnable=("lodo_nominal_0_05", "all"),
            fdr_lt_0_10_all_runnable=("lodo_FDR_0_10", "all"),
        )
    )
    summary["interpretation"] = "Direction preserved when each normal-reference donor is omitted; not independent replication."
    return diag, summary.sort_values("full_FDR")


def build_reporting_checklist() -> pd.DataFrame:
    rows = [
        ("GEO accessions named", "Complete", "Methods; metadata/dataset_registry.tsv", "GSE206528, GSE126005, and GSE146500 are named."),
        ("Sample/donor mapping", "Complete", "metadata/GSE206528/samples.tsv; Methods", "Donor IDs, ED subgroup, diabetes status, and tissue-source note are recorded."),
        ("Biological replicate definition", "Complete", "Methods", "Donor-level pseudobulk profiles are the discovery biological replicates."),
        ("Raw-count usage", "Complete", "Methods; data/processed/GSE206528/pseudobulk_*", "Raw counts were retained for pseudobulk edgeR modeling."),
        ("QC/filtering", "Complete", "Methods; results/tables/*qc*.tsv", "Exploratory gene/cell filtering and minimum pseudobulk cell thresholds are reported."),
        ("Clustering/annotation parameters", "Complete", "Methods; metadata/GSE206528/*annotations.tsv", "Leiden r0.5 and marker-based major/refined mural labels are documented."),
        ("Pseudobulk aggregation rule", "Complete", "Methods", "One profile per donor x major cell type or donor x refined mural subtype."),
        ("Minimum-cell threshold", "Complete", "Methods", f"Formal testing required at least {MIN_CELLS_FORMAL} cells per donor-compartment profile."),
        ("Statistical model", "Complete", "Methods", "edgeR quasi-likelihood negative-binomial models with TMM normalization and robust empirical Bayes moderation."),
        ("Covariates not included", "Complete", "Methods; Discussion", "No covariates were added because donor number and clinical metadata were limited."),
        ("Multiple-testing family", "Complete", "Complete all-tests table", "BH FDR family is listed for every edgeR result table."),
        ("Leave-one-out runnability", "Complete", "Supplementary normal-reference leave-one-out table", "Diabetic-donor omission is explicitly marked as not runnable for the key contrast."),
        ("Permutation sanity check", "Complete", "Supplementary exact permutation table", "All ten possible 2-vs-3 labels among the five donor profiles are listed."),
        ("Diabetic vs non-diabetic descriptive check", "Complete", "Supplementary descriptive 12-gene score table", "Normal, non-diabetic ED, and diabetic ED pericyte-like profiles are shown; ED2 is flagged below threshold."),
    ]
    return pd.DataFrame(rows, columns=["item", "status", "location", "note"])


def plot_permutation(permutation: pd.DataFrame) -> None:
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    plot_df = permutation.sort_values("score_difference", ascending=True)
    colors = np.where(plot_df["is_observed_diabetic_assignment"], "#B23A48", "#778899")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.barh(plot_df["putative_diabetic_pair"], plot_df["score_difference"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Mean score difference: putative diabetic pair minus remaining triplet")
    ax.set_ylabel("Two-donor assignment")
    ax.set_title("Exact 2-vs-3 label-permutation sanity check")
    ax.text(
        0.02,
        0.04,
        "Observed DMED1,DMED2 assignment highlighted",
        transform=ax.transAxes,
        fontsize=8,
        color="#333333",
    )
    fig.tight_layout()
    fig.savefig(OUT_FIG / "exact_label_permutation_12gene_directional_score.png", dpi=300)
    fig.savefig(OUT_FIG / "exact_label_permutation_12gene_directional_score.pdf")
    plt.close(fig)


def plot_descriptive_scores(scores: pd.DataFrame) -> None:
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    sub = scores[scores["module_id"].eq("candidate_directional_signature_12")].copy()
    order = ["normal", "non_diabetic_ED", "diabetic_ED"]
    labels = ["Tumor-margin\nreference", "Non-diabetic\norganic ED", "Diabetic ED"]
    colors = {"normal": "#4C78A8", "non_diabetic_ED": "#59A14F", "diabetic_ED": "#B23A48"}
    xmap = {group: i for i, group in enumerate(order)}
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    for group in order:
        g = sub[sub["ed_subgroup"].eq(group)].copy()
        xs = np.full(len(g), xmap[group], dtype=float)
        jitter = np.linspace(-0.07, 0.07, len(g)) if len(g) > 1 else np.array([0])
        alpha = np.where(g["passes_primary_min_cell_threshold"], 1.0, 0.45)
        for i, (_, row) in enumerate(g.iterrows()):
            marker = "o" if row["passes_primary_min_cell_threshold"] else "s"
            ax.scatter(
                xs[i] + jitter[i],
                row["mean_directional_reference_z"],
                s=55,
                color=colors[group],
                alpha=float(alpha[i]),
                marker=marker,
                edgecolor="black",
                linewidth=0.5,
                zorder=3,
            )
            ax.text(xs[i] + jitter[i], row["mean_directional_reference_z"] + 0.09, row["donor_id"], ha="center", fontsize=7)
        if len(g):
            ax.hlines(g["mean_directional_reference_z"].mean(), xmap[group] - 0.22, xmap[group] + 0.22, colors="black", linewidth=1.2)
    ax.axhline(0, color="#333333", linewidth=0.8, linestyle="--")
    ax.set_xticks(range(len(order)), labels)
    ax.set_ylabel("12-gene directional score\n(normal-reference z scale)")
    ax.set_title("Descriptive diabetic vs non-diabetic ED check")
    ax.text(
        0.02,
        0.04,
        "Square marks profile below the formal 50-cell testing threshold",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#333333",
    )
    fig.tight_layout()
    fig.savefig(OUT_FIG / "diabetic_vs_nondiabetic_descriptive_12gene_score.png", dpi=300)
    fig.savefig(OUT_FIG / "diabetic_vs_nondiabetic_descriptive_12gene_score.pdf")
    plt.close(fig)


def main() -> None:
    OUT_TABLE.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)

    permutation, permutation_summary = build_exact_permutation()
    write_table(permutation, "exact_label_permutation_12gene_directional_score")
    write_table(permutation_summary, "exact_label_permutation_12gene_directional_score_summary")
    plot_permutation(permutation)

    gene_level, scores, score_summary = build_descriptive_scores()
    write_table(gene_level, "diabetic_vs_nondiabetic_12gene_gene_level_reference_z")
    write_table(scores, "diabetic_vs_nondiabetic_12gene_module_scores")
    write_table(score_summary, "diabetic_vs_nondiabetic_12gene_module_score_group_summary")
    plot_descriptive_scores(scores)

    write_table(build_all_tests_table(), "all_edgeR_tests_complete")
    module_tests = read_tsv(TABLE_DIR / "robustness_checks" / "GSE206528_refined_mural_module_cameraPR_results.tsv")
    write_table(module_tests, "all_GSE206528_refined_mural_cameraPR_module_tests")

    lodo_diag, lodo_candidate_summary = build_lodo_tables()
    write_table(lodo_diag, "normal_reference_donor_leave_one_out_runnability")
    write_table(lodo_candidate_summary, "normal_reference_donor_leave_one_out_candidate_summary")

    write_table(build_reporting_checklist(), "secondary_scrnaseq_reanalysis_reporting_checklist")

    print(f"Wrote reproducibility-check tables to {OUT_TABLE}")
    print(f"Wrote reproducibility-check figures to {OUT_FIG}")


if __name__ == "__main__":
    main()
