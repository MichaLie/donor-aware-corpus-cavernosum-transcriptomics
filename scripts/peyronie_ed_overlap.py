#!/usr/bin/env python3
"""Direction-aware overlap between robust ED genes and Peyronie's validation contrasts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]


def read_table(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    return pd.read_csv(path, sep="\t", compression=compression)


def sign(value: float) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def load_ed() -> pd.DataFrame:
    ed = read_table(
        ROOT
        / "results/tables/GSE206528_pseudobulk_edger_lodo/robust_FDR_0_05_candidates_with_programs.tsv"
    )
    ed = ed.rename(
        columns={
            "comparison": "ed_comparison",
            "full_logFC": "ed_logFC",
            "full_FDR": "ed_FDR",
        }
    )
    ed["gene_symbol"] = ed["gene_symbol"].astype(str)
    ed["ed_sign"] = ed["ed_logFC"].map(sign)
    return ed


def load_peyronie() -> pd.DataFrame:
    paths = [
        ROOT / "results/tables/GSE126005/GSE126005_edgeR_all_results.tsv.gz",
        ROOT / "results/tables/GSE146500/GSE146500_edgeR_all_results.tsv.gz",
    ]
    frames = [read_table(path) for path in paths]
    pey = pd.concat(frames, ignore_index=True)
    pey["gene_symbol"] = pey["gene_symbol"].astype(str)
    pey["peyronie_sign"] = pey["logFC"].map(sign)
    pey["peyronie_FDR_0_05"] = pey["FDR"] < 0.05
    pey["peyronie_FDR_0_10"] = pey["FDR"] < 0.10
    return pey


def summarize_overlap(ed: pd.DataFrame, pey: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    details = []
    summary = []

    pey_groups = pey.groupby(["dataset", "comparison"], sort=True)
    ed_groups = ed.groupby(["cell_type", "ed_comparison"], sort=True)

    for (cell_type, ed_comparison), ed_group in ed_groups:
        ed_set = set(ed_group["gene_symbol"])
        for (dataset, pey_comparison), pey_group in pey_groups:
            pey_universe = set(pey_group["gene_symbol"])
            pey_sig = set(pey_group.loc[pey_group["peyronie_FDR_0_05"], "gene_symbol"])
            ed_in_universe = ed_set & pey_universe
            overlap_sig = ed_in_universe & pey_sig

            joined = ed_group.merge(
                pey_group,
                on="gene_symbol",
                how="inner",
                suffixes=("", "_peyronie"),
            )
            joined_sig = joined[joined["peyronie_FDR_0_05"]].copy()
            if not joined_sig.empty:
                joined_sig["same_direction"] = joined_sig["ed_sign"] == joined_sig["peyronie_sign"]
                joined_sig["opposite_direction"] = joined_sig["ed_sign"] == -joined_sig["peyronie_sign"]
                joined_sig["validation_mode"] = joined_sig.apply(
                    lambda row: classify_validation(row["comparison"], row["same_direction"], row["opposite_direction"]),
                    axis=1,
                )
                details.append(
                    joined_sig[
                        [
                            "cell_type",
                            "ed_comparison",
                            "dataset",
                            "comparison",
                            "gene_symbol",
                            "ed_logFC",
                            "ed_FDR",
                            "logFC",
                            "FDR",
                            "same_direction",
                            "opposite_direction",
                            "validation_mode",
                            "program_ids",
                            "has_curated_program",
                        ]
                    ].rename(columns={"comparison": "peyronie_comparison", "logFC": "peyronie_logFC", "FDR": "peyronie_FDR"})
                )

            a = len(overlap_sig)
            b = len(ed_in_universe) - a
            c = len(pey_sig - ed_in_universe)
            d = len(pey_universe) - a - b - c
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

            n_same = int((joined_sig["ed_sign"] == joined_sig["peyronie_sign"]).sum()) if not joined_sig.empty else 0
            n_opp = int((joined_sig["ed_sign"] == -joined_sig["peyronie_sign"]).sum()) if not joined_sig.empty else 0

            summary.append(
                {
                    "cell_type": cell_type,
                    "ed_comparison": ed_comparison,
                    "dataset": dataset,
                    "peyronie_comparison": pey_comparison,
                    "n_ed_genes": len(ed_set),
                    "n_ed_genes_in_peyronie_universe": len(ed_in_universe),
                    "n_peyronie_FDR_0_05": len(pey_sig),
                    "n_overlap_peyronie_FDR_0_05": a,
                    "n_same_direction": n_same,
                    "n_opposite_direction": n_opp,
                    "fisher_odds_ratio": odds_ratio,
                    "fisher_p_value": p_value,
                }
            )

    detail_df = pd.concat(details, ignore_index=True) if details else pd.DataFrame()
    summary_df = pd.DataFrame(summary).sort_values(
        ["n_overlap_peyronie_FDR_0_05", "n_same_direction", "cell_type"], ascending=[False, False, True]
    )
    summary_df["fisher_FDR"] = benjamini_hochberg(summary_df["fisher_p_value"])
    return summary_df, detail_df


def classify_validation(peyronie_comparison: str, same_direction: bool, opposite_direction: bool) -> str:
    if peyronie_comparison == "PFPC_vs_PF":
        if opposite_direction:
            return "treatment_opposes_ED_direction"
        if same_direction:
            return "treatment_matches_ED_direction"
        return "treatment_flat_or_ambiguous"
    if same_direction:
        return "disease_same_direction"
    if opposite_direction:
        return "disease_opposite_direction"
    return "disease_flat_or_ambiguous"


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    p = p_values.astype(float).reset_index(drop=True)
    order = p.sort_values().index
    ranked = p.loc[order]
    m = len(p)
    adjusted = ranked * m / pd.Series(range(1, m + 1), index=ranked.index)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    out = pd.Series(index=p.index, dtype=float)
    out.loc[order] = adjusted
    return out.values


def program_summary(details: pd.DataFrame) -> pd.DataFrame:
    if details.empty:
        return pd.DataFrame()
    rows = []
    exploded = details.dropna(subset=["program_ids"]).copy()
    exploded["program_id"] = exploded["program_ids"].str.split(r"[;,]")
    exploded = exploded.explode("program_id")
    exploded["program_id"] = exploded["program_id"].str.strip()
    exploded = exploded[exploded["program_id"] != ""]
    for keys, group in exploded.groupby(["dataset", "peyronie_comparison", "program_id"], sort=True):
        dataset, comparison, program_id = keys
        rows.append(
            {
                "dataset": dataset,
                "peyronie_comparison": comparison,
                "program_id": program_id,
                "n_overlap_rows": len(group),
                "n_unique_genes": group["gene_symbol"].nunique(),
                "genes": ",".join(sorted(group["gene_symbol"].unique())),
                "n_same_direction": int(group["same_direction"].sum()),
                "n_opposite_direction": int(group["opposite_direction"].sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["n_unique_genes", "dataset"], ascending=[False, True])


def main() -> None:
    out_dir = ROOT / "results/tables/peyronie_validation"
    out_dir.mkdir(parents=True, exist_ok=True)

    ed = load_ed()
    pey = load_peyronie()
    summary, details = summarize_overlap(ed, pey)
    programs = program_summary(details)

    summary.to_csv(out_dir / "ed_robust_gene_overlap_summary.tsv", sep="\t", index=False)
    details.to_csv(out_dir / "ed_robust_gene_overlap_details.tsv", sep="\t", index=False)
    programs.to_csv(out_dir / "ed_robust_program_overlap_summary.tsv", sep="\t", index=False)

    print(f"Wrote {len(summary)} overlap summary rows")
    print(f"Wrote {len(details)} overlap detail rows")
    print(f"Wrote {len(programs)} program summary rows")


if __name__ == "__main__":
    main()
