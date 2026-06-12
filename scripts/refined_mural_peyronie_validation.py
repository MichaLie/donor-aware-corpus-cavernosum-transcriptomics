#!/usr/bin/env python3
"""Direction-aware Peyronie's overlap for refined GSE206528 mural subtype signatures."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip" if path.suffix == ".gz" else None)


def sign(value: float) -> int:
    return 1 if value > 0 else -1 if value < 0 else 0


def classify(peyronie_comparison: str, same: bool, opposite: bool) -> str:
    if peyronie_comparison == "PFPC_vs_PF":
        if opposite:
            return "treatment_opposes_ED_direction"
        if same:
            return "treatment_matches_ED_direction"
        return "treatment_ambiguous"
    if same:
        return "disease_same_direction"
    if opposite:
        return "disease_opposite_direction"
    return "disease_ambiguous"


def bh(p_values: pd.Series) -> pd.Series:
    p = p_values.astype(float).reset_index(drop=True)
    order = p.sort_values().index
    ranked = p.loc[order]
    m = len(p)
    adjusted = ranked * m / pd.Series(range(1, m + 1), index=ranked.index)
    adjusted = adjusted.iloc[::-1].cummin().iloc[::-1].clip(upper=1.0)
    out = pd.Series(index=p.index, dtype=float)
    out.loc[order] = adjusted
    return out.values


def main() -> None:
    out_dir = ROOT / "results/tables/refined_mural_validation"
    out_dir.mkdir(parents=True, exist_ok=True)

    ed = read_tsv(ROOT / "results/tables/GSE206528_mural_subtype_edger/all_mural_subtype_edgeR_results.tsv.gz")
    ed = ed.loc[ed["FDR"] < 0.05].copy()
    ed["ed_sign"] = ed["logFC"].map(sign)
    ed["direction"] = ed["logFC"].map(lambda x: "up_in_ED" if x > 0 else "down_in_ED")

    pey = pd.concat(
        [
            read_tsv(ROOT / "results/tables/GSE126005/GSE126005_edgeR_all_results.tsv.gz"),
            read_tsv(ROOT / "results/tables/GSE146500/GSE146500_edgeR_all_results.tsv.gz"),
        ],
        ignore_index=True,
    )
    pey = pey.dropna(subset=["gene_symbol"]).copy()
    pey["peyronie_sign"] = pey["logFC"].map(sign)
    pey["peyronie_FDR_0_05"] = pey["FDR"] < 0.05

    summaries: list[dict[str, object]] = []
    detail_frames: list[pd.DataFrame] = []

    for (mural_subtype, ed_comparison), ed_group in ed.groupby(["mural_subtype", "comparison"]):
        ed_set = set(ed_group["gene_symbol"])
        for (dataset, pey_comp), pey_group in pey.groupby(["dataset", "comparison"]):
            universe = set(pey_group["gene_symbol"])
            pey_sig = set(pey_group.loc[pey_group["peyronie_FDR_0_05"], "gene_symbol"])
            ed_in_universe = ed_set & universe
            overlap = ed_in_universe & pey_sig

            joined = ed_group.merge(
                pey_group.loc[pey_group["peyronie_FDR_0_05"]],
                on="gene_symbol",
                suffixes=("_ed", "_peyronie"),
            )
            if not joined.empty:
                joined["same_direction"] = joined["ed_sign"] == joined["peyronie_sign"]
                joined["opposite_direction"] = joined["ed_sign"] == -joined["peyronie_sign"]
                joined["validation_mode"] = joined.apply(
                    lambda row: classify(pey_comp, row["same_direction"], row["opposite_direction"]), axis=1
                )
                detail_frames.append(
                    joined[
                        [
                            "mural_subtype",
                            "comparison_ed",
                            "dataset",
                            "comparison_peyronie",
                            "gene_symbol",
                            "logFC_ed",
                            "FDR_ed",
                            "logFC_peyronie",
                            "FDR_peyronie",
                            "same_direction",
                            "opposite_direction",
                            "validation_mode",
                        ]
                    ].rename(
                        columns={
                            "comparison_ed": "ed_comparison",
                            "comparison_peyronie": "peyronie_comparison",
                            "logFC_ed": "ed_logFC",
                            "FDR_ed": "ed_FDR",
                            "logFC_peyronie": "peyronie_logFC",
                            "FDR_peyronie": "peyronie_FDR",
                        }
                    )
                )

            a = len(overlap)
            b = len(ed_in_universe) - a
            c = len(pey_sig - ed_in_universe)
            d = len(universe) - a - b - c
            odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
            summaries.append(
                {
                    "mural_subtype": mural_subtype,
                    "ed_comparison": ed_comparison,
                    "dataset": dataset,
                    "peyronie_comparison": pey_comp,
                    "n_ed_FDR_0_05_genes": len(ed_set),
                    "n_ed_genes_in_peyronie_universe": len(ed_in_universe),
                    "n_peyronie_FDR_0_05": len(pey_sig),
                    "n_overlap_peyronie_FDR_0_05": a,
                    "n_same_direction": int((joined["ed_sign"] == joined["peyronie_sign"]).sum()) if not joined.empty else 0,
                    "n_opposite_direction": int((joined["ed_sign"] == -joined["peyronie_sign"]).sum()) if not joined.empty else 0,
                    "fisher_odds_ratio": odds,
                    "fisher_p_value": p_value,
                }
            )

    summary = pd.DataFrame(summaries).sort_values(
        ["n_overlap_peyronie_FDR_0_05", "n_same_direction"], ascending=[False, False]
    )
    summary["fisher_FDR"] = bh(summary["fisher_p_value"])
    details = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()

    summary.to_csv(out_dir / "refined_mural_peyronie_overlap_summary.tsv", sep="\t", index=False)
    details.to_csv(out_dir / "refined_mural_peyronie_overlap_details.tsv", sep="\t", index=False)
    print(f"Wrote {len(summary)} refined overlap summary rows")
    print(f"Wrote {len(details)} refined overlap detail rows")


if __name__ == "__main__":
    main()
