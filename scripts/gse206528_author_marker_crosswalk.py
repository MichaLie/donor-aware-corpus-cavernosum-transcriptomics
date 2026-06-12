#!/usr/bin/env python3
"""Compare Leiden clusters with author-reported GSE206528 marker tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP_DIR = ROOT / "metadata" / "GSE206528" / "publication_supplement"
RESULTS_DIR = ROOT / "results" / "tables"

OUR_MARKERS = RESULTS_DIR / "GSE206528_leiden_r0_5_marker_genes.tsv"
AUTHOR_MAJOR = SUPP_DIR / "41467_2022_31950_MOESM5_ESM.csv"
AUTHOR_SUBCLUSTER = SUPP_DIR / "41467_2022_31950_MOESM6_ESM.xlsx"


def top_marker_sets(df: pd.DataFrame, cluster_col: str, gene_col: str, n: int) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    for cluster, sub in df.groupby(cluster_col):
        sort_cols = [col for col in ["avg_log2FC", "scores"] if col in sub.columns]
        if sort_cols:
            sub = sub.sort_values(sort_cols[0], ascending=False)
        sets[str(cluster)] = set(sub[gene_col].dropna().astype(str).head(n))
    return sets


def make_crosswalk(author_sets: dict[str, set[str]], our_sets: dict[str, set[str]], label: str) -> pd.DataFrame:
    rows = []
    for our_cluster, our_genes in our_sets.items():
        for author_cluster, author_genes in author_sets.items():
            overlap = sorted(our_genes & author_genes)
            union = our_genes | author_genes
            rows.append(
                {
                    "comparison_level": label,
                    "leiden_r0_5": our_cluster,
                    "author_cluster": author_cluster,
                    "n_overlap_top100": len(overlap),
                    "jaccard_top100": len(overlap) / len(union) if union else 0,
                    "overlap_genes_top30_display": ",".join(overlap[:30]),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["leiden_r0_5", "n_overlap_top100", "jaccard_top100"],
        ascending=[True, False, False],
    )


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    our = pd.read_csv(OUR_MARKERS, sep="\t")
    our["group"] = our["group"].astype(str)
    our_sets = top_marker_sets(our, "group", "names", 100)

    major = pd.read_csv(AUTHOR_MAJOR)
    major_sets = top_marker_sets(major, "cluster", "gene", 100)

    subclusters = pd.read_excel(AUTHOR_SUBCLUSTER)
    subcluster_sets = top_marker_sets(subclusters, "cluster", "gene", 100)

    crosswalk = pd.concat(
        [
            make_crosswalk(major_sets, our_sets, "author_7_major_clusters"),
            make_crosswalk(subcluster_sets, our_sets, "author_11_subclusters"),
        ],
        ignore_index=True,
    )
    crosswalk.to_csv(
        RESULTS_DIR / "GSE206528_author_marker_crosswalk.tsv",
        sep="\t",
        index=False,
    )

    best = (
        crosswalk.sort_values(
            ["comparison_level", "leiden_r0_5", "n_overlap_top100", "jaccard_top100"],
            ascending=[True, True, False, False],
        )
        .groupby(["comparison_level", "leiden_r0_5"], as_index=False)
        .head(3)
    )
    best.to_csv(
        RESULTS_DIR / "GSE206528_author_marker_crosswalk_top3.tsv",
        sep="\t",
        index=False,
    )
    print("Wrote author marker crosswalk tables to:", RESULTS_DIR)


if __name__ == "__main__":
    main()
