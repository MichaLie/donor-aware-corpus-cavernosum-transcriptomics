#!/usr/bin/env python3
"""Annotate robust GSE206528 pseudobulk candidates with curated disease programs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
GENE_SET_FILE = ROOT / "metadata" / "gene_sets" / "curated_disease_programs.tsv"
LODO_DIR = ROOT / "results" / "tables" / "GSE206528_pseudobulk_edger_lodo"
LODO_SUMMARY = LODO_DIR / "candidate_lodo_summary.tsv"


def main() -> None:
    gene_sets = pd.read_csv(GENE_SET_FILE, sep="\t")
    gene_sets["gene_set"] = gene_sets["gene_symbols"].str.split(",").map(set)

    candidates = pd.read_csv(LODO_SUMMARY, sep="\t")
    robust = candidates.loc[
        (candidates["full_FDR"] < 0.05)
        & (candidates["stable_tested_direction_all_lodo"])
    ].copy()
    robust["direction"] = robust["full_logFC"].map(lambda x: "up_in_ED" if x > 0 else "down_in_ED")

    annotation_rows: list[dict[str, object]] = []
    for row in robust.itertuples(index=False):
        for program in gene_sets.itertuples(index=False):
            if row.gene_symbol in program.gene_set:
                annotation_rows.append(
                    {
                        "cell_type": row.cell_type,
                        "comparison": row.comparison,
                        "gene_symbol": row.gene_symbol,
                        "direction": row.direction,
                        "full_logFC": row.full_logFC,
                        "full_FDR": row.full_FDR,
                        "program_id": program.program_id,
                        "program_label": program.program_label,
                        "program_family": program.program_family,
                        "source_status": program.source_status,
                    }
                )

    annotated = robust.copy()
    program_lookup = {
        gene: sorted(
            gene_sets.loc[
                gene_sets["gene_set"].map(lambda genes: gene in genes),
                "program_id",
            ].tolist()
        )
        for gene in robust["gene_symbol"]
    }
    annotated["program_ids"] = annotated["gene_symbol"].map(
        lambda gene: ",".join(program_lookup.get(gene, []))
    )
    annotated["has_curated_program"] = annotated["program_ids"] != ""
    annotated.to_csv(
        LODO_DIR / "robust_FDR_0_05_candidates_with_programs.tsv",
        sep="\t",
        index=False,
    )

    overlaps = pd.DataFrame(annotation_rows)
    if overlaps.empty:
      overlaps = pd.DataFrame(
          columns=[
              "cell_type",
              "comparison",
              "gene_symbol",
              "direction",
              "full_logFC",
              "full_FDR",
              "program_id",
              "program_label",
              "program_family",
              "source_status",
          ]
      )
    overlaps.to_csv(
        LODO_DIR / "robust_FDR_0_05_program_gene_overlaps.tsv",
        sep="\t",
        index=False,
    )

    if overlaps.empty:
        overview = pd.DataFrame()
    else:
        overview = (
            overlaps.groupby(
                ["cell_type", "comparison", "program_id", "program_label", "program_family"],
                as_index=False,
            )
            .agg(
                n_genes=("gene_symbol", "nunique"),
                up_genes=("gene_symbol", lambda x: ",".join(sorted(set(overlaps.loc[x.index][overlaps.loc[x.index, "direction"] == "up_in_ED"]["gene_symbol"])))),
                down_genes=("gene_symbol", lambda x: ",".join(sorted(set(overlaps.loc[x.index][overlaps.loc[x.index, "direction"] == "down_in_ED"]["gene_symbol"])))),
            )
            .sort_values(["cell_type", "comparison", "n_genes"], ascending=[True, True, False])
        )
    overview.to_csv(
        LODO_DIR / "robust_FDR_0_05_program_overlap_overview.tsv",
        sep="\t",
        index=False,
    )

    print("Robust candidates:", robust.shape[0])
    print("Robust candidates with >=1 curated program:", int(annotated["has_curated_program"].sum()))
    print("Wrote program overlap outputs to:", LODO_DIR)


if __name__ == "__main__":
    main()
