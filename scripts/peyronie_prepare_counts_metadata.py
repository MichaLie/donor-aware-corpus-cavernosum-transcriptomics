#!/usr/bin/env python3
"""Prepare Peyronie's validation count matrices and sample metadata.

Inputs are GEO processed/raw count supplements for:
- GSE126005: bulk tunica albuginea, Peyronie's plaque versus control.
- GSE146500: primary tunica albuginea fibroblasts, normal/Peyronie's/PC-NV-treated.
"""

from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


GSE126005_SAMPLES = [
    ("GSM3587848", "Seq4@S1", "peyronies_plaque", "Peyronie's plaque"),
    ("GSM3587849", "Seq8@1@S2", "peyronies_plaque", "Peyronie's plaque"),
    ("GSM3587850", "s1@S3", "peyronies_plaque", "Peyronie's plaque"),
    ("GSM3587851", "s7@S4", "peyronies_plaque", "Peyronie's plaque"),
    ("GSM3587852", "s16@S7", "peyronies_plaque", "Peyronie's plaque"),
    ("GSM3587853", "Seq27@S12", "peyronies_plaque", "Peyronie's plaque"),
    ("GSM3587854", "s12@S5", "normal_control", "Control tunica albuginea"),
    ("GSM3587855", "s14@S6", "normal_control", "Control tunica albuginea"),
    ("GSM3587856", "Seq21@2@S8", "normal_control", "Control tunica albuginea"),
    ("GSM3587857", "Seq22@2@S9", "normal_control", "Control tunica albuginea"),
    ("GSM3587858", "Seq24@2@S10", "normal_control", "Control tunica albuginea"),
    ("GSM3587859", "Seq26@2@S11", "normal_control", "Control tunica albuginea"),
]


GSE146500_SAMPLES = [
    ("GSM4387232", "NF1", "Con-1", "normal_fibroblast", "Normal TA fibroblast"),
    ("GSM4387233", "NF2", "Con-2", "normal_fibroblast", "Normal TA fibroblast"),
    ("GSM4387234", "NF3", "Con-3", "normal_fibroblast", "Normal TA fibroblast"),
    ("GSM4387235", "NF4", "Con-4", "normal_fibroblast", "Normal TA fibroblast"),
    ("GSM4387236", "PF1", "PT-1", "peyronies_fibroblast", "Peyronie's TA fibroblast"),
    ("GSM4387237", "PF2", "PT-2", "peyronies_fibroblast", "Peyronie's TA fibroblast"),
    ("GSM4387238", "PF3", "PT-3", "peyronies_fibroblast", "Peyronie's TA fibroblast"),
    ("GSM4387239", "PF4", "PT-4", "peyronies_fibroblast", "Peyronie's TA fibroblast"),
    (
        "GSM4387240",
        "PFPC1",
        "PTPC-1",
        "peyronies_fibroblast_PC_NV",
        "Peyronie's TA fibroblast treated with PC-NVs",
    ),
    (
        "GSM4387241",
        "PFPC2",
        "PTPC-2",
        "peyronies_fibroblast_PC_NV",
        "Peyronie's TA fibroblast treated with PC-NVs",
    ),
    (
        "GSM4387242",
        "PFPC3",
        "PTPC-3",
        "peyronies_fibroblast_PC_NV",
        "Peyronie's TA fibroblast treated with PC-NVs",
    ),
    (
        "GSM4387243",
        "PFPC4",
        "PTPC-4",
        "peyronies_fibroblast_PC_NV",
        "Peyronie's TA fibroblast treated with PC-NVs",
    ),
]


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip")


def sample_qc(counts: pd.DataFrame, metadata: pd.DataFrame, feature_col: str) -> pd.DataFrame:
    sample_cols = metadata["sample_id"].tolist()
    matrix = counts[sample_cols]
    qc = pd.DataFrame(
        {
            "sample_id": sample_cols,
            "library_size": matrix.sum(axis=0).astype(int).values,
            "detected_genes": (matrix > 0).sum(axis=0).astype(int).values,
        }
    )
    qc = qc.merge(metadata, on="sample_id", how="left")
    qc["n_features_total"] = counts[feature_col].nunique()
    return qc


def prepare_gse126005() -> None:
    metadata = pd.DataFrame(
        GSE126005_SAMPLES,
        columns=["sample_id", "sample_title", "group", "group_label"],
    )
    metadata["dataset"] = "GSE126005"
    metadata["tissue"] = "tunica_albuginea"
    metadata["modality"] = "bulk_RNA_seq"
    metadata["raw_count_column"] = metadata["sample_title"]
    write_tsv(metadata, ROOT / "metadata/GSE126005/samples.tsv")

    raw_path = ROOT / "data/raw/GSE126005/GSE126005_RawCountsTable.txt.gz"
    raw = pd.read_csv(raw_path, sep="\t", compression="gzip")
    title_to_sample = dict(zip(metadata["sample_title"], metadata["sample_id"]))
    sample_titles = metadata["sample_title"].tolist()
    missing = sorted(set(sample_titles) - set(raw.columns))
    if missing:
        raise ValueError(f"GSE126005 raw count table missing sample columns: {missing}")

    counts = raw[["Gene ID", "Gene Name", *sample_titles]].copy()
    counts = counts.rename(columns={"Gene ID": "gene_id", "Gene Name": "gene_symbol", **title_to_sample})
    counts["gene_symbol"] = counts["gene_symbol"].astype(str)
    write_tsv_gz(counts, ROOT / "data/processed/GSE126005/GSE126005_counts_gene.tsv.gz")

    feature_cols = [
        "Gene ID",
        "Gene Name",
        "Chromosome",
        "StartPosition",
        "EndPosition",
        "Strand",
        "Length",
        "RNAType",
        "Description",
    ]
    feature_cols = [col for col in feature_cols if col in raw.columns]
    features = raw[feature_cols].rename(columns={"Gene ID": "gene_id", "Gene Name": "gene_symbol"})
    write_tsv_gz(features, ROOT / "data/processed/GSE126005/GSE126005_gene_annotation.tsv.gz")

    qc = sample_qc(counts, metadata, "gene_id")
    write_tsv(qc, ROOT / "results/tables/GSE126005/GSE126005_bulk_sample_qc.tsv")


def prepare_gse146500() -> None:
    metadata = pd.DataFrame(
        GSE146500_SAMPLES,
        columns=["sample_id", "sample_title", "file_label", "group", "group_label"],
    )
    metadata["dataset"] = "GSE146500"
    metadata["tissue"] = "tunica_albuginea_fibroblast"
    metadata["modality"] = "bulk_RNA_seq"
    metadata["agent"] = metadata["group"].map(
        {
            "normal_fibroblast": "none",
            "peyronies_fibroblast": "none",
            "peyronies_fibroblast_PC_NV": "PC-NV",
        }
    )
    write_tsv(metadata, ROOT / "metadata/GSE146500/samples.tsv")

    tar_path = ROOT / "data/raw/GSE146500/GSE146500_RAW.tar"
    sample_tables = []
    with tarfile.open(tar_path, "r") as tar:
        members = {Path(member.name).name: member for member in tar.getmembers()}
        for row in metadata.itertuples(index=False):
            suffix = f"{row.file_label}.txt.gz"
            matches = [name for name in members if name.endswith(suffix)]
            if len(matches) != 1:
                raise ValueError(f"Expected one file ending {suffix}, found {matches}")
            handle = tar.extractfile(members[matches[0]])
            if handle is None:
                raise ValueError(f"Could not extract {matches[0]}")
            with gzip.GzipFile(fileobj=handle) as gz:
                table = pd.read_csv(gz, sep="\t")
            if "Gene symbol" not in table.columns:
                raise ValueError(f"{matches[0]} does not contain a Gene symbol column")
            value_cols = [col for col in table.columns if col != "Gene symbol"]
            if len(value_cols) != 1:
                raise ValueError(f"{matches[0]} has unexpected count columns: {value_cols}")
            table = table.rename(columns={"Gene symbol": "gene_symbol", value_cols[0]: row.sample_id})
            table = table.groupby("gene_symbol", as_index=False)[row.sample_id].sum()
            sample_tables.append(table)

    counts = sample_tables[0]
    for table in sample_tables[1:]:
        counts = counts.merge(table, on="gene_symbol", how="outer")
    sample_ids = metadata["sample_id"].tolist()
    counts[sample_ids] = counts[sample_ids].fillna(0).astype(int)
    counts = counts[["gene_symbol", *sample_ids]].sort_values("gene_symbol")
    write_tsv_gz(counts, ROOT / "data/processed/GSE146500/GSE146500_counts_gene_symbol.tsv.gz")

    qc = sample_qc(counts, metadata, "gene_symbol")
    write_tsv(qc, ROOT / "results/tables/GSE146500/GSE146500_bulk_sample_qc.tsv")


def main() -> None:
    prepare_gse126005()
    prepare_gse146500()


if __name__ == "__main__":
    main()
