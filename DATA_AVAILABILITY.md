# Data Availability

All primary source datasets are public or otherwise identified below.

## Public GEO Datasets

- `GSE206528`: human corpus cavernosum single-cell RNA-seq discovery dataset.
- `GSE259348`: independent single-donor diabetic erectile dysfunction single-cell RNA-seq dataset.
- `GSE259299`: rat corpus cavernosum single-cell RNA-seq boundary dataset.
- `GSE228900`: human perivascular in-vitro diabetic-stimulus model.
- `GSE126005`: Peyronie's plaque versus control tunica albuginea bulk RNA-seq.
- `GSE146500`: Peyronie's fibroblast bulk RNA-seq.
- `GSE261085`: normal human corpus cavernosum spatial transcriptomics boundary dataset.
- `GSE131882`: diabetic nephropathy kidney snRNA-seq cross-tissue boundary dataset.

Dataset URLs and roles are listed in `metadata/dataset_registry.tsv`.

## Controlled-Access Dataset

The aging-related erectile dysfunction atlas is available through GSA-Human accession `HRA010374`. This repository does not redistribute or reprocess controlled-access raw data. The source article reports the accession in its data-availability statement and provides open supplementary material at `https://www.frontiersin.org/articles/10.3389/fendo.2025.1671482/full#supplementary-material`. This repository uses that open supplement only as a non-rerunnable external context check for whether the nominated genes appear among reported aging-related erectile dysfunction pericyte markers.

## Included Data Products

This repository includes compact processed products needed to inspect the reported analyses:

- donor-level pseudobulk matrices for `GSE206528` major cell types and refined mural subtypes;
- processed count tables for `GSE126005` and `GSE146500`;
- sample metadata and group mappings;
- marker and annotation tables;
- full edgeR result tables;
- reproducibility, robustness, sensitivity, and cross-context result tables;
- selected generated figure outputs.

## Excluded Large Files

Raw GEO archives and full `.h5ad` objects are excluded from this repository. They can be regenerated from source data with `scripts/fetch_external_datasets.sh` and the preprocessing scripts in `scripts/`. If a complete raw-to-results archival package is required, those large files should be deposited in a separate large-file archive and cited alongside this repository.
