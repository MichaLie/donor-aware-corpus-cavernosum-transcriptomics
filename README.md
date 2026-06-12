# Donor-Aware Corpus Cavernosum Transcriptomics Reanalysis

This repository contains the code, compact processed inputs, metadata, result tables, and selected figures needed to inspect and rerun a donor-aware reanalysis of public corpus cavernosum transcriptomic datasets in erectile dysfunction and diabetes.

## Scope

The primary analysis reuses public human corpus cavernosum single-cell RNA-seq data from `GSE206528`. Additional public datasets are used for directional replication, boundary checks, cross-context comparisons, and sensitivity analyses:

- `GSE206528`: human corpus cavernosum single-cell RNA-seq discovery dataset.
- `GSE259348`: independent single-donor diabetic erectile dysfunction single-cell RNA-seq dataset.
- `GSE259299`: rat corpus cavernosum single-cell RNA-seq boundary dataset.
- `GSE228900`: human perivascular in-vitro diabetic-stimulus model.
- `GSE126005` and `GSE146500`: Peyronie's disease bulk RNA-seq cross-context datasets.
- `GSE261085`: normal human corpus cavernosum spatial transcriptomics boundary dataset.
- `GSE131882`: diabetic nephropathy kidney snRNA-seq cross-tissue boundary dataset.
- `HRA010374`: controlled-access aging-related erectile dysfunction atlas; only open supplementary marker tables were used.

The repository excludes raw GEO archives and full AnnData objects to keep the archive lightweight. Raw data can be re-fetched from source repositories using `scripts/fetch_external_datasets.sh`, then processed with the scripts in `scripts/` and `scripts/repro/`.

## Repository Structure

- `data/processed/`: compact processed count tables and donor-level pseudobulk matrices.
- `metadata/`: dataset registry, sample metadata, annotation maps, and gene-set definitions.
- `results/tables/`: differential-expression results, robustness checks, sensitivity analyses, and cross-context outputs.
- `results/figures/`: selected generated figures.
- `scripts/`: preprocessing, annotation, pseudobulk, result-table, and figure-generation scripts.
- `scripts/r/`: R/edgeR and cameraPR analyses.
- `scripts/repro/`: focused reproducibility and boundary-check analyses.
- `docker/`, `docker-compose.yml`, `requirements.txt`: container and dependency definitions.
- `checks/`: file manifest with SHA256 hashes.

## Main Reproducibility Boundary

The key disease contrast uses donor-level pseudobulk profiles rather than treating cells as independent disease replicates. The diabetic erectile dysfunction group is small, so the archive includes exact label permutation, cross-engine checks, stress attenuation, cell-type specificity, downsampling, and boundary analyses to make the evidentiary limits transparent.

## Quick Start

Build the analysis container:

```bash
docker compose build
```

Run R analyses in the container:

```bash
docker run --rm -v "$PWD":/work -w /work donor-aware-corpus-cavernosum-transcriptomics:latest Rscript scripts/repro/05_edger_diabetes_specificity.R
```

Run Python analyses in the container:

```bash
docker run --rm -v "$PWD":/work -w /work donor-aware-corpus-cavernosum-transcriptomics:latest python scripts/repro/01_reframe_pericyte_de_gse206528.py
```

See `REPRODUCIBILITY.md` for the rerun map and `DATA_AVAILABILITY.md` for source-data details.

## License

The archived repository contents are released under CC BY 4.0; see `LICENSE`. This license applies to the files in this repository and does not change the terms of use for third-party source datasets.

## Citation

Use `CITATION.cff` for repository citation metadata. After GitHub-Zenodo archiving, cite the Zenodo DOI together with the GitHub repository URL.
