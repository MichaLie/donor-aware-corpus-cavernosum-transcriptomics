# Reproducibility

This repository is organized for reproducibility in two layers:

1. Compact processed inputs and output tables are included for lightweight inspection.
2. Raw-to-results reruns can be performed by downloading source data from GEO and rebuilding intermediate AnnData objects.

## Environment

The recommended environment is Docker-based:

- `docker/Dockerfile`
- `docker-compose.yml`
- `requirements.txt`
- `ENVIRONMENT_VERSIONS.md`

Build the container:

```bash
docker compose build
```

Run R scripts:

```bash
docker run --rm -v "$PWD":/work -w /work donor-aware-corpus-cavernosum-transcriptomics:latest Rscript scripts/repro/05_edger_diabetes_specificity.R
```

Run Python scripts:

```bash
docker run --rm -v "$PWD":/work -w /work donor-aware-corpus-cavernosum-transcriptomics:latest python scripts/repro/01_reframe_pericyte_de_gse206528.py
```

## Lightweight Rerun From Included Inputs

The following analyses use included compact pseudobulk/count inputs and can be rerun without full `.h5ad` files:

1. `scripts/repro/01_reframe_pericyte_de_gse206528.py`
2. `scripts/repro/02_exact_permutation_diabetes_specificity.py`
3. `scripts/repro/05_edger_diabetes_specificity.R`
4. `scripts/repro/08_third_engine_limma_voom.R`

Primary outputs are written to `results/tables/repro/`.

## Full Rerun From Raw Data

A full raw-to-results rerun requires source data downloads and intermediate AnnData generation:

1. Download public source data:

   ```bash
   bash scripts/fetch_external_datasets.sh
   ```

2. Build `GSE206528` AnnData objects:

   ```bash
   python scripts/gse206528_convert_to_h5ad.py
   python scripts/gse206528_stream_qc.py
   python scripts/gse206528_build_combined_h5ad.py
   python scripts/gse206528_scanpy_preprocess.py
   ```

3. Generate pseudobulk inputs:

   ```bash
   python scripts/gse206528_make_pseudobulk.py
   python scripts/gse206528_make_mural_subtype_pseudobulk.py
   ```

4. Run donor-aware edgeR models:

   ```bash
   docker run --rm -v "$PWD":/work -w /work donor-aware-corpus-cavernosum-transcriptomics:latest Rscript scripts/r/gse206528_pseudobulk_edger.R
   docker run --rm -v "$PWD":/work -w /work donor-aware-corpus-cavernosum-transcriptomics:latest Rscript scripts/r/gse206528_mural_subtype_edger.R
   ```

5. Run focused reproducibility analyses:

   ```bash
   python scripts/repro/03_stress_corrected_pericyte_GSE206528.py
   python scripts/repro/04_gse259348_annotate_and_replicate.py
   python scripts/repro/04b_gse259348_mural_refine_replicate.py
   python scripts/repro/05_invitro_gse228900_hyperglycemia.py
   python scripts/repro/05_rat_negative.py
   python scripts/repro/05_spatial_gse261085.py
   Rscript scripts/repro/05a_inspect_gse131882_rds.R
   Rscript scripts/repro/05b_gse131882_cluster_pericyte_pseudobulk.R
   Rscript scripts/repro/05c_gse131882_refine_mural.R
   python scripts/repro/05d_gse131882_donor_level.py
   python scripts/repro/06_celltype_specificity.py
   python scripts/repro/07_qc_sensitivity.py
   ```

Some steps require large raw files or generated `.h5ad`/RDS-derived intermediate objects that are intentionally excluded from the repository.

## Interpretation Boundary

The key contrast is exploratory because it uses a small number of donors. The repository therefore includes exact label permutation, cross-engine differential-expression checks, stress attenuation, cell-type specificity, equal-cell downsampling, independent-donor directionality, and negative/boundary analyses.
