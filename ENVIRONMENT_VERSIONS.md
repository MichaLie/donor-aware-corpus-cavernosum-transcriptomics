# Environment Versions

This file records the analysis environment used to generate the archived outputs.

## Docker Base Image

- Base image: `bioconductor/bioconductor_docker:RELEASE_3_23`
- OCI index digest: `sha256:b10002b39efa30c3779ad839549806ebdbb29b3266f0d2428478b04426e55929`
- Linux/amd64 manifest digest: `sha256:f2dfb497dce909e04380799225cf62e6736a607d5a284e8783640cea8b33bd27`
- Linux/arm64 manifest digest: `sha256:a0ae0f11745eeba64b18ca57a4d8f46d15a07aa405f0a11ea6b89e3396f7b172`

The Dockerfile pins the OCI index digest.

## R and Bioconductor

Versions verified in the analysis container:

| Component | Version |
|---|---:|
| R | 4.6.0 |
| Bioconductor | 3.23 |
| edgeR | 4.10.0 |
| limma | 3.68.2 |
| data.table | 1.18.4 |
| ggplot2 | 4.0.3 |
| SingleCellExperiment | 1.34.0 |
| SummarizedExperiment | 1.42.0 |
| scater | 1.40.0 |
| scran | 1.40.0 |
| scDblFinder | 1.26.0 |
| muscat | 1.26.0 |
| dreamlet | 1.10.0 |
| zellkonverter | 1.22.0 |
| batchelor | 1.28.0 |
| sva | 3.60.0 |
| ComplexHeatmap | 2.28.0 |
| fgsea | 1.38.0 |

## Python

Versions verified in the analysis container:

| Package | Version |
|---|---:|
| Python | 3.12.3 |
| pandas | 2.3.3 |
| numpy | 2.4.5 |
| scipy | 1.17.1 |
| anndata | 0.12.14 |
| scanpy | 1.12.1 |
| scikit-learn | 1.8.0 |
| statsmodels | 0.14.6 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| h5py | 3.16.0 |
| pyarrow | 24.0.0 |
| gseapy | 1.2.1 |
| decoupler | 2.1.6 |
| scikit-misc | 0.5.2 |
| python-igraph | 1.0.0 |
| leidenalg | 0.11.0 |
| harmonypy | 2.0.0 |
| scrublet | 0.2.3 |
| openpyxl | 3.1.5 |


## Python (focused reanalysis additions)

The donor-aware reanalysis, replication, specificity, QC, and third-engine scripts
(`scripts/repro/`) were run with:

| Component | Version |
|---|---:|
| Python | 3.12 |
| scanpy | 1.12.1 |
| anndata | 0.12.10 |
| pydeseq2 | 0.5.4 |
| numpy | 2.2.6 |
| pandas | 2.3.3 |
| scipy | 1.16.3 |
| statsmodels | 0.14.6 |
| python-igraph | 1.0.0 |
| leidenalg | 0.12.0 |

R-based confirmation (edgeR-QL and limma-voom) used the Bioconductor container above (edgeR 4.10.0,
limma 3.68.2).
