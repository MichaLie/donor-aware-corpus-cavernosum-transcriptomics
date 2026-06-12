#!/usr/bin/env bash
# Fetch public external source datasets from GEO. Run from the repository root.
# Large raw archives are not stored in this repository.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="${ROOT}/data/raw"
mkdir -p "${RAW}"

fetch() {
  local url="$1"
  local out="$2"
  mkdir -p "$(dirname "${out}")"
  if [[ -s "${out}" ]]; then
    echo "[skip] ${out}"
  else
    echo "[fetch] ${out}"
    curl -fL --retry 3 -o "${out}" "${url}"
  fi
}

fetch_raw_tar() {
  local acc="$1"
  local series="$2"
  local dir="${RAW}/${acc}"
  local tar="${dir}/${acc}_RAW.tar"
  mkdir -p "${dir}"
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/${series}/${acc}/suppl/${acc}_RAW.tar" "${tar}"
  tar -xf "${tar}" -C "${dir}"
}

fetch_raw_tar "GSE206528" "GSE206nnn"
fetch_raw_tar "GSE259348" "GSE259nnn"
fetch_raw_tar "GSE259299" "GSE259nnn"
fetch_raw_tar "GSE228900" "GSE228nnn"
fetch_raw_tar "GSE146500" "GSE146nnn"
fetch_raw_tar "GSE261085" "GSE261nnn"
fetch_raw_tar "GSE131882" "GSE131nnn"

fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126005/suppl/GSE126005_RawCountsTable.txt.gz" \
  "${RAW}/GSE126005/GSE126005_RawCountsTable.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126005/suppl/GSE126005_EdgeRTable.txt.gz" \
  "${RAW}/GSE126005/GSE126005_EdgeRTable.txt.gz"

echo "[note] HRA010374 is controlled access and is not downloaded by this script."
echo "[note] Open HRA010374-related supplement context: https://www.frontiersin.org/articles/10.3389/fendo.2025.1671482/full#supplementary-material"
echo "[note] See DATA_AVAILABILITY.md for source access details and included compact derivatives."
