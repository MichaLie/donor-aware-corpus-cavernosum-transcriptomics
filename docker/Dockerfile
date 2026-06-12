FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:b10002b39efa30c3779ad839549806ebdbb29b3266f0d2428478b04426e55929

ENV DEBIAN_FRONTEND=noninteractive
ENV RETICULATE_PYTHON=/usr/bin/python3

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('remotes', 'data.table', 'ggplot2', 'here', 'patchwork', 'cowplot', 'pheatmap', 'ggrepel'), repos='https://cloud.r-project.org', Ncpus=max(1, parallel::detectCores()-1))"

RUN R -q -e "BiocManager::install(c('SingleCellExperiment', 'SummarizedExperiment', 'scater', 'scran', 'scDblFinder', 'edgeR', 'limma', 'muscat', 'dreamlet', 'zellkonverter', 'batchelor', 'sva', 'ComplexHeatmap', 'fgsea'), ask=FALSE, update=FALSE, Ncpus=max(1, parallel::detectCores()-1))"

COPY requirements.txt /tmp/requirements.txt
RUN apt-get update && apt-get install -y --no-install-recommends \
    python-is-python3 \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"
ENV RETICULATE_PYTHON=/opt/venv/bin/python
RUN pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /work
