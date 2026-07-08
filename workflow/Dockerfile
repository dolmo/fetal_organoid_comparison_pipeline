# Start with a pure, stable Ubuntu Linux base
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    software-properties-common dirmngr wget ca-certificates \
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

# ALL the C++ system dependencies we discovered!
# (Notice we swapped pkg-config for pkgconf here)
RUN apt-get update && apt-get install -y \
    r-base r-base-dev \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libglpk-dev libgmp3-dev libgsl-dev libbz2-dev liblzma-dev \
    pkgconf libwebp-dev libcairo2-dev libmagick++-dev \
    python3-pip python3-dev git \
    && apt-get clean

RUN pip3 install --no-cache-dir anndata scanpy pandas numpy awscli

RUN R -e "install.packages(c('BiocManager', 'devtools', 'remotes', 'reticulate'))"

RUN R -e "BiocManager::install(c('DESeq2', 'limma', 'edgeR', 'ComplexHeatmap', \
    'clusterProfiler', 'GSVA', 'biomaRt', 'S4Vectors', 'GenomicRanges', \
    'SummarizedExperiment', 'SingleCellExperiment', 'miloR', 'rhdf5'))"

RUN R -e "remotes::install_github('scverse/anndataR')"
RUN R -e "devtools::install_github('ma-jacques/RIMA')"

ENV RETICULATE_PYTHON=/usr/bin/python3
WORKDIR /app