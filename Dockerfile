FROM rocker/tidyverse:4.4.0

# 1. System dependencies (including the bz2 and lzma fixes)
RUN apt-get update && apt-get install -y \
    libxml2-dev libssl-dev libcurl4-openssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libglpk-dev libgmp3-dev libgsl-dev libbz2-dev liblzma-dev \
    python3-pip python3-dev git \
    && apt-get clean

# 2. Python and AWS CLI
RUN pip3 install --no-cache-dir anndata scanpy pandas numpy awscli

# 3. Base R Package Managers
RUN R -e "install.packages(c('BiocManager', 'devtools', 'remotes', 'reticulate'))"

ENV RETICULATE_PYTHON=/usr/bin/python3
WORKDIR /apps