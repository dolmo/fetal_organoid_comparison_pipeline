# install_r_packages.R
# Installs the R dependencies for the RIMA stage (analysis/rima/). Run once to
# set up a local R environment:  Rscript env/install_r_packages.R
#
# (In the containerized pipeline these packages are baked into the Docker image
# instead -- see the Dockerfile.)

options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c("jsonlite", "rlang", "reticulate", "remotes"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.22")

# RIMA / plotting dependencies
install.packages("igraph")
install.packages("tidygraph")
install.packages("ggraph")
install.packages("graphlayouts")
install.packages("ggplot2")   # plotting (used throughout analysis/rima)
install.packages("dtw")       # DTW temporal alignment (dtw_alignment.R)
install.packages("miloR")

BiocManager::install(c(
    "anndataR",
    "rhdf5",
    "SingleCellExperiment",
    "miloR"
))

remotes::install_github("ma-jacques/RIMA")
