options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c("jsonlite", "rlang", "reticulate", "remotes"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.22")


# Install igraph and other RIMA dependencies first
install.packages("igraph") 
install.packages("tidygraph")
install.packages("ggraph")
install.packages("graphlayouts")
install.packages("miloR")

BiocManager::install(c(
    "anndataR",
    "rhdf5",
    "SingleCellExperiment",
    "miloR"
))



remotes::install_github("ma-jacques/RIMA")