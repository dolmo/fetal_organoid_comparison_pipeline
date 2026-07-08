# data_io.R
# Load the scVI-harmonized .h5ad files (produced by run_harmonization.py) as
# SingleCellExperiment objects with a logcounts assay for RIMA.

load_sce <- function(h5ad_path) {
  print(paste("Loading", h5ad_path))
  adata <- anndataR::read_h5ad(h5ad_path)

  sce <- SingleCellExperiment(
    assays = list(counts = t(adata$X)),
    colData = as.data.frame(adata$obs),
    rowData = as.data.frame(adata$var),
    reducedDims = list(X_scVI = adata$obsm[["X_scVI"]])
  )

  # RIMA gene-expression analysis works on logcounts.
  assays(sce)$logcounts <- log1p(assays(sce)$counts)

  print(paste("  dims:", paste(dim(sce), collapse = " x "),
              "| reducedDims:", paste(reducedDimNames(sce), collapse = ", ")))
  sce
}

# Load both harmonized datasets from an input directory. Filenames match the
# outputs written by the Python harmonization stage.
load_datasets <- function(input_dir) {
  fetal_path <- file.path(input_dir, "metaatlas_subsample.h5ad")
  organoid_path <- file.path(input_dir, "rnh027_rnh030_subsample.h5ad")

  print("Loading datasets...")
  list(
    fetal = load_sce(fetal_path),
    organoid = load_sce(organoid_path)
  )
}
