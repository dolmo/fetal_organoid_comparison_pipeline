# neighborhoods.R
# Milo neighbourhood construction and cleanup on the scVI latent space.

# Build a Milo object and define refined neighbourhoods on the KNN graph.
define_neighbourhoods <- function(sce, prop_seeds, knn = 10, reduced.dim = "X_scVI") {
  n_components <- ncol(reducedDim(sce, reduced.dim))
  mi <- Milo(sce)

  mi <- miloR::buildGraph(mi, k = knn, d = n_components, reduced.dim = reduced.dim)
  mi <- miloR::makeNhoods(mi, prop = prop_seeds, k = knn, d = n_components,
                          reduced_dims = reduced.dim, refined = TRUE)
  mi
}

# Diagnostic print of the preprocessed Milo objects (cells, genes, nhoods, dims).
inspect_milos <- function(milos) {
  for (i in seq_along(milos)) {
    print(paste("=== Post-preprocess Milo", i, "==="))
    print(paste("  class:", class(milos[[i]])))
    print(paste("  ncol (cells):", ncol(milos[[i]])))
    print(paste("  nrow (genes):", nrow(milos[[i]])))
    print(paste("  colData rows:", nrow(colData(milos[[i]]))))
    print(paste("  nhoods dim:", paste(dim(nhoods(milos[[i]])), collapse = " x ")))
    print(paste("  reducedDimNames:", paste(reducedDimNames(milos[[i]]), collapse = ", ")))
    if ("X_scVI" %in% reducedDimNames(milos[[i]])) {
      print(paste("  X_scVI dim:", paste(dim(reducedDim(milos[[i]], "X_scVI")), collapse = " x ")))
    } else {
      print("  WARNING: X_scVI not found in reducedDims!")
    }
  }
}

# Replace missing/NA/empty cell-type labels with "Unknown" so the significance
# scrambling step does not crash. Returns the modified list of Milo objects.
patch_missing_celltypes <- function(milos, col = "unified_celltype") {
  print("Patching missing cell types to prevent scrambling crash...")
  for (i in seq_along(milos)) {
    ct <- as.character(colData(milos[[i]])[[col]])
    n_na <- sum(is.na(ct) | ct == "nan" | ct == "")
    print(paste("  Milo", i, "- patching", n_na, "missing/nan/empty cell types"))
    ct[is.na(ct) | ct == "nan" | ct == ""] <- "Unknown"
    colData(milos[[i]])[[col]] <- ct
  }
  milos
}
