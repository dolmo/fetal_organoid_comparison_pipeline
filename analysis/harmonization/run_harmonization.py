#!/usr/bin/env python3
"""
scVI harmonization stage (Phase 1).

Loads the raw fetal MetaAtlas and brain-organoid scRNA-seq datasets, harmonizes
their metadata, selects shared highly-variable genes, and trains an independent
scVI model per dataset to produce a batch-corrected latent space (``X_scVI``).
The two harmonized .h5ad files are written to ``config.OUTPUT_DIR`` for the R
RIMA stage (run_rima.R) to consume.

Run:
    python3 run_harmonization.py

Paths and hyperparameters are configured in config.py.
"""

import os

import numpy as np

import config
from preprocessing import (
    load_datasets,
    subsample,
    align_metadata,
    intersect_and_filter_genes,
)
from scvi_training import select_hvgs, train_scvi


def main():
    os.makedirs(config.OUTPUT_DIR, exist_ok=True)
    np.random.seed(config.RANDOM_SEED)

    print("Loading datasets...")
    adata_fetal, adata_organoid = load_datasets(config.FETAL_PATH, config.ORGANOID_PATH)

    if config.SUBSAMPLE:
        print(f"Subsampling datasets to {config.SUBSAMPLE_N} cells each...")
        adata_fetal = subsample(adata_fetal, config.SUBSAMPLE_N)
        adata_organoid = subsample(adata_organoid, config.SUBSAMPLE_N)
    else:
        print("Using full datasets for processing...")

    print("Aligning metadata and creating unified cell types...")
    adata_fetal, adata_organoid = align_metadata(adata_fetal, adata_organoid)

    print("Finding common genes and filtering rare genes...")
    adata_fetal, adata_organoid = intersect_and_filter_genes(adata_fetal, adata_organoid)

    print("Selecting highly variable genes...")
    adata_fetal, adata_organoid = select_hvgs(adata_fetal, adata_organoid)

    print("Training scVI model for FETAL data...")
    adata_fetal = train_scvi(adata_fetal)

    print("Training scVI model for ORGANOID data...")
    adata_organoid = train_scvi(adata_organoid)

    fetal_out = os.path.join(config.OUTPUT_DIR, config.FETAL_OUT_NAME)
    organoid_out = os.path.join(config.OUTPUT_DIR, config.ORGANOID_OUT_NAME)

    print(f"Saving harmonized fetal data to {fetal_out}...")
    adata_fetal.write_h5ad(fetal_out)
    print(f"Saving harmonized organoid data to {organoid_out}...")
    adata_organoid.write_h5ad(organoid_out)

    print("Harmonization complete! Ready for the RIMA R stage (run_rima.R).")


if __name__ == "__main__":
    main()
