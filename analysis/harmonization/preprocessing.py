"""Data loading, metadata alignment, and gene filtering for the harmonization stage."""

import numpy as np
import pandas as pd
import scanpy as sc

import config


def load_datasets(fetal_path, organoid_path):
    """Read the raw fetal MetaAtlas and organoid .h5ad files."""
    adata_fetal = sc.read_h5ad(fetal_path)
    adata_organoid = sc.read_h5ad(organoid_path)
    return adata_fetal, adata_organoid


def subsample(adata, n_cells, seed=None):
    """Randomly downsample an AnnData to ``n_cells`` (used for fast test runs)."""
    if seed is not None:
        np.random.seed(seed)
    n = min(n_cells, adata.n_obs)
    idx = np.random.choice(adata.n_obs, n, replace=False)
    return adata[idx].copy()


# Organoid cell-type labels -> shared MetaAtlas (Type.v1) vocabulary. Labels not
# in this mapping (e.g. "panNeuronal") keep their original name.
LABEL_MAPPING = {
    "Excitatory-Neurons": "Excitatory Neuron",
    "Interneurons": "Inhibitory Neuron",
    "Radial-Glia": "RG",
    "Intermediate-Progenitor-Cells": "IPC",
    "Astrocytes": "Astrocyte",
    "Cajal-Retzius-Neurons": "CR",
    "Neural-Crest": "Neural Crest",
}


def align_metadata(adata_fetal, adata_organoid):
    """Harmonize ``obs`` columns across datasets so both share ``dataset``,
    ``Indvd``, ``unified_celltype``, age fields, and ECM condition for the
    downstream RIMA stage."""
    adata_organoid.obs["dataset"] = "Organoid"
    adata_fetal.obs["dataset"] = "MetaAtlas"

    # Rename organoid columns to the names the R stage expects.
    adata_organoid.obs.rename(columns={"basal-lamina": "ecm", "age": "day"}, inplace=True)

    # --- Organoid prep ---
    # 'sampleID' (e.g. rnh30-h9-003) is the batch key for scVI correction.
    adata_organoid.obs["Indvd"] = adata_organoid.obs["sampleID"].astype(str)
    adata_organoid.obs["unified_celltype"] = (
        adata_organoid.obs["cell-type"].astype(str).map(LABEL_MAPPING)
        .fillna(adata_organoid.obs["cell-type"].astype(str))
    )
    # Clean the organoid 'day' column to be purely numeric ("day084" -> 84.0).
    cleaned_days = adata_organoid.obs["day"].astype(str).str.replace("day", "", regex=False)
    adata_organoid.obs["day_num"] = pd.to_numeric(cleaned_days, errors="coerce")

    keep_org = ["dataset", "Indvd", "ecm", "day", "day_num", "unified_celltype", "cell-type"]
    keep_org = [col for col in keep_org if col in adata_organoid.obs.columns]
    adata_organoid.obs = adata_organoid.obs[keep_org]
    for col in ["Gestational_week", "Gest_week_num", "Type.v1", "State"]:
        if col not in adata_organoid.obs.columns:
            adata_organoid.obs[col] = np.nan

    # --- Fetal prep ---
    adata_fetal.obs["unified_celltype"] = adata_fetal.obs["Type.v1"].astype(str)
    cleaned_gw = adata_fetal.obs["Gestational_week"].astype(str).replace("postnatal", "40")
    adata_fetal.obs["Gest_week_num"] = pd.to_numeric(cleaned_gw, errors="coerce")

    keep_fetal = ["dataset", "Indvd", "Gestational_week", "Gest_week_num",
                  "Type.v1", "State", "unified_celltype"]
    keep_fetal = [col for col in keep_fetal if col in adata_fetal.obs.columns]
    adata_fetal.obs = adata_fetal.obs[keep_fetal]
    for col in ["ecm", "day", "day_num"]:
        if col not in adata_fetal.obs.columns:
            adata_fetal.obs[col] = np.nan

    return adata_fetal, adata_organoid


def intersect_and_filter_genes(adata_fetal, adata_organoid, min_cells=None):
    """Restrict both datasets to their shared genes and drop rarely-expressed
    genes. RIMA requires both datasets to have the exact same gene set."""
    if min_cells is None:
        min_cells = config.MIN_CELLS_PER_GENE

    common_genes = adata_fetal.var_names.intersection(adata_organoid.var_names)
    adata_fetal = adata_fetal[:, common_genes].copy()
    adata_organoid = adata_organoid[:, common_genes].copy()

    sc.pp.filter_genes(adata_fetal, min_cells=min_cells)
    sc.pp.filter_genes(adata_organoid, min_cells=min_cells)

    # Re-intersect in case filtering dropped different genes in each dataset.
    common_genes_filtered = adata_fetal.var_names.intersection(adata_organoid.var_names)
    adata_fetal = adata_fetal[:, common_genes_filtered].copy()
    adata_organoid = adata_organoid[:, common_genes_filtered].copy()

    return adata_fetal, adata_organoid
