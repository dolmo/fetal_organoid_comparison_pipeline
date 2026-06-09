import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import matplotlib.pyplot as plt
import os

# ==========================================
# 1. CONFIGURATION
# ==========================================
FETAL_PATH = "/workspace/input/metaatlas_final_raw.h5ad"
ORGANOID_PATH = "/workspace/input/rnh027_rnh030_complete.h5ad"
OUTPUT_DIR = "/workspace/input"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# 2. DATA LOADING & SUBSAMPLING
# ==========================================
print("Loading datasets...")
adata_fetal = sc.read_h5ad(FETAL_PATH)
adata_organoid = sc.read_h5ad(ORGANOID_PATH)

print("Subsampling datasets to 10,000 cells each...")
np.random.seed(42)

"""# Fetal Subsample
n_fetal = min(10000, adata_fetal.n_obs) # Updated to 10k to match your print statement
idx_fetal = np.random.choice(adata_fetal.n_obs, n_fetal, replace=False)
adata_fetal_sub = adata_fetal[idx_fetal].copy()

# Organoid Subsample
n_organoid = min(10000, adata_organoid.n_obs) # Updated to 10k to match your print statement
idx_organoid = np.random.choice(adata_organoid.n_obs, n_organoid, replace=False)
adata_organoid_sub = adata_organoid[idx_organoid].copy()"""

print("Copying full datasets for processing...")
adata_fetal_sub = adata_fetal.copy()
adata_organoid_sub = adata_organoid.copy()
# ==========================================
# 3. METADATA ALIGNMENT & UNIFICATION
# ==========================================
print("Aligning metadata and creating unified cell types...")
adata_organoid_sub.obs["dataset"] = "Organoid"
adata_fetal_sub.obs["dataset"] = "MetaAtlas"

# We are renaming your new columns back to the old names so the R script 
# doesn't crash and requires zero changes!
adata_organoid_sub.obs.rename(columns={"basal-lamina": "ecm", "age": "day"}, inplace=True)

# Based on your new metadata.txt, here is the updated mapping to Fetal Type.v1
label_mapping = {
    "Excitatory-Neurons": "Excitatory Neuron",
    "Interneurons": "Inhibitory Neuron",
    "Radial-Glia": "RG",
    "Intermediate-Progenitor-Cells": "IPC",
    "Astrocytes": "Astrocyte",
    "Cajal-Retzius-Neurons": "CR",
    "Neural-Crest": "Neural Crest"
}

# --- Organoid Prep ---
# Your new metadata has a perfect 'sampleID' column (e.g., rnh30-h9-003), 
# which is perfect for scVI batch correction.
adata_organoid_sub.obs["Indvd"] = adata_organoid_sub.obs["sampleID"].astype(str)

# Apply the mapping. If a cell type isn't in the dictionary (like 'panNeuronal'), 
# it just keeps its original name.
adata_organoid_sub.obs["unified_celltype"] = (
    adata_organoid_sub.obs["cell-type"].astype(str).map(label_mapping)
    .fillna(adata_organoid_sub.obs["cell-type"].astype(str))
)

# Clean the organoid 'day' column to be purely numeric ("day084" -> 84.0)
cleaned_days = adata_organoid_sub.obs["day"].astype(str).str.replace("day", "", regex=False)
adata_organoid_sub.obs["day_num"] = pd.to_numeric(cleaned_days, errors='coerce')

keep_org = ["dataset", "Indvd", "ecm", "day", "day_num", "unified_celltype", "cell-type"]
keep_org = [col for col in keep_org if col in adata_organoid_sub.obs.columns]
adata_organoid_sub.obs = adata_organoid_sub.obs[keep_org]

for col in ["Gestational_week", "Gest_week_num", "Type.v1", "State"]:
    if col not in adata_organoid_sub.obs.columns:
        adata_organoid_sub.obs[col] = np.nan

# --- Fetal Prep ---
adata_fetal_sub.obs["unified_celltype"] = adata_fetal_sub.obs["Type.v1"].astype(str)

cleaned_gw = adata_fetal_sub.obs["Gestational_week"].astype(str).replace("postnatal", "40")
adata_fetal_sub.obs["Gest_week_num"] = pd.to_numeric(cleaned_gw, errors='coerce')

keep_fetal = ["dataset", "Indvd", "Gestational_week", "Gest_week_num", "Type.v1", "State", "unified_celltype"]
keep_fetal = [col for col in keep_fetal if col in adata_fetal_sub.obs.columns]
adata_fetal_sub.obs = adata_fetal_sub.obs[keep_fetal]

for col in ["ecm", "day", "day_num"]:
    if col not in adata_fetal_sub.obs.columns:
        adata_fetal_sub.obs[col] = np.nan

# ==========================================
# 4. GENE INTERSECTION & FILTERING (RIMA Prep)
# ==========================================
print("Finding common genes across Fetal and Organoid datasets...")
# RIMA needs both datasets to have the exact same genes to run Spearman correlations
common_genes = adata_fetal_sub.var_names.intersection(adata_organoid_sub.var_names)

adata_fetal_sub = adata_fetal_sub[:, common_genes].copy()
adata_organoid_sub = adata_organoid_sub[:, common_genes].copy()

print("Filtering genes expressed in fewer than 10 cells independently...")
sc.pp.filter_genes(adata_fetal_sub, min_cells=10)
sc.pp.filter_genes(adata_organoid_sub, min_cells=10)

# Re-intersect in case filtering dropped different genes in different datasets
common_genes_filtered = adata_fetal_sub.var_names.intersection(adata_organoid_sub.var_names)
adata_fetal_sub = adata_fetal_sub[:, common_genes_filtered].copy()
adata_organoid_sub = adata_organoid_sub[:, common_genes_filtered].copy()

# ==========================================
# 5. INDEPENDENT HARMONIZATION WITH scVI
# ==========================================
print("Pre-processing: Selecting Highly Variable Genes (HVGs)...")
# FIX: Removed batch_key="Indvd" to avoid LOESS singularity crashes on small subsampled batches.
# Changed flavor to "cell_ranger" to handle float matrices gracefully.
sc.pp.highly_variable_genes(adata_fetal_sub, n_top_genes=3000, flavor="cell_ranger")
sc.pp.highly_variable_genes(adata_organoid_sub, n_top_genes=3000, flavor="cell_ranger")

# Combine HVGs into a master list so we capture drivers of variation from BOTH systems
hvg_fetal = adata_fetal_sub.var[adata_fetal_sub.var['highly_variable']].index
hvg_organoid = adata_organoid_sub.var[adata_organoid_sub.var['highly_variable']].index
hvg_union = hvg_fetal.union(hvg_organoid)

# Subset down to the unified highly variable gene list
adata_fetal_sub = adata_fetal_sub[:, hvg_union].copy()
adata_organoid_sub = adata_organoid_sub[:, hvg_union].copy()

print("Setting up and training scVI model for FETAL data...")
# Batch key is Indvd to fix internal variation, leaving cross-species alone
scvi.model.SCVI.setup_anndata(adata_fetal_sub, batch_key="Indvd")
model_fetal = scvi.model.SCVI(adata_fetal_sub, n_latent=30, n_layers=2, dispersion="gene-batch")
model_fetal.train(max_epochs=200, batch_size=1024, plan_kwargs={"lr": 1e-3}, early_stopping=True, early_stopping_patience=20)
adata_fetal_sub.obsm["X_scVI"] = model_fetal.get_latent_representation()

print("Setting up and training scVI model for ORGANOID data...")
# Batch key is Indvd to fix internal variation, leaving cross-species alone
scvi.model.SCVI.setup_anndata(adata_organoid_sub, batch_key="Indvd")
model_organoid = scvi.model.SCVI(adata_organoid_sub, n_latent=30, n_layers=2, dispersion="gene-batch")
model_organoid.train(max_epochs=200, batch_size=1024, plan_kwargs={"lr": 1e-3}, early_stopping=True, early_stopping_patience=20)
adata_organoid_sub.obsm["X_scVI"] = model_organoid.get_latent_representation()

# ==========================================
# 6. SAVE FOR RIMA (R SCRIPT)
# ==========================================
fetal_out_path = os.path.join(OUTPUT_DIR, "metaatlas_subsample.h5ad")
organoid_out_path = os.path.join(OUTPUT_DIR, "rnh027_rnh030_subsample.h5ad")

print(f"Saving independent Fetal data to {fetal_out_path}...")
adata_fetal_sub.write_h5ad(fetal_out_path)

print(f"Saving independent Organoid data to {organoid_out_path}...")
adata_organoid_sub.write_h5ad(organoid_out_path)

print("Python pre-processing complete! You can now run your RIMA R script.")