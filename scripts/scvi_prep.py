import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import matplotlib.pyplot as plt
import os

# ==========================================
# 1. CONFIGURATION
# ==========================================
# The bash script downloads your raw data here:
FETAL_PATH = "/workspace/input/metaatlas_final_raw.h5ad"
ORGANOID_PATH = "/workspace/input/rnh027_complete_08072025.h5ad"

# The bash script tells the R script to look for inputs here, 
# so we will save the Python outputs back into the input folder!
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

# Fetal Subsample
n_fetal = min(50000, adata_fetal.n_obs) # Updated to 10k to match your print statement
idx_fetal = np.random.choice(adata_fetal.n_obs, n_fetal, replace=False)
adata_fetal_sub = adata_fetal[idx_fetal].copy()

# Organoid Subsample
n_organoid = min(50000, adata_organoid.n_obs) # Updated to 10k to match your print statement
idx_organoid = np.random.choice(adata_organoid.n_obs, n_organoid, replace=False)
adata_organoid_sub = adata_organoid[idx_organoid].copy()

# ==========================================
# 3. METADATA ALIGNMENT & UNIFICATION
# ==========================================
print("Aligning metadata and creating unified cell types...")
adata_organoid_sub.obs["dataset"] = "Organoid"
adata_fetal_sub.obs["dataset"] = "MetaAtlas"

# Maps organoid specific labels to the fetal 'Type.v1' labels
label_mapping = {
    "Excitatory_Neurons": "Excitatory Neuron",
    "Interneurons_(FOXG1-)": "Inhibitory Neuron",
    "Interneurons_(FOXG1+)": "Inhibitory Neuron",
    "RadialGlia": "RG",
    "IPCs": "IPC",
    "Newborn": "Newborn Neuron",
    "Dividing": "Div",
    "Astrocytes_oRG": "Astrocyte",
    "PAX3+": "PAX3+",
    "Progenitors_(FOXG1-)": "Progenitors_(FOXG1-)",
    "ChP_Hem": "ChP_Hem"
}

# --- Organoid Prep ---
adata_organoid_sub.obs["Indvd"] = (
    "org_" + 
    adata_organoid_sub.obs["replicate"].astype(str) + "_" + 
    adata_organoid_sub.obs["ecm"].astype(str)
)

adata_organoid_sub.obs["unified_celltype"] = (
    adata_organoid_sub.obs["type_rnh"].astype(str).map(label_mapping)
    .fillna(adata_organoid_sub.obs["type_rnh"].astype(str))
)

keep_org = ["dataset", "Indvd", "ecm", "day", "unified_celltype", "type_rnh", "cell_pred"]
keep_org = [col for col in keep_org if col in adata_organoid_sub.obs.columns]
adata_organoid_sub.obs = adata_organoid_sub.obs[keep_org]

for col in ["Gestational_week", "Type.v1", "State"]:
    if col not in adata_organoid_sub.obs.columns:
        adata_organoid_sub.obs[col] = np.nan

# --- Fetal Prep ---
# Now pulling from Type.v1!
adata_fetal_sub.obs["unified_celltype"] = adata_fetal_sub.obs["Type.v1"].astype(str)

keep_fetal = ["dataset", "Indvd", "Gestational_week", "Type.v1", "State", "unified_celltype"]
keep_fetal = [col for col in keep_fetal if col in adata_fetal_sub.obs.columns]
adata_fetal_sub.obs = adata_fetal_sub.obs[keep_fetal]

for col in ["ecm", "day"]:
    if col not in adata_fetal_sub.obs.columns:
        adata_fetal_sub.obs[col] = np.nan

# ==========================================
# 4. CONCATENATION & FILTERING
# ==========================================
print("Concatenating datasets...")
adata = adata_fetal_sub.concatenate(
    adata_organoid_sub, 
    batch_key="dataset_concat", 
    index_unique="-" 
)

adata.obs_names_make_unique() 

print("Filtering genes expressed in fewer than 10 cells...")
sc.pp.filter_genes(adata, min_cells=10)

# ==========================================
# 5. HARMONIZATION WITH scVI
# ==========================================
print("Pre-processing: Selecting Highly Variable Genes (HVGs)...")
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=3000, 
    flavor="seurat_v3", 
    batch_key="dataset", 
    subset=True
)

print("Setting up scVI model...")
scvi.model.SCVI.setup_anndata(adata, batch_key="dataset")

model = scvi.model.SCVI(
    adata, 
    n_latent=30, 
    n_layers=2,
    dispersion="gene-batch" 
)

print("Training scVI model with early stopping...")
model.train(
    max_epochs=400,
    batch_size=1024, 
    plan_kwargs={"lr": 1e-3},
    early_stopping=True,
    early_stopping_patience=20
)

print("Extracting X_scVI latent representation...")
adata.obsm["X_scVI"] = model.get_latent_representation()

# Optional: Keep the UMAPs to verify scVI worked well before splitting
print("Computing UMAP for visualization...")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
sc.tl.umap(adata) 

sc.pl.umap(adata, color=["dataset"], title="scVI Integrated Latent Space", show=False)
plt.savefig(f"{OUTPUT_DIR}/umap_dataset.png", bbox_inches='tight', dpi=300)
plt.close()

sc.pl.umap(adata, color=["unified_celltype"], title="Unified Cell Types", show=False)
plt.savefig(f"{OUTPUT_DIR}/umap_celltype.png", bbox_inches='tight', dpi=300)
plt.close()

# ==========================================
# 6. SPLIT AND SAVE FOR RIMA (R SCRIPT)
# ==========================================
print("Splitting the harmonized dataset back into Fetal and Organoid...")

# Subset the combined anndata object based on the dataset column
adata_fetal_out = adata[adata.obs['dataset'] == 'MetaAtlas'].copy()
adata_organoid_out = adata[adata.obs['dataset'] == 'Organoid'].copy()

# Construct output paths
fetal_out_path = os.path.join(OUTPUT_DIR, "metaatlas_subsample.h5ad")
organoid_out_path = os.path.join(OUTPUT_DIR, "rnh027_subsample.h5ad")

print(f"Saving Fetal data to {fetal_out_path}...")
adata_fetal_out.write_h5ad(fetal_out_path)

print(f"Saving Organoid data to {organoid_out_path}...")
adata_organoid_out.write_h5ad(organoid_out_path)

print("Python pre-processing complete! You can now run your RIMA R script.")