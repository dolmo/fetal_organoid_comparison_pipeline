"""Configuration for the scVI harmonization stage (Phase 1).

All file paths and model hyperparameters used by ``run_harmonization.py`` live
here so they can be changed in one place without touching the pipeline logic.
"""

# --- Input / output paths -------------------------------------------------
# These default to the container mount points used by the Kubernetes job
# (see ../../workflow/run_pipeline_local.sh). Override for local runs.
FETAL_PATH = "/workspace/input/metaatlas_final_raw.h5ad"
ORGANOID_PATH = "/workspace/input/rnh027_rnh030_complete.h5ad"
OUTPUT_DIR = "/workspace/input"

# Output filenames (consumed by the R RIMA stage, run_rima.R).
FETAL_OUT_NAME = "metaatlas_subsample.h5ad"
ORGANOID_OUT_NAME = "rnh027_rnh030_subsample.h5ad"

# --- Reproducibility ------------------------------------------------------
RANDOM_SEED = 42

# --- Optional subsampling (used during pipeline development/testing) ------
# The full analysis runs on the complete datasets (SUBSAMPLE = False). Set
# SUBSAMPLE = True to randomly downsample each dataset to SUBSAMPLE_N cells
# for a fast test run.
SUBSAMPLE = False
SUBSAMPLE_N = 10000

# --- Gene selection -------------------------------------------------------
MIN_CELLS_PER_GENE = 10
N_TOP_HVGS = 3000
HVG_FLAVOR = "cell_ranger"

# --- scVI model hyperparameters ------------------------------------------
SCVI_BATCH_KEY = "Indvd"
N_LATENT = 30
N_LAYERS = 2
DISPERSION = "gene-batch"
MAX_EPOCHS = 200
BATCH_SIZE = 1024
LEARNING_RATE = 1e-3
EARLY_STOPPING = True
EARLY_STOPPING_PATIENCE = 20
