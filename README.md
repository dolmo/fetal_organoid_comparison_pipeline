# Brain Organoid ↔ Fetal Brain Atlas Matching Pipeline

Computational pipeline for the project *"Mapping Extracellular Matrix and Chemical
Signaling Influence on Brain Organoid Maturation using Single-Cell Transcriptomics"*
(David Olmo Marchal).

It quantifies how far an extracellular-matrix (ECM) + basal-lamina protocol pushes
brain-organoid cells toward true fetal-brain maturation, by harmonizing an organoid
scRNA-seq dataset (Rh027–30) against a fetal-brain reference (MetaAtlas) and rigorously
matching cell-state neighbourhoods across the *in vitro / in vivo* barrier.

## Pipeline overview

```
raw .h5ad (fetal + organoid)
        │
        ▼
[Phase 1] scVI harmonization        analysis/harmonization/   (Python)
   metadata alignment, shared HVGs, per-dataset scVI batch correction (X_scVI)
        │  → metaatlas_subsample.h5ad, rnh027_rnh030_subsample.h5ad
        ▼
[Phase 2] RIMA cross-atlas analysis analysis/rima/            (R)
   Milo neighbourhoods → similarity → significance → matching → COPE
   → maturation age transfer → DTW temporal alignment + figures
        │
        ▼
[Phase 3] DE validation             analysis/differential_expression/  (Python)
   Wilcoxon test (ECM vs control) on targeted maturation genes + dot plot
```

## Repository layout

```
analysis/
  harmonization/            Phase 1 — scVI batch correction (Python)
    config.py                 paths + model hyperparameters
    preprocessing.py          load, (optional) subsample, metadata align, gene filter
    scvi_training.py          HVG selection + scVI training
    run_harmonization.py      entry point
  rima/                     Phase 2 — RIMA atlas matching (R)
    data_io.R                 load harmonized .h5ad as SingleCellExperiment
    neighborhoods.R           Milo neighbourhood construction + cleanup
    matching.R                similarity significance, matching, conserved genes (COPE)
    plots.R                   embedding / similarity-heatmap / paired-expression figures
    cell_extraction.R         pull cells from a matched neighbourhood pair
    maturation.R              transfer fetal age onto organoid neighbourhoods + plot
    dtw_alignment.R           per-cell-type Dynamic Time Warping alignment
    run_rima.R                entry point (sources the files above)
  differential_expression/  Phase 3 — DE validation (Python)
    run_de.py                 entry point
env/
  install_r_packages.R      R dependencies (RIMA, miloR, SingleCellExperiment, ...)
  requirements.txt          Python dependencies
workflow/                   Deployment (Docker + NRP/Kubernetes + S3); see below
logs/                       Saved outputs from previous runs (figures, tables)
Dockerfile                  Analysis container base image
```

## Environment setup

Python (Phase 1 + Phase 3):

```bash
pip install -r env/requirements.txt
```

R (Phase 2):

```bash
Rscript env/install_r_packages.R
```

A prebuilt Docker image with both stacks is defined in the `Dockerfile` (see also
`workflow/Dockerfile`), which is what the cluster jobs use.

## Running the pipeline

Edit the paths at the top of `analysis/harmonization/config.py` and
`analysis/differential_expression/run_de.py` to point at your local `.h5ad` files,
then run each phase in order:

```bash
# Phase 1 — harmonize (writes two *_subsample.h5ad files to config.OUTPUT_DIR)
python3 analysis/harmonization/run_harmonization.py

# Phase 2 — RIMA analysis (reads the harmonized files, writes tables + figures)
Rscript analysis/rima/run_rima.R --input <harmonized_dir> --output <results_dir>

# Phase 3 — differential expression validation
python3 analysis/differential_expression/run_de.py
```

Each R stage file is sourced by `run_rima.R` relative to its own location, so the
files in `analysis/rima/` must stay together.

## Cluster deployment

The pipeline was run on the NRP Nautilus Kubernetes cluster with data in AWS S3.
The orchestration lives in `workflow/`:

- `workflow/run_pipeline_local.sh` — downloads input from S3, runs Phase 1 (Python)
  then Phase 2 (R), uploads results.
- `workflow/diffanalysis/run_pipeline_local.sh` — runs the Phase 3 DE validation.
- `workflow/run_job.sh` + `jobdefinition.yaml` — submit the Kubernetes job.

The job mounts the analysis scripts flat into `/opt/scripts` via a ConfigMap. When
rebuilding that ConfigMap, include every file from the relevant `analysis/` stage,
e.g.:

```bash
kubectl create configmap rima-scripts-all \
  --from-file=analysis/harmonization/ --from-file=analysis/rima/
```

## Data availability

- **Rh027–30** (day-30 organoids ± ECM/basal lamina): generated internally, available
  on request to the PI.
- **MetaAtlas** (fetal-brain reference): compendium of 7 published fetal-brain scRNA-seq
  studies (Fan 2018/2020, Polioudakis 2019, Nowakowski 2017, Zhong 2018, Smith 2021,
  Bhaduri 2021).

## Notes

- `workflow/config.sh` holds your local AWS key env vars and is gitignored (never
  committed). It contains real credentials in plaintext on disk — rotate them if they
  may have been exposed.
- `workflow/` was originally a separate cloned template repo (`JesusGF1/nrp_template`);
  it has been folded into this repository, and the unused template files (Allen NWB
  downloader, boto3 S3 helpers, boilerplate docs) have been removed.
