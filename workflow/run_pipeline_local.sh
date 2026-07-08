#!/bin/bash

# Hardcode the paths directly here
PROJECT_ROOT="/workspace"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"
INPUT_S3_PATH="s3://braingeneersdev/dolmo/original/data/"
S3_ENDPOINT_URL_RESULTS="${S3_ENDPOINT_URL_RESULTS:-s3://braingeneersdev/dolmo/original/outputs/}"
S3_URL="https://s3-west.nrp-nautilus.io"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PROJECT_ROOT}/input"

LOG_FILE="${OUTPUT_DIR}/pipeline_crash_report.txt"

echo "Starting Full Integration Pipeline on Cluster..." > "$LOG_FILE"

echo "Installing awscli to bypass GitHub..."
pip install awscli


# ---------------------------------------------------------
# PHASE 0: DOWNLOAD
# ---------------------------------------------------------
echo "Downloading raw input data from S3 using aws cli..." | tee -a "$LOG_FILE"

# 1. Download the Organoid dataset
aws s3 cp "${INPUT_S3_PATH}rnh027_rnh030_complete.h5ad" "${PROJECT_ROOT}/input/" --endpoint-url "${S3_URL}" >> "$LOG_FILE" 2>&1
DL_STATUS=$?

# 2. Download the Fetal Brain MetaAtlas (only if the first download succeeded)
if [ $DL_STATUS -eq 0 ]; then
    aws s3 cp "${INPUT_S3_PATH}metaatlas_final_raw.h5ad" "${PROJECT_ROOT}/input/" --endpoint-url "${S3_URL}" >> "$LOG_FILE" 2>&1
    DL_STATUS=$?
fi

# If either download fails, STOP the script immediately. Do not run Python on broken data.
if [ $DL_STATUS -ne 0 ]; then
    echo "ERROR: S3 Download failed or timed out! Uploading logs and exiting." | tee -a "$LOG_FILE"
    aws s3 cp "$LOG_FILE" "${S3_ENDPOINT_URL_RESULTS}" --endpoint-url "${S3_URL}"
    exit $DL_STATUS
fi

# ---------------------------------------------------------
# PHASE 1: PYTHON scVI HARMONIZATION
# ---------------------------------------------------------
echo "Installing Python dependencies..." | tee -a "$LOG_FILE"
pip3 install scvi-tools scanpy pandas numpy matplotlib scikit-misc zarr >> "$LOG_FILE" 2>&1

echo "Running Phase 1: Python scVI Harmonization..." | tee -a "$LOG_FILE"
# Point to the injected ConfigMap path
python3 /opt/scripts/run_harmonization.py >> "$LOG_FILE" 2>&1
PY_STATUS=$?

# If Python crashes, upload the logs and stop the job immediately
if [ $PY_STATUS -ne 0 ]; then
    echo "ERROR: Python script crashed! Uploading logs to S3..." | tee -a "$LOG_FILE"
    s5cmd --endpoint-url "${S3_URL}" cp "${OUTPUT_DIR}/*" "${S3_ENDPOINT_URL_RESULTS}"
    exit $PY_STATUS
fi

# Move the Python UMAP plots to the output folder so they get uploaded
mv ${PROJECT_ROOT}/input/*.png ${OUTPUT_DIR}/ 2>/dev/null

# ---------------------------------------------------------
# PHASE 2: RIMA R ANALYSIS
# ---------------------------------------------------------
echo "Running Phase 2: RIMA R Analysis..." | tee -a "$LOG_FILE"
Rscript /opt/scripts/run_rima.R --input "${PROJECT_ROOT}/input" --output "${OUTPUT_DIR}" >> "$LOG_FILE" 2>&1
R_STATUS=$?

if [ $R_STATUS -ne 0 ]; then
    echo "ERROR: R script crashed! Uploading logs to S3..." | tee -a "$LOG_FILE"
else
    echo "SUCCESS: Entire pipeline finished perfectly!" | tee -a "$LOG_FILE"
fi

# ---------------------------------------------------------
# PHASE 3: UPLOAD RESULTS
# ---------------------------------------------------------
echo "Uploading all results, plots, and logs back to S3..." | tee -a "$LOG_FILE"

aws s3 sync "${OUTPUT_DIR}/" "${S3_ENDPOINT_URL_RESULTS}" --endpoint-url "${S3_URL}" >> "$LOG_FILE" 2>&1

# Tell Kubernetes if the pod succeeded or failed
exit $R_STATUS