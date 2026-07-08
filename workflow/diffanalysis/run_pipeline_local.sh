#!/bin/bash

# Define paths (using the high-speed /data mount)
PROJECT_ROOT="/data"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"
INPUT_S3_PATH="s3://braingeneersdev/dolmo/original/data/"
S3_ENDPOINT_URL_RESULTS="s3://braingeneersdev/dolmo/original/outputs/validation/"
S3_URL="https://s3-west.nrp-nautilus.io"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PROJECT_ROOT}/input"

LOG_FILE="${OUTPUT_DIR}/validation_log.txt"
echo "Starting SLUG-noid Validation..." > "$LOG_FILE"

# ---------------------------------------------------------
# PHASE 1: DOWNLOAD ONLY THE ORGANOID DATA
# ---------------------------------------------------------
echo "Downloading organoid data from S3..." | tee -a "$LOG_FILE"
# Directly copy just the one file we need to save bandwidth and time
aws s3 cp "${INPUT_S3_PATH}rnh027_rnh030_complete.h5ad" "${PROJECT_ROOT}/input/rnh027_complete_08072025.h5ad" --endpoint-url "${S3_URL}" >> "$LOG_FILE" 2>&1
DL_STATUS=$?

if [ $DL_STATUS -ne 0 ]; then
    echo "ERROR: S3 Download failed!" | tee -a "$LOG_FILE"
    aws s3 cp "${OUTPUT_DIR}" "${S3_ENDPOINT_URL_RESULTS}" --recursive --endpoint-url "${S3_URL}"
    exit $DL_STATUS
fi

# ---------------------------------------------------------
# PHASE 2: LIGHTWEIGHT PYTHON ANALYSIS
# ---------------------------------------------------------
echo "Installing Python dependencies..." | tee -a "$LOG_FILE"
pip3 install scanpy pandas matplotlib >> "$LOG_FILE" 2>&1

echo "Running Validation Script..." | tee -a "$LOG_FILE"
python3 /opt/scripts/run_de.py >> "$LOG_FILE" 2>&1
PY_STATUS=$?

if [ $PY_STATUS -ne 0 ]; then
    echo "ERROR: Python script crashed!" | tee -a "$LOG_FILE"
else
    echo "SUCCESS: Validation complete!" | tee -a "$LOG_FILE"
fi

# ---------------------------------------------------------
# PHASE 3: UPLOAD RESULTS
# ---------------------------------------------------------
echo "Uploading results to S3..."
aws s3 cp "${OUTPUT_DIR}" "${S3_ENDPOINT_URL_RESULTS}" --recursive --endpoint-url "${S3_URL}"

exit $PY_STATUS