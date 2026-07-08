#!/bin/bash
# Deploy Kubernetes job for SLUG-noid Validation

#=============================================================================
# JOB IDENTIFICATION
#=============================================================================
export NAME="slugnoid-validation"
export JOB_PREFIX="dolmo"
export CONTAINER_NAME="rima-image"

#=============================================================================
# DOCKER IMAGE
#=============================================================================
export DOCKER_IMAGE="dolmomar/rima-public:v6"

#=============================================================================
# AWS/S3 CREDENTIALS
#=============================================================================
export AWS_ACCESS_KEY_ID="${AWS_ACCESS_KEY_ID:-}"   
export AWS_SECRET_ACCESS_KEY="${AWS_SECRET_ACCESS_KEY:-}"
export S3_ENDPOINT="https://s3-west.nrp-nautilus.io"

#=============================================================================
# PATHS
#=============================================================================
export PROJECT_ROOT="/data"
export OUTPUT_DIR="/data/outputs"
export INPUT_S3_PATH="s3://braingeneersdev/dolmo/original/data/"
# Saving to a new folder so it doesn't overwrite your RIMA outputs
export S3_ENDPOINT_URL_RESULTS="s3://braingeneersdev/dolmo/original/outputs/validation/" 

#=============================================================================
# RUN COMMAND
#=============================================================================
export RUN_COMMAND="bash /opt/scripts/run_pipeline_local.sh"

#=============================================================================
# RESOURCE LIMITS & REQUESTS (Reduced for lightweight DE analysis)
#=============================================================================
export MEMORY_LIMIT="64Gi"
export CPU_LIMIT="8"
export STORAGE_LIMIT="100Gi"

export MEMORY_REQUEST="64Gi"
export CPU_REQUEST="8"
export STORAGE_REQUEST="100Gi"

#=============================================================================
# VALIDATION
#=============================================================================
if [ -z "$AWS_ACCESS_KEY_ID" ] || [ -z "$AWS_SECRET_ACCESS_KEY" ]; then
    echo "ERROR: AWS credentials not set!"
    exit 1
fi

#=============================================================================
# DEPLOY
#=============================================================================
echo "Deploying Kubernetes job: ${JOB_PREFIX}-${NAME}"
envsubst < jobdefinition.yaml | kubectl apply -f -
echo "Job created successfully!"