#!/bin/bash
# Deploy Kubernetes job for pipeline
# All configuration is done via environment variables below

#=============================================================================
# JOB IDENTIFICATION
#=============================================================================
export NAME="rima-job-gpu-moreepochs"                    # Job name identifier
export JOB_PREFIX="dolmo"                 # Prefix for job naming
export CONTAINER_NAME="rima-image"        # Container name

#=============================================================================
# DOCKER IMAGE
#=============================================================================
export DOCKER_IMAGE="dolmomar/rima-public:v6"
#=============================================================================
# AWS/S3 CREDENTIALS (required - set these before running)
#=============================================================================
export AWS_ACCESS_KEY_ID="${AWS_ACCESS_KEY_ID:-}"   
export AWS_SECRET_ACCESS_KEY="${AWS_SECRET_ACCESS_KEY:-}"

# UPDATED: Pointing to the external NRP S3 endpoint you verified
export S3_ENDPOINT="https://s3-west.nrp-nautilus.io"

#=============================================================================
# PATHS
#=============================================================================
export PROJECT_ROOT="/workspace"
export OUTPUT_DIR="/workspace/outputs"

# NEW: Your exact S3 input bucket
export INPUT_S3_PATH="s3://braingeneersdev/dolmo/original/data/"

# NEW: A suggested S3 output bucket for your results
export S3_ENDPOINT_URL_RESULTS="s3://braingeneersdev/dolmo/original/outputs/new_metadata/high_epochs/"

#=============================================================================
# RUN COMMAND
#=============================================================================
export RUN_COMMAND="bash /opt/scripts/run_pipeline_local.sh"

#=============================================================================
# RESOURCE LIMITS
#=============================================================================
export MEMORY_LIMIT="400Gi"
export CPU_LIMIT="20"
export STORAGE_LIMIT="400Gi"

#=============================================================================
# RESOURCE REQUESTS
#=============================================================================
export MEMORY_REQUEST="400Gi"
export CPU_REQUEST="20"
export STORAGE_REQUEST="400Gi"

#=============================================================================
# VALIDATION
#=============================================================================
if [ -z "$AWS_ACCESS_KEY_ID" ] || [ -z "$AWS_SECRET_ACCESS_KEY" ]; then
    echo "ERROR: AWS credentials not set!"
    echo "Please set AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY environment variables"
    exit 1
fi

#=============================================================================
# DEPLOY
#=============================================================================
echo "Deploying Kubernetes job: ${JOB_PREFIX}-${NAME}-${CONTAINER_NAME}"

# Apply the job definition with envsubst
envsubst < jobdefinition.yaml | kubectl apply -f -

echo "Job created successfully!"
echo ""
echo "Monitor with:"
echo "  kubectl get jobs | grep ${JOB_PREFIX}-${NAME}"
echo "  kubectl get pods | grep ${JOB_PREFIX}-${NAME}"
echo "  kubectl logs -f <pod-name>"
