#!/usr/bin/env bash
# AlphaFold 3 per-chunk wrapper (TEMPLATE — copy + edit for your cluster).
#
# Called as: af3.sh <INPUT_DIR> <OUTPUT_DIR>
#   <INPUT_DIR>  contains one or more AF3 *.json files for this chunk
#   <OUTPUT_DIR> will receive the per-prediction result directories
#
# Required env vars:
#   AF3_CONTAINER  path to the AF3 Apptainer/Singularity image (.sif)
#   AF3_DB_PATH    path to the AF3 sequence database directory
#   AF3_ROOT       path to the AF3 source clone (needed by the 7-decimal wrapper)
#
# Optional env vars:
#   COPAL_AF3_HOST_CONFIDENCE_TYPES
#       Path to the patched confidence_types.py. If unset, the default
#       bundled copy in third_party/alphafold3-patches/ is used.
#   AF3_WRAPPER
#       Path to the runtime-injection wrapper. Defaults to the shipped
#       third_party/alphafold3-patches/run_af3_with_host_confidence.py.
#
# Override this script entirely by setting COPAL_AF3_SCRIPT=/your/path/af3.sh
# before calling `copal predict`.

set -euo pipefail

INDIR=$1
OUTDIR=$2

: "${AF3_CONTAINER:?AF3_CONTAINER must point at the AF3 apptainer/singularity image}"
: "${AF3_DB_PATH:?AF3_DB_PATH must point at the AF3 sequence database}"
: "${AF3_ROOT:?AF3_ROOT must point at your alphafold3 source clone}"

# Locate the bundled AF3 patch wrapper, two levels up from this script:
#   src/copal/runner/scripts/af3.sh  →  third_party/alphafold3-patches/
COPAL_ROOT="$(cd "$(dirname "$(readlink -f "$0")")/../../../.." && pwd)"
: "${AF3_WRAPPER:=${COPAL_ROOT}/third_party/alphafold3-patches/run_af3_with_host_confidence.py}"
: "${COPAL_AF3_HOST_CONFIDENCE_TYPES:=${COPAL_ROOT}/third_party/alphafold3-patches/confidence_types.py}"

export AF3_ROOT COPAL_AF3_HOST_CONFIDENCE_TYPES

PRED_TMPDIR="$(mktemp -d)"
trap 'rm -rf "$PRED_TMPDIR"' EXIT

unset APPTAINERENV_PYTHONPATH
apptainer exec --nv "$AF3_CONTAINER" python "$AF3_WRAPPER" \
    --input_dir "$INDIR" \
    --run_data_pipeline true \
    --output_dir "$PRED_TMPDIR" \
    --db_dir "$AF3_DB_PATH"

mkdir -p "$OUTDIR"
cp -a "$PRED_TMPDIR"/. "$OUTDIR"/
