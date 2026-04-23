#!/usr/bin/env bash
# Custom AF3 wrapper: passes --model_dir and --run_data_pipeline false
# so AF3 uses pre-computed MSAs (unpairedMsaPath) and finds model params.
#
# Required env vars:
#   AF3_CONTAINER, AF3_DB_PATH, AF3_ROOT, AF3_MODEL_DIR

set -euo pipefail

mkdir -p "$1" "$2"
INDIR=$(readlink -f "$1")
OUTDIR=$(readlink -f "$2")

: "${AF3_CONTAINER:?AF3_CONTAINER must be set}"
: "${AF3_DB_PATH:?AF3_DB_PATH must be set}"
: "${AF3_ROOT:?AF3_ROOT must be set}"
: "${AF3_MODEL_DIR:?AF3_MODEL_DIR must be set (directory containing af3.bin)}"

COPAL_ROOT="$(cd "$(dirname "$(readlink -f "$0")")/../.." && pwd)"
: "${AF3_WRAPPER:=${COPAL_ROOT}/third_party/alphafold3-patches/run_af3_with_host_confidence.py}"
: "${COPAL_AF3_HOST_CONFIDENCE_TYPES:=${COPAL_ROOT}/third_party/alphafold3-patches/confidence_types.py}"

export AF3_ROOT COPAL_AF3_HOST_CONFIDENCE_TYPES

PRED_TMPDIR="$(mktemp -d)"
trap 'rm -rf "$PRED_TMPDIR"' EXIT

unset APPTAINERENV_PYTHONPATH
apptainer exec --nv \
    --bind "$AF3_ROOT:$AF3_ROOT" \
    --bind "$AF3_MODEL_DIR:$AF3_MODEL_DIR" \
    --bind "$AF3_DB_PATH:$AF3_DB_PATH" \
    --bind "$COPAL_ROOT:$COPAL_ROOT" \
    --bind "$INDIR:$INDIR" \
    --bind "$PRED_TMPDIR:$PRED_TMPDIR" \
    "$AF3_CONTAINER" python "$AF3_WRAPPER" \
        --input_dir "$INDIR" \
        --run_data_pipeline false \
        --output_dir "$PRED_TMPDIR" \
        --db_dir "$AF3_DB_PATH" \
        --model_dir "$AF3_MODEL_DIR"

cp -a "$PRED_TMPDIR"/. "$OUTDIR"/
