#!/usr/bin/env bash
# AF3 wrapper — native (non-container) install.
# Uses the uv venv at $AF3_ROOT/.venv and calls run_alphafold.py directly.
#
# Required env vars:
#   AF3_ROOT       AF3 source clone (contains .venv/, run_alphafold.py)
#   AF3_DB_PATH    sequence DB dir (unused when --run_data_pipeline false)
#   AF3_MODEL_DIR  dir containing af3.bin / af3.bin.zst
#
# Called as: af3_native.sh <INPUT_DIR> <OUTPUT_DIR>

set -euo pipefail

mkdir -p "$1" "$2"
INDIR=$(readlink -f "$1")
OUTDIR=$(readlink -f "$2")

: "${AF3_ROOT:?AF3_ROOT must be set}"
: "${AF3_DB_PATH:?AF3_DB_PATH must be set}"
: "${AF3_MODEL_DIR:?AF3_MODEL_DIR must be set}"

PRED_TMPDIR="$(mktemp -d)"
trap 'rm -rf "$PRED_TMPDIR"' EXIT

# shellcheck disable=SC1091
source "$AF3_ROOT/.venv/bin/activate"

python "$AF3_ROOT/run_alphafold.py" \
    --input_dir "$INDIR" \
    --run_data_pipeline false \
    --output_dir "$PRED_TMPDIR" \
    --db_dir "$AF3_DB_PATH" \
    --model_dir "$AF3_MODEL_DIR"

cp -a "$PRED_TMPDIR"/. "$OUTDIR"/
