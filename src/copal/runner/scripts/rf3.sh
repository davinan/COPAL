#!/usr/bin/env bash
# RosettaFold3 (Foundry) per-chunk wrapper (TEMPLATE — copy + edit for your cluster).
#
# Called as: rf3.sh <INPUT_DIR> <OUTPUT_DIR>
#
# Required env vars:
#   RF3_INSTALL  path to your Foundry/RF3 install (must contain .venv/bin/activate)
#   RF3_CKPT     path to the RF3 weights checkpoint
#
# Optional env vars:
#   RF3_EXTRA_ARGS   additional `rf3 fold` args, appended to the command line
#
# Override this script entirely by setting COPAL_RF3_SCRIPT=/your/path/rf3.sh
# before calling `copal predict`.

set -euo pipefail

INDIR=$1
OUTDIR=$2

: "${RF3_INSTALL:?RF3_INSTALL must point at your Foundry install (containing .venv/bin/activate)}"
: "${RF3_CKPT:?RF3_CKPT must point at the RF3 weights file}"

PRED_TMPDIR="$(mktemp -d)"
trap 'rm -rf "$PRED_TMPDIR"' EXIT

# shellcheck disable=SC1091
source "$RF3_INSTALL/.venv/bin/activate"
rf3 fold inputs="$INDIR" out_dir="$PRED_TMPDIR" ckpt_path="$RF3_CKPT" \
    early_stopping_plddt_threshold=0.0 ${RF3_EXTRA_ARGS:-}
deactivate

mkdir -p "$OUTDIR"
cp -a "$PRED_TMPDIR"/. "$OUTDIR"/
