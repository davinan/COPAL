#!/usr/bin/env bash
# Boltz2 per-chunk wrapper (TEMPLATE — copy + edit for your cluster).
#
# Called as: boltz2.sh <INPUT_DIR> <OUTPUT_DIR>
#
# Either BOLTZ2_BIN (absolute path to a `boltz` or `run_boltz.sh` launcher)
# must be set, or `boltz` must be on PATH.
#
# Optional env vars:
#   BOLTZ2_DIFFUSION_SAMPLES   number of samples per input (default: 5)
#
# Override this script entirely by setting COPAL_BOLTZ2_SCRIPT=/your/path/boltz2.sh
# before calling `copal predict`.

set -euo pipefail

INDIR=$1
OUTDIR=$2

: "${BOLTZ2_DIFFUSION_SAMPLES:=5}"

PRED_TMPDIR="$(mktemp -d)"
trap 'rm -rf "$PRED_TMPDIR"' EXIT

if [[ -n "${BOLTZ2_BIN:-}" ]]; then
    "$BOLTZ2_BIN" predict "$INDIR" --out_dir "$PRED_TMPDIR" --diffusion_samples "$BOLTZ2_DIFFUSION_SAMPLES"
else
    boltz predict "$INDIR" --out_dir "$PRED_TMPDIR" --diffusion_samples "$BOLTZ2_DIFFUSION_SAMPLES"
fi

mkdir -p "$OUTDIR"
cp -a "$PRED_TMPDIR"/. "$OUTDIR"/
