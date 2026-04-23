#!/usr/bin/env bash
# Custom Boltz2 wrapper: adds --no_kernels to work around Boltz 2.2.1
# torch-cu130 vs cuequivariance-ops-cu12 ABI mismatch.

set -euo pipefail

INDIR=$1
OUTDIR=$2

: "${BOLTZ2_DIFFUSION_SAMPLES:=5}"

PRED_TMPDIR="$(mktemp -d)"
trap 'rm -rf "$PRED_TMPDIR"' EXIT

if [[ -n "${BOLTZ2_BIN:-}" ]]; then
    "$BOLTZ2_BIN" predict "$INDIR" --out_dir "$PRED_TMPDIR" \
        --diffusion_samples "$BOLTZ2_DIFFUSION_SAMPLES" --no_kernels
else
    boltz predict "$INDIR" --out_dir "$PRED_TMPDIR" \
        --diffusion_samples "$BOLTZ2_DIFFUSION_SAMPLES" --no_kernels
fi

mkdir -p "$OUTDIR"
cp -a "$PRED_TMPDIR"/. "$OUTDIR"/
