#!/usr/bin/env bash
# Quickstart: schedule predictions for 2 toy LuxR-family receptors against 3 AHL ligands.
#
# Prerequisites:
#   - `uv sync` from the repo root
#   - AF3, Boltz2, RF3 set up per ../../INSTALLATION.md (env vars for wrappers)
#
# This emits the SLURM array script + per-model inputs without submitting.
# Add `--submit` to the copal call to actually queue the array.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/out}"

uv run copal predict \
    --receptors "$SCRIPT_DIR/receptors.csv" \
    --ligands "$SCRIPT_DIR/ligands.csv" \
    --out-dir "$OUT_DIR" \
    --models af3,boltz2,rf3 \
    --batch-size 6 \
    --skip-done-check
