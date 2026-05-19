#!/usr/bin/env bash
# Quickstart: 2 toy LuxR-family receptors × 3 AHL ligands, end-to-end.
#
#   Phase 1  copal predict       schedule + run the 6 predictions
#   Phase 2  copal organize      gather raw outputs into a hierarchical store
#   Phase 3  collect+postprocess build per_model_rankings.csv (legacy entrypoint)
#   Phase 4  copal full-ranking  emit stylized label-free ranking workbooks
#
# Phase 1 needs AF3 / Boltz2 / RF3 installed (see ../../INSTALLATION.md).
# Phases 2-4 are CPU-only. This toy set has no curated cognate-AHL ground
# truth, so the terminal step is `copal full-ranking set=inference` (the
# label-free path) rather than `copal analyze`. For the benchmark path with
# ground-truth labels, see ../../docs/reproduce_paper.md.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/out}"
export COPAL_REPO_ROOT="$REPO_ROOT"

# --- Phase 1: schedule + run inference -------------------------------------
# --executor local runs the predictions inline. On a SLURM cluster, drop
# `--executor local` and submit the emitted array instead:
#     sbatch "$OUT_DIR/sarray/copal.sbatch"
uv run copal predict \
    --receptors "$SCRIPT_DIR/receptors.csv" \
    --ligands   "$SCRIPT_DIR/ligands.csv" \
    --out-dir   "$OUT_DIR" \
    --models    af3,boltz2,rf3 \
    --batch-size 6 \
    --executor  local \
    --skip-done-check

# --- Phase 2: organize raw outputs -----------------------------------------
uv run copal organize --pred-dir "$OUT_DIR/outs" --dst-dir "$OUT_DIR/store"

# --- Phase 3: collect + postprocess ----------------------------------------
# The Hydra `copal analyze` / `full-ranking` commands consume
# per_model_rankings.csv; building it from raw outputs (collect + postprocess)
# currently lives in the legacy entrypoint. With a non-existent --labels-path
# it runs label-free and just emits the ranking CSVs.
uv run python -m copal.analysis.pipeline \
    --basepath     "$OUT_DIR/store/all" \
    --workflow-dir "$OUT_DIR/workflow" \
    --output-root  "$OUT_DIR/workflow" \
    --labels-path  "$OUT_DIR/no_labels.csv" \
    --ligands-path "$SCRIPT_DIR/ligands.csv" \
    --expected-ligand-count 3 \
    --skip-generate

# --- Phase 4: label-free ranking workbooks ---------------------------------
uv run copal full-ranking \
    set=inference \
    paths.workflow_dir="$OUT_DIR/workflow" \
    paths.output_root="$OUT_DIR/results"

echo
echo "Done. Stylized ranking workbooks:"
echo "  $OUT_DIR/results/inference/full_rankings/"
