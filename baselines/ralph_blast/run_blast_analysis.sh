#!/usr/bin/env bash
# Run the RALPH BLAST baseline.
#
# Required env vars (or pass them via flags):
#   BENCHMARK_CSV     labels CSV (default: ../../data/paper_luxrs.csv)
#   RALPH_DATABASE    dir containing built BLAST DB + ralph_scores.csv
#                     (default: ./ralph_database; see README.md to build)
#   LIGANDS_PATH      annotated ligand CSV
#                     (default: ../../data/curated_ligands_annotated.csv)
#   OUTPUT_ROOT       results directory (default: ./outputs/results)
#   WORKFLOW_DIR      intermediate dir (default: ./workflow_ralph)
#   MIN_PIDENT        BLAST percent-identity floor (default: 30)
#   MAX_HITS          BLAST max_target_seqs (default: 50)
#
# Requires `blastp` on PATH and `pip install -e ../../[inference]` for the
# `copal` package, plus the `ralph` package (see README.md).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

BENCHMARK_CSV="${BENCHMARK_CSV:-$REPO_ROOT/data/paper_luxrs.csv}"
RALPH_DATABASE="${RALPH_DATABASE:-$SCRIPT_DIR/ralph_database}"
LIGANDS_PATH="${LIGANDS_PATH:-$REPO_ROOT/data/curated_ligands_annotated.csv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$SCRIPT_DIR/outputs/results}"
WORKFLOW_DIR="${WORKFLOW_DIR:-$SCRIPT_DIR/workflow_ralph}"
MIN_PIDENT="${MIN_PIDENT:-30}"
MAX_HITS="${MAX_HITS:-50}"

python "$SCRIPT_DIR/run_ralph_blast.py" \
    --benchmark-csv "$BENCHMARK_CSV" \
    --ralph-database "$RALPH_DATABASE" \
    --ligands-path "$LIGANDS_PATH" \
    --output-root "$OUTPUT_ROOT" \
    --workflow-dir "$WORKFLOW_DIR" \
    --min-pident "$MIN_PIDENT" \
    --max-hits "$MAX_HITS"
