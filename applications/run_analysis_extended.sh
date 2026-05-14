#!/usr/bin/env bash
# Reproduce the COPAL paper analysis from the bundled predictions table,
# emitting the *extended* layout (per-conf no-leakage variants for conf
# {1,2,3}; mpKT pooled across conf 1+2+3; ablation on conf 1+2+3).
#
# Steps:
#   1. Decompress data/predictions.csv.xz -> data/predictions.csv (if needed).
#   2. `copal analyze --skip-generate` to build per_model_rankings.csv,
#      ensemble_rankings.csv, ablation_table.csv, done_pairs.json under
#      $WORKFLOW_DIR.
#   3. `applications/copal_analysis_extended.py` to write the extended
#      generate-phase outputs under $OUTPUT_ROOT.
#
# Overrides (env vars):
#   PREDICTIONS_CSV    default: data/predictions.csv
#   WORKFLOW_DIR       default: outputs_official_workflow
#   OUTPUT_ROOT        default: $WORKFLOW_DIR/results
#   EXCLUDE_SUFFIXES   space-separated refahl_id suffixes to drop from the
#                      no-leakage subsets before the extended generate phase
#                      (default: "_R1 _R2"; use "" to keep all curated rows).
#   AGGREGATION        ensemble aggregator: sum (default) | mean | median |
#                      rrf | min | max. Threaded into both `copal analyze`
#                      and the extended generate phase.
#   RRF_K              RRF damping constant (default 60; only used when
#                      AGGREGATION=rrf).
#
# Sweep example (one run per aggregator into separate result dirs):
#   for AGG in rrf sum mean median min max; do
#       AGGREGATION=$AGG \
#       OUTPUT_ROOT="outputs_official_workflow/results_$AGG" \
#         bash applications/run_analysis_extended.sh
#   done
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

PREDICTIONS_CSV="${PREDICTIONS_CSV:-data/predictions.csv}"
WORKFLOW_DIR="${WORKFLOW_DIR:-outputs_official_workflow}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$WORKFLOW_DIR/results}"
EXCLUDE_SUFFIXES="${EXCLUDE_SUFFIXES:-_R1 _R2}"
AGGREGATION="${AGGREGATION:-sum}"
RRF_K="${RRF_K:-60}"

if [[ ! -f "$PREDICTIONS_CSV" ]]; then
    if [[ -f "${PREDICTIONS_CSV}.xz" ]]; then
        echo "Decompressing ${PREDICTIONS_CSV}.xz -> ${PREDICTIONS_CSV}"
        xz -dk "${PREDICTIONS_CSV}.xz"
    else
        echo "ERROR: neither $PREDICTIONS_CSV nor ${PREDICTIONS_CSV}.xz exists." >&2
        exit 1
    fi
fi

echo "[1/2] Stock postprocess (rankings, ensemble, ablation, done_pairs)"
uv run copal analyze \
    --workflow-dir   "$WORKFLOW_DIR" \
    --output-root    "$OUTPUT_ROOT" \
    --labels         data/benchmark_data_all_corrected.csv \
    --ligands        data/Final_confirmed_ligands_58_annotated.csv \
    --input-csv      "$PREDICTIONS_CSV" \
    --no-leakage-csv data/curated_conf1_no_leakage.csv \
    --aggregation    "$AGGREGATION" \
    --rrf-k          "$RRF_K" \
    --skip-generate

FILTERED_DIR="$WORKFLOW_DIR/_no_leakage_filtered"
mkdir -p "$FILTERED_DIR"
filter_curated() {
    local src="$1" dst="$2"
    if [[ -z "$EXCLUDE_SUFFIXES" ]]; then
        cp "$src" "$dst"
        return
    fi
    EXCLUDE_SUFFIXES="$EXCLUDE_SUFFIXES" python -c '
import os, sys
suffixes = tuple(s for s in os.environ["EXCLUDE_SUFFIXES"].split() if s)
src, dst = sys.argv[1], sys.argv[2]
dropped = []
with open(src) as fh, open(dst, "w") as out:
    out.write(fh.readline())
    for line in fh:
        rid = line.split(",", 1)[0]
        if rid.endswith(suffixes):
            dropped.append(rid)
        else:
            out.write(line)
tail = " -> " + ", ".join(dropped) if dropped else ""
print("  " + os.path.basename(src) + ": dropped " + str(len(dropped)) + tail)
' "$src" "$dst"
}

echo "Filtering no-leakage CSVs (excluding suffixes: ${EXCLUDE_SUFFIXES:-<none>})"
filter_curated data/curated_conf1_no_leakage.csv "$FILTERED_DIR/curated_conf1_no_leakage.csv"
filter_curated data/curated_conf2_no_leakage.csv "$FILTERED_DIR/curated_conf2_no_leakage.csv"
filter_curated data/curated_conf3_no_leakage.csv "$FILTERED_DIR/curated_conf3_no_leakage.csv"

echo "[2/2] Extended generate phase (aggregation=$AGGREGATION rrf_k=$RRF_K)"
uv run python applications/copal_analysis_extended.py \
    --workflow-dir     "$WORKFLOW_DIR" \
    --output-root      "$OUTPUT_ROOT" \
    --no-leakage-conf1 "$FILTERED_DIR/curated_conf1_no_leakage.csv" \
    --no-leakage-conf2 "$FILTERED_DIR/curated_conf2_no_leakage.csv" \
    --no-leakage-conf3 "$FILTERED_DIR/curated_conf3_no_leakage.csv" \
    --aggregation      "$AGGREGATION" \
    --rrf-k            "$RRF_K"

echo
echo "Done. Extended results under: $OUTPUT_ROOT"
