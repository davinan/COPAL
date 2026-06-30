# Sequence-identity "nearest-leak" baseline

A sequence-only baseline for the COPAL AHL-ranking benchmark, analogous to the RALPH
BLAST baseline but with no external binaries. It builds a sequence → AHL database from
the **leaked** receptors (those excluded from the curated no-leakage set) and, for each
**leak-free** benchmark receptor, finds the closest leaked sequence by global-alignment
% identity (BLOSUM62 — the same metric that defined the no-leakage set), takes that
receptor's cognate AHL as an anchor, and ranks every candidate AHL by its **MCS
bond-penalized** chemical similarity to the anchor (the precomputed `mcs_penalized`
matrix). It then writes the standard COPAL benchmarking CSVs by reusing
`copal.analysis.pipeline`.

## Run

```bash
uv run python baselines/seqid_nn/run_seqid_baseline.py \
    --labels data/paper_luxrs.csv \
    --no-leakage-csv data/curated_luxR_RefAHL_no_leakage.csv \
    --ligands-path data/curated_ligands_annotated.csv \
    --sim-matrix data/chem_sim_matrices/mcs_penalized.csv \
    --output-root outputs_official_workflow/results/seqid_no_leakage \
    --db-scope refahl     # or: all
```

`--db-scope refahl` (default) uses only the RefAHL-split leaked receptors as the
database (strict "leaked vs training PDB" examples); `--db-scope all` uses every
non-no-leakage receptor with a usable in-pool AHL (a larger nearest-neighbour pool).

## Outputs

- `<workflow-dir>/per_model_rankings.csv`, `ensemble_rankings.csv`
  (default `<output-root>/../workflow_seqid/`)
- `<output-root>/ranking_performance/ranking_performance_conf{1,2,3}.csv`
- `<output-root>/full_rankings/all_splits.csv` and `SeqID.xlsx`
