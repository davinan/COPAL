# Benchmark example

This directory is a thin pointer for reproducing the paper benchmark.
The canonical recipe lives in
[`docs/reproduce_paper.md`](../../docs/reproduce_paper.md); the two
higher-level reproduction paths are in
[`REPRODUCE.md`](../../REPRODUCE.md) at the repo root.

In short:

1. Download `predictions.csv` from the public data archive (DOI in the
   repo root README) into `data/`.
2. From the repo root, build the postprocessed ranking CSVs, then run the
   Hydra-driven analysis on the published `final_n96` set:

   ```bash
   export COPAL_REPO_ROOT=$PWD

   # predictions.csv is the raw collected table -> run postprocess
   uv run copal analyze \
       set=final_n96 \
       analysis.skip.postprocess=false \
       paths.input_csv=data/predictions.csv

   # similarity matrix + the same analysis with the paper's mcs_penalized metric
   uv run copal sim-matrix sim_matrix.method=mcs_penalized
   uv run copal analyze \
       set=final_n96 \
       analysis.chem_sim_matrix=data/chem_sim_matrices/mcs_penalized.csv
   ```

See [`docs/reproduce_paper.md`](../../docs/reproduce_paper.md) for the
stylized-workbook step, the RRF-aggregation variant, and the case-study
set.
