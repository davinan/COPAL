# Benchmark example

This directory is a thin pointer for reproducing the paper benchmark.
The full reproduction recipe lives in [REPRODUCE.md](../../REPRODUCE.md)
at the repo root.

In short:

1. Download `predictions.csv` from the public data archive (DOI in the
   repo root README) into `data/`.
2. From the repo root:

   ```bash
   uv run copal analyze \
       --workflow-dir outputs_official_workflow \
       --output-root  outputs_official_workflow/results \
       --labels       data/benchmark_data_all_corrected.csv \
       --ligands      data/Final_confirmed_ligands_58_annotated.csv \
       --input-csv    data/predictions.csv \
       --no-leakage-csv data/curated_conf1_no_leakage.csv
   ```
