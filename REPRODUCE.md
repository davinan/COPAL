# Reproducing the paper

Two paths, depending on whether you want to re-run inference. For the exact
canonical command sequence (similarity matrix → stylized workbooks →
analysis artifacts), see [`docs/reproduce_paper.md`](docs/reproduce_paper.md);
this file is the higher-level overview.

`copal analyze` is **Hydra-driven** — it takes `key=value` config overrides,
not `--flags`. The benchmark set, I/O paths, and aggregation are all config
leaves. See [`USAGE.md`](USAGE.md) §3–4.

## Path A — analysis only (no GPUs, ~minutes)

The fastest way to reproduce the figures and tables from the paper.

1. Download `predictions.csv` (~27 MB) from the public data archive
   (host TBD; `<DOI_ANALYSIS>`, link will be added at publication time)
   and place it under `data/`:

   ```bash
   wget <DATA_URL>/predictions.csv -O data/predictions.csv
   ```

   `data/predictions.csv` is the pre-collected raw prediction table; it is
   not shipped in the repo (it is git-ignored).

2. Install the package and run the analysis. `predictions.csv` is the raw
   collected table, so postprocess must run (it builds
   `per_model_rankings.csv` + `ensemble_rankings.csv` under
   `paths.workflow_dir`):

   ```bash
   uv sync
   uv run copal analyze \
       set=final_n96 \
       analysis.skip.postprocess=false \
       paths.input_csv=data/predictions.csv
   ```

   Already have `per_model_rankings.csv` in `outputs_official_workflow/`
   from an earlier run? Drop the two overrides and just run
   `uv run copal analyze set=final_n96` — postprocess is skipped by default.

3. Outputs land under `outputs_official_workflow/results/final_n96/` in
   seven artifact groups:

   ```
   results/final_n96/
   ├── ranking_performance/   # true-AHL rank per model + ensemble
   ├── mpkt/                  # cross-model agreement (mean pairwise Kendall τ)
   ├── acyl_moiety/           # true-AHL rank stratified by acyl subcategory
   ├── chemical_similarity/   # RDKit MCS-based similarity to the true AHL
   ├── score_vs_rank/         # per-model confidence vs true-AHL rank
   ├── full_rankings/         # wide rank matrix (+ stylized .xlsx via full-ranking)
   └── ablation/              # model-combination ablation
   ```

   Each group's schema is documented in [`docs/artifacts.md`](docs/artifacts.md).

## Path B — full inference (GPU cluster, days–weeks)

Re-run AlphaFold 3, Boltz2, and RosettaFold3 on the benchmark from
scratch, then re-collect → analyze.

1. Set up AF3, Boltz2, and RF3 per [`INSTALLATION.md`](INSTALLATION.md).
   Verify the env vars the bundled wrappers read:

   ```bash
   echo $AF3_CONTAINER $AF3_DB_PATH $AF3_ROOT $BOLTZ2_BIN $RF3_INSTALL $RF3_CKPT
   ```

2. Extract the bundled MSA archive (shipped in-tree, 31 MB):

   ```bash
   mkdir -p data/benchmark_msas
   tar -xJf data/benchmark_msas.tar.xz -C data/benchmark_msas/
   ```

3. Schedule and submit:

   ```bash
   uv run copal predict \
       --receptors data/paper_luxrs.csv \
       --ligands   data/curated_ligands_annotated.csv \
       --out-dir   /scratch/$USER/copal_benchmark \
       --models    af3,boltz2,rf3 \
       --batch-size 58 \
       --redo-schedule
   sbatch /scratch/$USER/copal_benchmark/sarray/copal.sbatch
   ```

   The full benchmark is ~12k predictions × 3 models. Plan accordingly.

4. After the array finishes, organize the raw outputs:

   ```bash
   uv run copal organize \
       --pred-dir  /scratch/$USER/copal_benchmark/outs \
       --dst-dir   /scratch/$USER/copal_benchmark/store \
       --models    af3,boltz2,rf3
   ```

5. **Collect** the organized outputs into `predictions.csv` and postprocess
   them into `per_model_rankings.csv` + `ensemble_rankings.csv`. The
   Hydra-driven `copal analyze` runs postprocess + analysis but **not**
   raw-output collection — collection currently lives only in the legacy
   (deprecated) entrypoint:

   ```bash
   uv run python -m copal.analysis.pipeline \
       --basepath     /scratch/$USER/copal_benchmark/store/all \
       --workflow-dir outputs_official_workflow \
       --output-root  outputs_official_workflow/results \
       --labels-path  data/paper_luxrs.csv \
       --ligands-path data/curated_ligands_annotated.csv \
       --skip-generate
   ```

   `--skip-generate` stops after postprocess; ipTM source defaults to
   `json` (override with `--af3-iptm-source base` for the base-folder CSV).

6. Run the Hydra analysis (postprocess is now done, so it is skipped):

   ```bash
   uv run copal analyze set=final_n96
   ```

## Re-running just the leakage filter

```bash
uv run copal seqid-filter \
    --queries    data/paper_luxrs.csv \
    --references path/to/pdb_receptors.csv \
    --out-table  outputs/seqid_table.csv \
    --out-leakage data/curated_luxR_RefAHL_no_leakage.csv \
    --threshold  0.40
```

The references CSV needs `refahl_id` (or `receptor_id`), `seq`, and
`cocrystal_ligand` columns to enable the same-ligand check.

## Reproducing the RALPH BLAST baseline

See [`baselines/ralph_blast/README.md`](baselines/ralph_blast/README.md)
for BLAST+ install, RALPH database build, and the one-line invocation.
