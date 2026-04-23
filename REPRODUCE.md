# Reproducing the paper

Two paths, depending on whether you want to re-run inference.

## Path A — analysis only (no GPUs, ~minutes)

The fastest way to reproduce the figures and tables from the paper.

1. Download `predictions.csv` (~27 MB) from Zenodo (`<DOI_ANALYSIS>`,
   link will be added at publication time) and place it under `data/`:

   ```bash
   wget <ZENODO_URL>/predictions.csv -O data/predictions.csv
   ```

2. Install the package and run the analysis:

   ```bash
   uv sync
   uv run copal analyze \
       --workflow-dir outputs_official_workflow \
       --output-root  outputs_official_workflow/results \
       --labels       data/benchmark_data_all_corrected.csv \
       --ligands      data/Final_confirmed_ligands_58_annotated.csv \
       --input-csv    data/predictions.csv \
       --no-leakage-csv data/curated_conf1_no_leakage.csv
   ```

3. Outputs:

   ```
   outputs_official_workflow/results/
   ├── ranking_performance/   # top-K accuracy, Kendall τ per model + ensemble
   ├── full_rankings/         # complete rankings + per-model XLSX
   ├── ablation/              # model-combination ablation
   ├── chemical_similarity/   # RDKit MCS / Morgan similarity (requires rdkit)
   ├── acyl_moiety/           # acyl-chain match metrics
   └── mpkt/                  # mpKT chemoinformatic tables
   ```

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
       --receptors data/benchmark_data_all_corrected.csv \
       --ligands   data/benchmark_ligands.csv \
       --out-dir   /scratch/$USER/copal_benchmark \
       --models    af3,boltz2,rf3 \
       --batch-size 58 \
       --redo-schedule
   sbatch /scratch/$USER/copal_benchmark/sarray/copal.sbatch
   ```

   The full benchmark is ~12k predictions × 3 models. Plan accordingly.

4. After the array finishes, organize then analyze (this re-creates the
   `predictions.csv` you'd otherwise download):

   ```bash
   uv run copal organize \
       --pred-dir  /scratch/$USER/copal_benchmark/outs \
       --dst-dir   /scratch/$USER/copal_benchmark/store \
       --models    af3,boltz2,rf3

   uv run copal analyze \
       --workflow-dir outputs_official_workflow \
       --output-root  outputs_official_workflow/results \
       --labels       data/benchmark_data_all_corrected.csv \
       --ligands      data/Final_confirmed_ligands_58_annotated.csv \
       --basepath     /scratch/$USER/copal_benchmark/store/all \
       --af3-iptm-source json \
       --no-leakage-csv data/curated_conf1_no_leakage.csv
   ```

## Re-running just the leakage filter

```bash
uv run copal seqid-filter \
    --queries    data/benchmark_data_all_corrected.csv \
    --references path/to/pdb_receptors.csv \
    --out-table  outputs/seqid_table.csv \
    --out-leakage data/curated_conf1_no_leakage.csv \
    --threshold  0.40
```

The references CSV needs `refahl_id` (or `receptor_id`), `seq`, and
`cocrystal_ligand` columns to enable the same-ligand check.

## Reproducing the RALPH BLAST baseline

See [`baselines/ralph_blast/README.md`](baselines/ralph_blast/README.md)
for BLAST+ install, RALPH database build, and the one-line invocation.
