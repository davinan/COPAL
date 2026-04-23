# Using COPAL on your own data

COPAL is designed so that running predictions on a new set of receptors
and ligands is a one-line change: swap two CSVs.

## 1. Prepare the input CSVs

### `receptors.csv`

Required columns:

| Column | Type | Notes |
| ------ | ---- | ----- |
| `receptor_id` | string | Unique ID. |
| `seq` | string | Amino-acid sequence (single chain). |

Optional columns:

| Column | Notes |
| ------ | ----- |
| `msa_path` | Path to a pre-built MSA (`.a3m`). If absent, the model runner generates one at inference time. |
| `chain_count` | Integer; for homo-oligomers (default 1). |

### `ligands.csv`

Required columns:

| Column | Type | Notes |
| ------ | ---- | ----- |
| `ligand_id` | string | Unique ID. |
| `smiles` | string | SMILES. |

Optional metadata columns (`inchi`, `name`, `formula`, …) are passed
through but unused by inference.

### Backwards-compatible aliases

`copal` will auto-rename a few legacy column names so the paper benchmark
files work as-is:

| Receptor alias | → |
| -------------- | -- |
| `LuxR_name`, `refahl_id`, `target` | `receptor_id` |
| `sequence` | `seq` |

| Ligand alias | → |
| ------------ | -- |
| `AHLid`, `ahl_id`, `abbreviation` | `ligand_id` |
| `Isomeric SMILES`, `smi`, `SMILES` | `smiles` |

See [`docs/input_formats.md`](docs/input_formats.md) for full schema details.

## 2. Run inference

COPAL supports two execution backends — pick with `--executor`:

### SLURM cluster (default)

```bash
uv run copal predict \
    --receptors my_proteins.csv \
    --ligands   my_ligands.csv \
    --out-dir   ./run \
    --models    af3,boltz2,rf3 \
    --batch-size 50
```

Outputs:

```
run/
├── scheduled_ensemble.csv      # one row per (receptor, ligand) pair
├── inputs/{AF3,Boltz2,RF3}/    # per-model JSON / YAML
├── sarray/copal.{cmds,sbatch}  # SLURM array script (run with `sbatch`)
├── submit_sarray.sh            # convenience launcher
└── outs/                       # populated as predictions complete
```

Add `--submit` to auto-submit the array, or run
`bash run/submit_sarray.sh` manually.

### Single machine / no SLURM

```bash
uv run copal predict \
    --receptors my_proteins.csv \
    --ligands   my_ligands.csv \
    --out-dir   ./run \
    --models    af3,boltz2,rf3 \
    --batch-size 50 \
    --executor  local \
    --workers   1
```

The `.cmds` file is still written (handy for inspection and manual
re-runs), but instead of a SLURM script COPAL executes each command on
the current shell. Per-task stdout/stderr stream to
`run/logs/task-NNN.{out,err}`.

Use `--workers N` to run up to N commands concurrently. This is usually
only worth it if you have multiple GPUs or models that can share a GPU
(AF3 and RF3 both want a full GPU; Boltz2 is lighter).

Use `--dry-run` (either executor) to write the schedule + inputs +
`.cmds` file and stop without executing anything.

> The pipeline is **resumable**: re-running `copal predict` with the same
> `--out-dir` checks `outs/` for existing predictions and only emits
> work for the missing ones.

## 3. Rank the predictions

### Prediction-only mode (no ground truth)

```bash
uv run copal organize --pred-dir ./run/outs --dst-dir ./run/store
uv run copal analyze \
    --workflow-dir ./run/workflow \
    --output-root  ./run/results \
    --labels       my_proteins.csv \
    --ligands      my_ligands.csv \
    --basepath     ./run/store/all
```

You'll get per-receptor ranked ligand lists (`ranking_performance/…`,
`full_rankings/…`) computed from each model's confidence scores plus an
ensemble.

### Benchmark mode (with ground-truth labels)

If your `labels` CSV has a `true_ahl_abbrev` column (or your equivalent
"correct ligand" column), `analyze` additionally computes top-K accuracy,
Kendall τ / Spearman ρ, ablation tables, and chemical-similarity stats —
exactly as in the paper.

To replicate the paper benchmark numbers, see
[`REPRODUCE.md`](REPRODUCE.md).

## 4. Optional: filter by sequence similarity

To exclude query sequences that "leak" from the test set into a known
reference set (e.g. PDB structures), use:

```bash
uv run copal seqid-filter \
    --queries     queries.csv \
    --references  references.csv \
    --out-table   pairwise.csv \
    --out-leakage filtered.csv \
    --threshold   0.40
```

The reference CSV may carry a `cocrystal_ligand` column; if both the
queries CSV (`true_ahl_abbrev`) and the references CSV have ligand
annotations, "leakage" is defined as
`(identity ≥ threshold) AND (same cognate ligand)`.

## 5. Common knobs

| Flag | Default | Effect |
| ---- | ------- | ------ |
| `--models af3,boltz2` | `af3,boltz2,rf3` | Run a subset. |
| `--batch-size N` | 20 | Predictions per chunk (= per SLURM task, or per local command). |
| `--executor {slurm,local}` | `slurm` | Dispatch mode. |
| `--workers N` | 1 | Concurrent tasks when `--executor local`. |
| `--dry-run` | off | Write inputs + `.cmds` but don't execute. |
| `--submit` | off | Auto-submit the sbatch (slurm mode only). |
| `--skip-done-check` | off | Skip the `outs/` scan (faster, but doesn't resume). |
| `--redo-schedule` | off | Force-rebuild `scheduled_ensemble.csv`. |
| `--done-lookup PATH` | none | JSON `{ "RECEPTOR_LIGAND": ["af3", "boltz2"] }` of pairs to mark already-done. |
