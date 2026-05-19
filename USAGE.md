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

First collect the raw model outputs into a hierarchical store:

```bash
uv run copal organize --pred-dir ./run/outs --dst-dir ./run/store
```

Ranking is then done with one of two subcommands, **both Hydra-driven** —
you select a benchmark *set* and override config leaves with `key=value`
syntax. There are no `--flags`.

### `copal full-ranking` — label-free ranking (your own data)

If you have no curated ground truth, `copal full-ranking set=inference`
ranks every ligand per receptor using each model's confidence score and an
ensemble, and emits stylized `.xlsx` workbooks:

```bash
uv run copal full-ranking set=inference
```

Output: `results/inference/full_rankings/{AF3,Boltz2,Boltz2Affinity,RF3,ENSEMBLE}.xlsx`
plus `full_rankings.csv` (the programmatic wide-form table). See Section 5.

### `copal analyze` — benchmark analysis (with ground-truth labels)

When the selected set's labels CSV carries cognate-ligand ground truth,
`copal analyze` computes top-K accuracy, Kendall τ, ablation tables,
chemical-similarity, and per-model score-vs-rank correlations:

```bash
uv run copal analyze set=refAHL_all_conf_no_leakage
```

`copal analyze` **requires** a non-empty labelled set — it exits with an
error on a set with no ground truth. For label-free data, use
`copal full-ranking set=inference` instead.

Artifacts land under `<paths.output_root>/<set.name>/` in seven groups:
`ranking_performance/`, `mpkt/`, `acyl_moiety/`, `chemical_similarity/`,
`score_vs_rank/`, `full_rankings/`, `ablation/`. Each is documented in
[`docs/artifacts.md`](docs/artifacts.md).

To reproduce the paper benchmark numbers, see
[`docs/reproduce_paper.md`](docs/reproduce_paper.md).

## 4. Benchmark sets

A *set* is a Hydra config group under `src/copal/config/set/*.yaml`. Each
set resolves to a filtered view of a labels CSV and declares `labels_path`,
`ligands_path`, and a `filter` block with five independent axes:

| Filter axis | Effect (empty / missing = no filter) |
| ----------- | ------------------------------------ |
| `splits` | Keep only rows in these split labels (case-insensitive). |
| `confidences` | Keep only these confidence tiers (1 / 2 / 3). |
| `no_leakage_csvs` | Keep only LuxRs present in the union of these curated CSVs. |
| `exclude_suffixes` | Drop LuxRs whose name ends with any of these (case-insensitive). |
| `include_names` | Keep **only** these LuxRs — a case-insensitive allowlist. |

Select a set with `set=<name>`; override one leaf with Hydra list syntax,
e.g. `'set.filter.confidences=[1,2]'`. Shipped sets:

| Set | What it is |
| --- | ---------- |
| `refAHL_all_conf_no_leakage` | Default. RefAHL, confidence {1,2,3}, no-leakage union, replicate/litter suffixes dropped. |
| `refAHL_all_conf` | RefAHL, all confidences, no leakage filter. |
| `refAHL_conf1` / `_conf2` / `_conf3` | RefAHL, a single confidence tier. |
| `refAHL_conf1_no_leakage` | RefAHL, confidence 1, no-leakage filtered. |
| `final_n96` | The published benchmark set (N=96): `refAHL_all_conf_no_leakage` minus `Vib_001c`. |
| `paper_case_studies` | Hand-picked LuxRs for the paper case studies (`include_names` allowlist). |
| `all_data` | No filters — every row in the labels CSV. |
| `inference` | Label-free mode: `labels_path`/`ligands_path` are null; used by `copal full-ranking`. |

`set.name` controls the output subdirectory, so you can fork a run without
editing a YAML: `copal analyze set=final_n96 set.name=final_n96_variant`.

## 5. `copal full-ranking` — stylized ranking workbooks

Emits one `.xlsx` workbook per model (`AF3`, `Boltz2`, `Boltz2Affinity`,
`RF3`) plus `ENSEMBLE`, with the true-AHL rank cell highlighted yellow and
the mpKT column shaded orange (< 0.3) / green (≥ 0.3) — alongside the
programmatic `full_rankings.csv`.

```bash
uv run copal full-ranking set=final_n96     # benchmark: uses labels + ligand metadata
uv run copal full-ranking set=inference     # inference: no labels needed
```

In `set=inference` mode there is no ground truth, so the yellow highlight
is absent and the meta columns are blank — the rankings and mpKT colouring
still render.

## 6. `copal sim-matrix` — precomputed ligand-similarity matrices

`copal analyze`'s `chemical_similarity` artifact compares predicted ligands
to the true ligand by maximum-common-substructure (MCS) similarity.
Computing MCS on the fly is slow, so `copal sim-matrix` precomputes the full
pairwise matrix once:

```bash
uv run copal sim-matrix sim_matrix.method=mcs_penalized
```

`sim_matrix.method` selects the similarity definition — `mcs_order`,
`mcs_any`, or `mcs_penalized` (see [`docs/artifacts.md`](docs/artifacts.md)
for the differences). The matrix is written to
`data/chem_sim_matrices/<method>.csv`. Feed it back to analysis with
`copal analyze … analysis.chem_sim_matrix=data/chem_sim_matrices/mcs_penalized.csv`;
if that override is omitted, `copal analyze` falls back to on-the-fly
`mcs_order` computation.

## 7. Optional: filter by sequence similarity

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

## 8. Common knobs

### `copal predict` flags

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

### `copal analyze` / `full-ranking` Hydra overrides

These are `key=value` config overrides, not `--flags`.

| Override | Default | Effect |
| -------- | ------- | ------ |
| `set=<name>` | `refAHL_all_conf_no_leakage` | Which benchmark set (see Section 4). |
| `set.name=<custom>` | the set's own name | Rename the output subdirectory. |
| `analysis.aggregation=<m>` | `sum` | Ensemble aggregation. `rrf` (reciprocal rank fusion) is the main alternative; also `mean`, `median`, `min`, `max`. |
| `analysis.rrf_k=N` | 60 | RRF constant (only relevant when `aggregation=rrf`). |
| `analysis.chem_sim_matrix=<path>` | null | Precomputed ligand-similarity matrix; null = on-the-fly `mcs_order`. |
| `analysis.skip.ablation=true` | false | Skip the (slow) model-combination ablation. |
| `paths.workflow_dir=<dir>` | `outputs_official_workflow` | Where `per_model_rankings.csv` + `ensemble_rankings.csv` are read from. |
| `paths.output_root=<dir>` | `<workflow_dir>/results` | Where artifact groups are written. |
| `sim_matrix.method=<m>` | `mcs_order` | (`copal sim-matrix` only) similarity method. |
