# Pipeline overview

COPAL is a six-subcommand CLI. A typical run flows:

```
copal predict  →  copal organize  →  copal analyze     (benchmark)
                                  └→  copal full-ranking (label-free)

copal sim-matrix   — one-off: precompute a ligand-similarity matrix
copal seqid-filter — one-off: sequence-identity leakage filter
```

- **`copal predict`** — schedule every receptor × ligand pair and emit
  per-model inputs + a SLURM array (or run locally).
- **`copal organize`** — gather raw model outputs into a hierarchical store.
- **`copal analyze`** — Hydra-driven benchmark analysis (needs ground-truth
  labels). Emits seven artifact groups.
- **`copal full-ranking`** — Hydra-driven stylized ranking workbooks. Works
  with or without labels (`set=inference` for the label-free path).
- **`copal sim-matrix`** — precompute a pairwise ligand chemical-similarity
  matrix that `copal analyze` can look up instead of recomputing.
- **`copal seqid-filter`** — BLAST + global-alignment sequence-identity
  leakage filter.

```
                   ┌───────────────────────────────────────────────┐
                   │                copal predict                  │
   receptors.csv ──┤  • cross-product → scheduled_ensemble.csv      │
   ligands.csv  ───┤  • per-batch resume scan (outs/)               │
                   │  • make_inputs() → inputs/{af3,boltz2,rf3}/    │
                   │  • emit copal.cmds                             │
                   │  --executor slurm → copal.sbatch               │
                   │  --executor local → run inline, logs/task-*    │
                   └───────────────┬───────────────────────────────┘
                                   │ sbatch OR direct shell
                                   ▼
                   ┌───────────────────────────────────────────────┐
                   │  Per-chunk shell wrappers                      │
                   │   af3.sh / boltz2.sh / rf3.sh → outs/<MODEL>/  │
                   └───────────────┬───────────────────────────────┘
                                   ▼
                          copal organize
                       (raw outs/ → hierarchical store/)
                                   │
                   ┌───────────────┴─── collect + postprocess ──────┐
                   │  per_model_rankings.csv + ensemble_rankings.csv │
                   │  in paths.workflow_dir                          │
                   └───────────────┬───────────────────────────────┘
                                   │
                       ┌───────────┴────────────┐
                       ▼                        ▼
                copal full-ranking         copal analyze
                (stylized .xlsx;           (7 artifact groups under
                 set=inference OK)          results/<set.name>/)
```

## The analysis source of truth

Every analysis artifact derives from two CSVs in `paths.workflow_dir`:

- **`per_model_rankings.csv`** — long form: one row per (LuxR, AHL, model)
  with `rank`, `iptm`, `affinity_probability_binary`.
- **`ensemble_rankings.csv`** — the aggregated ensemble ranking.

These are produced by the **postprocess** step. `copal analyze` runs
postprocess only if `analysis.skip.postprocess=false` (default: `true`,
i.e. it assumes the CSVs already exist). The **collect** step (raw outputs →
`predictions.csv`, the input to postprocess) is *not* part of the Hydra
`copal analyze`; it currently lives only in the deprecated
`python -m copal.analysis.pipeline` entrypoint — see
[`REPRODUCE.md`](../REPRODUCE.md) Path B.

## Working directory layout

`copal predict` writes to its `--out-dir`. The analysis commands read/write
under `paths.workflow_dir` (default `outputs_official_workflow/`) and
`paths.output_root` (default `<workflow_dir>/results`). A typical tree:

```
run/                              # copal predict --out-dir
├── scheduled_ensemble.csv        # one row per (receptor, ligand)
├── inputs/{af3,boltz2,rf3}/      # per-model JSON / YAML
├── sarray/copal.{cmds,sbatch}    # SLURM array script
└── outs/                         # raw model outputs
    ├── AF3/<chunk>/<name>/...
    ├── Boltz2/<chunk>/...
    └── RF3/<chunk>/...

outputs_official_workflow/         # paths.workflow_dir
├── predictions.csv               # collected raw table (collect step)
├── per_model_rankings.csv        # postprocess output  ─┐ analysis
├── ensemble_rankings.csv         # postprocess output  ─┘ source of truth
└── results/                      # paths.output_root
    └── <set.name>/               # e.g. final_n96/
        ├── ranking_performance/
        ├── mpkt/
        ├── acyl_moiety/
        ├── chemical_similarity/
        ├── score_vs_rank/
        ├── full_rankings/        # full_rankings.csv (+ .xlsx via full-ranking)
        └── ablation/

data/chem_sim_matrices/           # copal sim-matrix output
└── <method>.csv                  # 58×58 pairwise ligand-similarity matrix
```

The seven artifact groups under `results/<set.name>/` are documented in
[`artifacts.md`](artifacts.md).

## Configuration (Hydra)

`copal analyze`, `copal full-ranking`, and `copal sim-matrix` are all
Hydra-driven: arguments are `key=value` config overrides, not `--flags`.
The config tree lives at `src/copal/config/`:

- `set/*.yaml` — benchmark-set definitions (which LuxRs, which labels CSV).
- `analysis/default.yaml` — aggregation method, skip flags, `chem_sim_matrix`.
- `paths/repo.yaml` — `workflow_dir`, `output_root`, `data_root` (rooted at
  `COPAL_REPO_ROOT`, defaulting to the current directory).
- `sim_matrix/default.yaml` — similarity method + output path.

See [`../USAGE.md`](../USAGE.md) §3–6 for the override syntax and the set
catalog.

## Modes

- **Benchmark mode** — `copal analyze set=<labelled-set>`. The set's labels
  CSV carries cognate-AHL ground truth; analyze emits all seven artifact
  groups (true-AHL rank, mpKT, ablation, chemical similarity, …). `copal
  analyze` **errors out** if the resolved set has no labelled rows.
- **Label-free mode** — `copal full-ranking set=inference`. No ground truth
  required; ranks every ligand per receptor by each model's confidence plus
  an ensemble, and emits stylized `.xlsx` workbooks. This is the path for
  arbitrary user data with no known cognate ligand.
