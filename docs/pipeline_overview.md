# Pipeline overview

COPAL has two stages — **inference** (`copal predict`) and **analysis**
(`copal analyze`) — connected through a working directory whose layout is
the only contract between them.

```
                   ┌───────────────────────────────────────────────┐
                   │                copal predict                  │
                   │                                               │
   receptors.csv ──┤  • cross-product → scheduled_ensemble.csv     │
   ligands.csv  ──┤  • per-batch resume scan (outs/)              │
                   │  • input_manager.make_inputs()                │
                   │      → inputs/{af3,boltz2,rf3}/*.{json,yaml} │
                   │  • emit copal.cmds                            │
                   │                                               │
                   │  --executor slurm  →  make_sbatch_file()      │
                   │                        → copal.sbatch +       │
                   │                          submit_sarray.sh     │
                   │  --executor local  →  run_commands_locally()  │
                   │                        → logs/task-NNN.out    │
                   └───────────────┬───────────────────────────────┘
                                   │ sbatch OR direct shell
                                   ▼
                   ┌───────────────────────────────────────────────┐
                   │  Per-chunk shell wrappers                     │
                   │   ├── af3.sh    → outs/AF3/<chunk>/...        │
                   │   ├── boltz2.sh → outs/Boltz2/<chunk>/...     │
                   │   └── rf3.sh    → outs/RF3/<chunk>/...        │
                   └───────────────┬───────────────────────────────┘
                                   │
                                   ▼
                   ┌───────────────────────────────────────────────┐
                   │                copal analyze                  │
                   │                                               │
                   │  Phase 1: collect (optional)                  │
                   │    outs/  → predictions.csv                   │
                   │    (skippable via --input-csv)                │
                   │                                               │
                   │  Phase 2: postprocess                         │
                   │    processed CSV + labels →                   │
                   │      per_model_rankings.csv                   │
                   │      ensemble_rankings.csv                    │
                   │      ablation_table.csv                       │
                   │      done_pairs.json                          │
                   │                                               │
                   │  Phase 3: generate                            │
                   │    workflow outputs + labels + ligands →      │
                   │      ranking_performance/                     │
                   │      full_rankings/                           │
                   │      ablation/                                │
                   │      chemical_similarity/  (if rdkit)         │
                   │      acyl_moiety/                             │
                   │      mpkt/                                    │
                   └───────────────────────────────────────────────┘
```

## Working directory layout

A single `--out-dir` (predict) plus `--workflow-dir` and `--output-root`
(analyze) carry the pipeline. A typical layout:

```
run/
├── scheduled_ensemble.csv        # one row per (receptor, ligand)
├── inputs/
│   ├── af3/    *.json
│   ├── boltz2/ *.yaml
│   └── rf3/    *.json
├── sarray/
│   ├── copal.cmds                # one shell line per batch
│   └── copal.sbatch              # SLURM array script
├── outs/                         # raw model outputs
│   ├── AF3/<chunk>/<name>/...
│   ├── Boltz2/<chunk>/.../boltz_results_<chunk>/predictions/<name>/...
│   └── RF3/<chunk>/<name>/{<name>_model.cif, seed-*_sample-*/, ranking_scores.csv}
├── workflow/                     # postprocess outputs
│   ├── predictions.csv           (collect)
│   ├── per_model_rankings.csv
│   ├── ensemble_rankings.csv
│   ├── ablation_table.csv
│   └── done_pairs.json
└── results/                      # generate outputs
    ├── ranking_performance/...
    ├── full_rankings/...
    ├── ablation/...
    ├── chemical_similarity/...
    ├── acyl_moiety/...
    └── mpkt/...
```

## Modes

`copal analyze` automatically switches between two modes based on what's
in `--labels`:

- **Benchmark mode** — labels CSV has `true_ahl_abbrev`. Phase 3 emits
  full ranking-performance metrics (top-K accuracy, Kendall τ, Spearman ρ),
  ablation tables, and chemical-similarity diagnostics.
- **Prediction mode** — labels CSV omits the ground-truth column. Phase 3
  still produces the per-receptor ranked ligand lists; performance
  metrics are simply skipped.
