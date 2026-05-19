# Quickstart example

Two demo receptors × three demo AHL ligands (the cross-product gives 6
predictions, batched into one chunk).

```
receptors.csv  →  receptor_id, seq
ligands.csv    →  ligand_id, smiles
run.sh         →  predict → organize → collect+postprocess → full-ranking
```

## Run it

```bash
cd examples/quickstart
./run.sh
```

`run.sh` is a four-phase pipeline:

1. **`copal predict`** — schedules the 6 predictions and runs them inline
   (`--executor local`). This phase needs AF3 / Boltz2 / RF3 installed; see
   [`../../INSTALLATION.md`](../../INSTALLATION.md). On a SLURM cluster, drop
   `--executor local` from `run.sh` and `sbatch out/sarray/copal.sbatch`.
2. **`copal organize`** — gathers the raw model outputs from `out/outs/`
   into the hierarchical store `out/store/`.
3. **collect + postprocess** — builds `out/workflow/per_model_rankings.csv`
   and `ensemble_rankings.csv` from the organized outputs. (This step uses
   the legacy `python -m copal.analysis.pipeline` entrypoint — the Hydra
   `copal analyze` consumes those CSVs but does not build them.)
4. **`copal full-ranking set=inference`** — emits stylized per-model
   ranking workbooks.

Output tree:

```
out/
├── scheduled_ensemble.csv
├── inputs/{af3,boltz2,rf3}/
├── sarray/copal.{cmds,sbatch}
├── outs/                          # raw model outputs
├── store/                         # hierarchical store (copal organize)
├── workflow/                      # per_model_rankings.csv, ensemble_rankings.csv
└── results/inference/full_rankings/
    ├── AF3.xlsx  Boltz2.xlsx  Boltz2Affinity.xlsx  RF3.xlsx  ENSEMBLE.xlsx
    └── full_rankings.csv
```

## Why `full-ranking set=inference`, not `copal analyze`?

The two demo receptors have **no curated cognate-AHL ground truth**, and
`copal analyze` requires a labelled benchmark set (it computes top-K
accuracy, Kendall τ, etc. against known answers). `copal full-ranking
set=inference` is the label-free path: it ranks every ligand per receptor
by each model's confidence score plus an ensemble, with no ground truth
needed.

To run the full benchmark analysis (`copal analyze`) on the labelled paper
dataset, see [`../../docs/reproduce_paper.md`](../../docs/reproduce_paper.md).

## Your own data

Copy `receptors.csv` / `ligands.csv`, replace the contents (any number of
rows), and re-run `run.sh`. Column schemas — including the optional
`msa_path` — are in [`../../docs/input_formats.md`](../../docs/input_formats.md).
