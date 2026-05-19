# SLURM notes

`copal predict` doesn't run inference itself — it emits a `.sbatch` array
script that you submit. A few notes on what that script expects.

## Resource defaults

The defaults baked into `copal.pipeline.run_pipeline` are tuned for our
cluster; override them by editing the call site or extending the CLI.

| Argument | Default | Notes |
| -------- | ------- | ----- |
| `partition` | `gpu` | SLURM partition. |
| `gres` | `gpu` | Generic resource. AF3 and Boltz2 want one full GPU. |
| `time` | `10:00:00` | Per-task wall time. |
| `mem` | `32G` | Per-task memory. |
| `constraint` | `[A100\|H200\|L40\|A4000\|A5000\|A6000]` | Supported GPU SKUs. |
| `max_parallel` | `500` | `--array=...%N` throttle. |

## Per-prediction expected runtime (rough)

| Model | Time | GPU |
| ----- | ---- | --- |
| AlphaFold 3 (no MSA gen) | 5–10 min | ≥A4000 |
| Boltz2 | 2–5 min | ≥A4000 |
| RosettaFold3 | 2–4 min | ≥A4000 |

With `--batch-size 50`, expect 2–4 hours per array task.

## Resuming

`copal predict` rescans `outs/` on every invocation and only emits
inputs for the missing predictions. If you have a JSON record of
known-good (receptor, ligand) pairs, pass it via `--done-lookup` to skip
the disk scan. The format is:

```json
{ "RECEPTOR_LIGAND": ["af3", "boltz2"], ... }
```

(IDs uppercased.)

## Organizing outputs before analyze

```bash
uv run copal organize --pred-dir run/outs --dst-dir run/store
```

This migrates the per-model raw outputs into the hierarchical
`<RECEPTOR>/{structures,confidences}/<name>/<model>/` store. That store is
the input to the **collect** step, which builds `predictions.csv`
(see [`../REPRODUCE.md`](../REPRODUCE.md) Path B). The Hydra-driven
`copal analyze` then reads the postprocessed `per_model_rankings.csv` —
it does not take a `--basepath` flag.
