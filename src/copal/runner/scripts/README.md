# `copal/runner/scripts/`

Shell wrappers invoked by `copal predict` once per (model, chunk). They
are **templates** — they wire together the generic "read inputs → run
model → stage outputs" flow, but they rely on a handful of env vars that
point at your local AlphaFold 3 / Boltz2 / RosettaFold3 installations.

## Default behaviour

`copal predict` emits SLURM commands that call these scripts directly by
their path inside the installed package:

```
<python_package>/copal/runner/scripts/{af3,boltz2,rf3}.sh <chunk_inputs> <chunk_outputs>
```

Each script checks for its required env vars and errors out clearly if
they're missing.

## How to customise

You have two options — pick whichever fits your cluster:

1. **Use the bundled scripts, just set env vars.** This is the simplest
   path. Before running `copal predict` (or inside the sbatch launcher),
   export the variables documented in each script header:

   ```bash
   export AF3_CONTAINER=/path/to/alphafold3.sif
   export AF3_DB_PATH=/databases/af3_DB
   export AF3_ROOT=/opt/alphafold3

   export BOLTZ2_BIN=/path/to/run_boltz.sh    # or leave unset if `boltz` is on PATH

   export RF3_INSTALL=/opt/foundry
   export RF3_CKPT=$RF3_INSTALL/rf3_latest.pt
   ```

2. **Use your own wrappers.** If your cluster needs something custom
   (module loads, Singularity flags, a different launcher), copy the
   template somewhere writable, edit it, and point COPAL at it:

   ```bash
   export COPAL_AF3_SCRIPT=/path/to/my/af3.sh
   export COPAL_BOLTZ2_SCRIPT=/path/to/my/boltz2.sh
   export COPAL_RF3_SCRIPT=/path/to/my/rf3.sh
   ```

   These overrides are read by `copal.runner.base.resolve_script_path`
   and win over the bundled copies.

## Contract

Every wrapper must accept two positional arguments — `<INPUT_DIR>` and
`<OUTPUT_DIR>` — and leave AF3/Boltz2/RF3 output files inside the output
directory in the shape the model natively produces. The analysis
pipeline expects:

| Model | Output layout (relative to `<OUTPUT_DIR>`) |
| ----- | ------------------------------------------- |
| AF3 | `<name>/<name>_model.cif` and `<name>/seed-*_sample-*/` |
| Boltz2 | `boltz_results_<chunk>/predictions/<name>/{model_*,pae_*,confidence_*,affinity_*}` |
| RF3 | `<name>.score`, `<name>_model_*.cif.gz` |

> `<name>` is the per-prediction complex name (`{receptor_id}_{ligand_id}`).
