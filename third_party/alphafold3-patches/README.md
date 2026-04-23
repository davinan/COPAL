# AlphaFold 3 7-decimal confidence patch

The default AlphaFold 3 source rounds the per-prediction confidence
metrics (ipTM, pTM, chain_pair_iptm, …) to a few decimals when writing
`*_summary_confidences.json`. For ensemble ranking — where many
predictions can sit within a small confidence delta — that rounding
discards signal.

This directory ships a one-file modification (`confidence_types.py`) that
rounds those values to **7 decimal places** instead, plus a wrapper script
that loads the modified file at runtime so you don't have to touch your
AlphaFold 3 install.

## Files

| File | License | Purpose |
| ---- | ------- | ------- |
| `confidence_types.py` | CC-BY-NC-SA 4.0 | Modified copy of AF3's `src/alphafold3/model/confidence_types.py` (7-decimal rounding in `to_json`). |
| `run_af3_with_host_confidence.py` | MIT | Runtime injector: imports the host `confidence_types.py` *before* AF3 emits JSON, then invokes `run_alphafold.main`. |
| `LICENSE` | — | CC-BY-NC-SA 4.0 license text covering `confidence_types.py`. |
| `NOTICE` | — | Provenance and modification summary. |

## Two install routes

### Route A — runtime injection (recommended)

No edit to your AF3 clone. Point env vars at your install and call the
wrapper instead of `run_alphafold` directly inside the container:

```bash
export AF3_ROOT=/opt/alphafold3                 # your AF3 source clone
export COPAL_AF3_HOST_CONFIDENCE_TYPES=$(realpath confidence_types.py)

apptainer exec /path/to/alphafold3.sif \
    python run_af3_with_host_confidence.py \
        --input_dir  $INPUT_DIR \
        --output_dir $OUTPUT_DIR \
        --db_dir     $AF3_DB_PATH \
        --run_data_pipeline true
```

Make sure both `AF3_ROOT` and the host `confidence_types.py` are inside a
path the container has bind-mounted.

### Route B — in-place edit

If you control your AF3 clone, copy the patched file over the original:

```bash
cp confidence_types.py "$AF3_ROOT/src/alphafold3/model/confidence_types.py"
```

After that you can call `run_alphafold` normally; you can ignore the
wrapper script and the env vars.

> **Reminder**: AlphaFold 3 source is CC-BY-NC-SA 4.0 (non-commercial),
> and the AlphaFold 3 model parameters are governed by separate Google
> DeepMind terms of use. Both apply regardless of which install route you
> pick. See `LICENSE`, `NOTICE`, and the repo-root `NOTICE`.
