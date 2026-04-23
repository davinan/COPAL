# Task: smoke-test the COPAL pipeline end-to-end on a local GPU

You are running on a machine that has:

- A CUDA GPU suitable for AlphaFold 3 / Boltz2 / RosettaFold3 inference
  (A100-80GB confirmed available; AF3 wants ≥24 GB VRAM for full pipelines).
- `uv`, `git`, `bash`, `apptainer` (or `singularity`), and `hmmer` binaries
  at `/software/hmmer/latest/bin/` available.
- Institute-provided AF3 apptainer image at `/software/containers/af3.sif`
  with AF3 databases at `/databases/lab/af3_DB/` and model params (the
  `af3.bin`/`af3.bin.zst` pair) at `/databases/alphafold/`. No Boltz2 or
  Foundry installs — you will install those fresh per `INSTALLATION.md`.

Your job is to prove that COPAL runs the **full pipeline** locally
(no SLURM) on **1 receptor × 2 ligands** from the paper benchmark,
using a pre-computed MSA, and produces analysis outputs. **Report what
worked, what broke, and what you had to change to make it work.**

The pipeline has several *known* pre-existing issues (see the
**Expected landmines** section below). Your job is to discover them,
fix them, and proceed — not to work around silently. Treat every bug as
a patch opportunity; leave the repo in a better state than you found it.

## Repository

```
/home/dnan/projects/luxr/src/publication_repo
```

Start by reading:
- `README.md` — top-level overview
- `INSTALLATION.md` — env-var table, per-model uv-venv layout
- `USAGE.md` — CLI reference + local-executor walkthrough
- `src/copal/runner/scripts/README.md` — shell-wrapper contract
- `docs/pipeline_overview.md` — data flow (`outs/` → `organize` → `store/` → `analyze`)

## Install

Follow `INSTALLATION.md` end-to-end. The "Quick install of all external
models" subsection in §1 is a single copy-paste block that clones AF3
(+ applies the 7-decimal ipTM patch Route B), creates a dedicated uv
venv for Boltz2, and a dedicated uv venv for Foundry / RF3 with the
`[rf3]` extras, then downloads the RF3 weights via `foundry install`.
After it succeeds, run the per-venv smoke tests at the bottom of
`INSTALLATION.md` (§10) before proceeding.

## Exact commands to execute (edit paths as needed)

```bash
cd /path/to/copal
REPO=$PWD

# 1. Environment — point the bundled wrappers at the installs.
#    On the IPD cluster these should work as-is; elsewhere override the
#    container / DB / model-dir paths. AF3_MODEL_DIR must contain the
#    af3.bin (or af3.bin.zst) model-parameter file.
#    NOTE: the bundled `af3.sh` is known-broken (no --model_dir, and
#    --run_data_pipeline true by default); you'll need a custom wrapper
#    pointed at via COPAL_AF3_SCRIPT. See the landmines table below.
export AF3_CONTAINER=/software/containers/af3.sif
export AF3_DB_PATH=/databases/lab/af3_DB
export AF3_ROOT=$REPO/external/alphafold3
export AF3_MODEL_DIR=/databases/alphafold       # dir containing af3.bin
export BOLTZ2_BIN=$REPO/external/boltz/.venv/bin/boltz
export RF3_INSTALL=$REPO/external/foundry
export RF3_CKPT=$RF3_INSTALL/rf3_latest.pt

# 2a. Extract the bundled MSA archive (31 MB → 263 MB of .a3m files).
mkdir -p data/benchmark_msas
tar -xJf data/benchmark_msas.tar.xz -C data/benchmark_msas/

# 2b. Stage 1 LuxR + 2 AHL ligands FROM THE BENCHMARK, with the MSA path.
#     NOTE: msa_path in the benchmark CSV is repo-root-relative. Convert to
#     absolute so the model wrappers can resolve it regardless of CWD.
mkdir -p copal_run_test
uv run python - <<'PY'
import pandas as pd
from pathlib import Path
repo = Path.cwd().resolve()
out = Path('copal_run_test')
labels = pd.read_csv('data/benchmark_data_all_corrected.csv')
recv = (labels[labels['refahl_id'] == 'Pse_001_3var']
        [['refahl_id','seq','msa_path']]
        .rename(columns={'refahl_id':'receptor_id'}))
recv['msa_path'] = recv['msa_path'].apply(lambda p: str(repo / p))
recv.to_csv(out/'receptors.csv', index=False)

ligs = pd.read_csv('data/Final_confirmed_ligands_58_annotated.csv').head(2)
ligs = ligs[['AHLid','smi']].rename(columns={'AHLid':'ligand_id','smi':'smiles'})
ligs.to_csv(out/'ligands.csv', index=False)

print('receptor:', recv.iloc[0]['receptor_id'])
print('msa_path:', recv.iloc[0]['msa_path'])
print('msa exists?', Path(recv.iloc[0]['msa_path']).exists())
print('ligands:', list(ligs['ligand_id']))
PY

# 3. Run inference LOCALLY (not SLURM), serial.
#    NOTE: the bundled af3.sh (a) omits --model_dir so AF3 falls back to
#    $HOME/models and crashes, and (b) hardcodes --run_data_pipeline true
#    which runs HMMer searches even though the input JSON already has
#    unpairedMsaPath. You need a custom wrapper that passes --model_dir
#    and `--run_data_pipeline false`. Write it under external/wrappers/
#    and set COPAL_AF3_SCRIPT before running predict.
#
#    NOTE: Boltz 2.2.1 ships torch 2.11.0+cu130 but pulls
#    cuequivariance-ops-torch-cu12 (CUDA-12 wheels) — at runtime the
#    triangular-mult kernel fails with `undefined symbol: cublasGemmGroupedBatchedEx`.
#    Either write a custom boltz2.sh that passes `--no_kernels` (upstream
#    Boltz flag, ~1.5-2x slower but numerically identical) or patch the
#    venv to pin torch to a cu12 build. Point COPAL_BOLTZ2_SCRIPT at your
#    custom wrapper.

uv run copal predict \
    --receptors copal_run_test/receptors.csv \
    --ligands   copal_run_test/ligands.csv \
    --out-dir   copal_run_test/pipeline \
    --models    af3,boltz2,rf3 \
    --batch-size 2 \
    --executor  local \
    --workers   1

# 4. Organize raw outputs into the hierarchical store copal analyze expects.
#    NOTE: `copal organize` has known bugs that need fixing before this
#    will succeed on receptor names like `Pse_001_3var`:
#    - organizer.py `get_range_string` does `int(name[5:])` which explodes
#      on any non-`p###` receptor ID. Fix: fall back to a flat chunk-range
#      (e.g. "all") when suffix isn't numeric, or drop the chunking entirely.
#    - organizer.py `name.split("_")[:2]` splits luxr/ligand assuming a
#      single underscore. Fix: use `name.rsplit("_AHL", 1)` to separate
#      receptor (may contain underscores) from ligand.
#    - AF3PathManager calls `.lower()` on the name because older AF3
#      builds lowercased output dirs; modern AF3 preserves case so the
#      path manager finds zero files and silently skips AF3. Fix: drop the
#      three `.lower()` calls.
#    - RF3PathManager is hard-coded for an older RF3 output layout
#      (`*_model_4.cif.gz` + top-level `.score`/`_metrics.csv`). Current
#      production-branch RF3 writes `<name>/<name>_model.cif` with
#      per-sample `seed-*_sample-*/` subdirs and a
#      `<name>_ranking_scores.csv`. Rewrite `RF3PathManager` to walk the
#      new shape and synthesize the legacy `.score`/`_metrics.csv` from
#      per-sample `*_summary_confidences.json` so downstream analyze
#      keeps reading the same schema.

uv run copal organize \
    --pred-dir copal_run_test/pipeline/outs \
    --dst-dir  copal_run_test/pipeline/store \
    --models   af3,boltz2,rf3

# 5. Collect + analyze.
#    NOTE: analyze's `compute_all_rankings` hard-requires exactly 58
#    unique `luxr_ahl` pairs per receptor — for a 1×2 smoke test it
#    raises `RuntimeError: No rankings were produced`. Two separate fixes
#    are in scope:
#    (a) Plumb an `expected_ligand_count` through the analyze CLI that
#        defaults to the ligands-CSV row count, so the gate becomes
#        dynamic ("all scheduled predictions are present") rather than
#        hardcoded 58.
#    (b) The Typer `analyze` wrapper in cli.py doesn't expose the
#        underlying argparse `--skip-rankings`/`--skip-ensemble`/
#        `--skip-ablation`/`--skip-done-json` flags, and the downstream
#        `run_generate` blindly reads the rankings CSVs without checking
#        they exist (so even if you skip the postprocess step, generate
#        crashes). Add a `--skip-generate` flag and a pre-check for the
#        rankings CSVs.

uv run copal analyze \
    --workflow-dir copal_run_test/pipeline/workflow \
    --output-root  copal_run_test/pipeline/results \
    --labels       data/benchmark_data_all_corrected.csv \
    --ligands      data/Final_confirmed_ligands_58_annotated.csv \
    --basepath     copal_run_test/pipeline/store/all \
    --af3-iptm-source json \
    --expected-ligand-count 2
```

## What to verify

After step 3 (`copal predict`):

- `copal_run_test/pipeline/scheduled_ensemble.csv` has 2 rows (one per
  ligand) with `is_done_*` all False initially.
- `copal_run_test/pipeline/inputs/{AF3,Boltz2,RF3}/0/` contains
  `Pse_001_3var_AHL01` and `Pse_001_3var_AHL02` input files
  (`.json` for AF3/RF3, `.yaml` for Boltz2). Open each and confirm the
  MSA path in `unpairedMsaPath` (AF3), `msa` (Boltz2), and `msa_path`
  (RF3 component) points at the correct `.a3m` file — **not null**.
- `copal_run_test/pipeline/sarray/copal.cmds` has 3 lines (one per model).
- `copal_run_test/pipeline/logs/task-{001,002,003}.{out,err}` show the
  model runs. `[task-NNN] OK` printed to the parent terminal for each.
- `copal_run_test/pipeline/outs/{AF3,Boltz2,RF3}/0/` contains the raw
  prediction files (CIFs, confidence JSONs, Boltz2 affinity JSON, RF3
  `ranking_scores.csv`). See `docs/pipeline_overview.md` for the exact
  layout.
- Open one AF3 `*_summary_confidences.json` and confirm ipTM values are
  written with 7 decimal places (e.g. `0.9702999`, not `0.97`). This
  confirms the Route-B patch is live inside the container via AF3_ROOT.

After step 4 (`copal organize`):

- `copal_run_test/pipeline/store/all/PSE_001_3VAR/{structures,confidences}/Pse_001_3var_AHL0{1,2}/{af3,boltz2,rf3}/`
  should exist and be non-empty.

After step 5 (`copal analyze`):

- `copal_run_test/pipeline/workflow/predictions.csv` exists with one row per
  (receptor, ligand, model, batch) prediction with an `iptm` column.
- `copal_run_test/pipeline/workflow/per_model_rankings.csv` (8 rows for 1×2 × 4 model keys).
- `copal_run_test/pipeline/workflow/ensemble_rankings.csv` (2 rows).
- `copal_run_test/pipeline/results/` contains `ranking_performance/`,
  `full_rankings/`, `acyl_moiety/`, `chemical_similarity/`, `mpkt/`,
  `ablation/` subdirs. `full_rankings/<MODEL>.xlsx` should be present
  for AF3, Boltz2, Boltz2Affinity, RF3, ENSEMBLE.
- `ablation_table.csv` may be empty — expected when the scheduled ligand
  set doesn't contain the receptor's cognate ligand (AHL04 for
  `Pse_001_3var`; our 2-ligand subset is AHL01 + AHL02).

## Expected landmines (summary)

The pipeline has several pre-existing bugs. Fix each in source and
commit the fix alongside your patches — do not work around silently.

| Where | Symptom | Root cause / fix sketch |
|-------|---------|------------------------|
| `runner/scripts/af3.sh` | AF3 crashes with `FileNotFoundError: '/home/dnan/models'` | Wrapper omits `--model_dir`. Also `--run_data_pipeline true` is wrong when you have precomputed MSAs. Write a custom wrapper; `COPAL_AF3_SCRIPT` points to it. |
| Boltz2 install | `ImportError: libcue_ops.so: undefined symbol: cublasGemmGroupedBatchedEx, version libcublas.so.12` | Boltz 2.2.1 wheel mismatch (torch cu130 vs cuequivariance-ops cu12). Custom boltz2.sh passing `--no_kernels` is the pragmatic fix. |
| `INSTALLATION.md` §6 | `rf3 fold` crashes with `ModuleNotFoundError: cuequivariance_torch` | Doc missed `[rf3]` extras. Update the doc + verify. |
| Foundry `rf3 fold` + modelhub-style ckpt | `TypeError: FabricTrainer.__init__() got an unexpected keyword argument 'symmetry_resolution_in_training'` | The ckpt was trained on a branch that subclasses FabricTrainer with an extra kwarg the `production` branch's base class rejects. Download the matching production-branch ckpt via `foundry install rf3` (NOT a modelhub ckpt). |
| `organizer.py:35` | `ValueError: invalid literal for int() with base 10: ''` on `Pse_001_3var` | `get_range_string` assumes `p###`-style IDs. Flatten chunk layout or guard the cast. |
| `organizer.py` name split | Luxr/ligand split is wrong for multi-underscore receptor IDs | `name.split("_")[:2]` is brittle. Use `rsplit("_AHL", 1)` or a single rsplit. |
| `prediction_path_managers.AF3PathManager` | AF3 predictions silently skipped by organize | `.lower()` on name; modern AF3 preserves case. Drop the `.lower()` calls. |
| `prediction_path_managers.RF3PathManager` | RF3 predictions silently skipped | Expects `*_model_4.cif.gz` + `.score`/`_metrics.csv`. Rewrite to walk the current production-branch layout and synthesize the legacy compat CSVs. |
| `analysis/pipeline.compute_all_rankings` | `RuntimeError: No rankings were produced` on 1×2 | Hardcoded 58 ligand gate. Plumb `--expected-ligand-count` through from the ligands CSV row count. |
| `cli.py` `analyze` | `--skip-rankings`/`--skip-ensemble`/`--skip-ablation`/`--skip-done-json` missing from Typer wrapper (they exist on argparse) | Expose them; also add `--skip-generate` and auto-skip generate when rankings CSV missing. |
| `analysis/pipeline.build_ablation_table` | `KeyError: 'top 6 success rate'` when results are empty | Guard the empty-df path; return empty DataFrame with the expected schema. |

## Key source files to grep when debugging

```
src/copal/pipeline.py                    — schedule + dispatch
src/copal/runner/base.py                 — make_inputs, script path resolution
src/copal/runner/{af3,boltz2,rf3}.py     — per-model input construction + commands
src/copal/runner/local.py                — run_commands_locally (local executor)
src/copal/runner/scripts/*.sh            — bundled wrappers (needs overrides)
src/copal/organizer.py                   — raw outs → organized store
src/copal/prediction_path_managers.py    — per-model path managers
src/copal/analysis/pipeline.py           — analyze (collect/postprocess/generate)
src/copal/io/schema.py                   — receptor/ligand CSV validation + aliases
third_party/alphafold3-patches/          — AF3 7-decimal confidence mod + injection wrapper
```

## Ground truth for the test receptor

`Pse_001_3var` has `true_ahl_abbrev = 3OC10` (AHL04) in the benchmark.
Your 2-ligand subset (AHL01, AHL02) does NOT include the cognate —
that's fine for a smoke test; ranking metrics won't be meaningful but
the pipeline should still write rankings/ensemble CSVs. Ablation
metrics will be empty (no cognate in the scheduled set) — this is
expected, not a failure.

## What to report back

A ≤400-word summary covering:

1. **Pass/fail verdict** for install (INSTALLATION.md quick-install block) plus each of steps 1–5 in the command block.
2. **Every source-code / wrapper / doc patch** you applied, as a
   concise diff list (file path + one-line summary per change). Be
   explicit about which landmines from the table above you hit.
3. **Wall-clock time** for each of the 3 per-model tasks (from
   `[task-NNN] OK` log lines or `outs/*/0/*` mtimes).
4. **Any surprising output** — empty CSVs, unexpected warnings, ipTM
   values that look wrong (e.g. all zeros, or >1), ligand SMILES
   handling issues, etc.
5. If the analyze step fails, **include the full traceback** and the
   first 30 lines of the relevant task-NNN.err.

Do NOT skip failing steps silently — debug them, fix them, retry. If a
wrapper script doesn't fit your install, write your own, point the env
var at it, and include the wrapper in the report.

You are authorized to modify files anywhere in `src/copal/` and any
documentation (`*.md`) to fix pre-existing bugs. Keep edits minimal and
surgical — don't refactor beyond what a fix requires. Do not touch
`data/`, `third_party/alphafold3-patches/`, or your own `external/`
clones (those are model-repo upstream code).
