# Installation

This repository assumes Linux + SLURM + a CUDA-capable GPU partition for
inference. The analysis pipeline runs CPU-only.

Tested on Python 3.11, CUDA 13.1 (NVIDIA driver 590.48.01), RTX A4000
(16 GB VRAM). The paper benchmark was run on our lab's SLURM cluster;
a single A4000 is sufficient for all three models at the sequence
lengths for this benchmark.

## 1. Python environment

```bash
git clone https://github.com/davinan/COPAL.git
cd COPAL
uv sync                                    # creates .venv/, installs copal + deps
uv run copal --help                        # smoke test
```

### Per-model venvs

Each external model ships in its own `uv` venv under
`external/<model>/.venv/` so that AF3 / Boltz2 / RF3 dependency graphs
(PyTorch versions, CUDA pins, JAX, `cuequivariance*`, etc.) don't
collide with each other or with the top-level `copal` env. Layout:

```
external/
├── alphafold3/       # source clone only; runtime is the .sif container
├── boltz/.venv/      # Boltz2 uv venv  (python -m boltz)
└── foundry/.venv/    # RF3 uv venv     (rf3 fold ...)
```

The sections below install each model into its own venv; point the
`AF3_ROOT` / `BOLTZ2_BIN` / `RF3_INSTALL` env vars at the matching path.

### Quick install of all external models

Run this from the repo root after `uv sync`. It performs every step
§§3, 5, 6 prescribe in one pass. See those sections for the full detail
(licenses, patch routes, troubleshooting); the commentary here is kept
minimal.

```bash
REPO=$PWD
mkdir -p external

# AF3: clone source, then apply the 7-decimal ipTM patch in place
# (Route B from third_party/alphafold3-patches/README.md).
# Runtime uses a container image supplied via $AF3_CONTAINER — see §3.
(cd external && git clone https://github.com/google-deepmind/alphafold3.git)
cp third_party/alphafold3-patches/confidence_types.py \
   external/alphafold3/src/alphafold3/model/confidence_types.py

# Boltz2: dedicated uv venv
mkdir -p external/boltz && (cd external/boltz && \
    uv venv .venv --python 3.12 && \
    source .venv/bin/activate && \
    uv pip install boltz)

# RF3 / Foundry: dedicated uv venv + install, then download weights.
# Two subtleties:
#   - install with the `[rf3]` extras so `cuequivariance*` is pulled in
#     (plain `uv pip install -e .` installs cleanly but `rf3 fold`
#     crashes at runtime with ModuleNotFoundError: cuequivariance_torch).
#   - the `foundry` CLI entry point sometimes doesn't register after
#     the extras-install; a second plain `uv pip install -e .` makes
#     it show up so we can run `foundry install rf3` to fetch weights.
(cd external && git clone -b production https://github.com/RosettaCommons/foundry.git)
(cd external/foundry && \
    uv venv .venv --python 3.12 && \
    source .venv/bin/activate && \
    uv pip install -e ".[rf3]" && \
    uv pip install -e . && \
    mkdir -p checkpoints && \
    foundry install rf3 --checkpoint-dir "$PWD/checkpoints" && \
    ln -sf checkpoints/rf3_foundry_01_24_latest_remapped.ckpt rf3_latest.pt)
```

After this completes, set the env vars listed in §2 and run the
smoke tests in §10 to verify each venv works before trying a prediction.

## 2. Environment variables

COPAL vendors the AF3/Boltz2/RF3 orchestration layer — there is no
separate `model_runner` dependency. The only external pieces you bring
are the model installations themselves and, optionally, your own shell
wrappers (see step 6).

| Variable | Used by | Purpose |
| -------- | ------- | ------- |
| `AF3_CONTAINER` | `runner/scripts/af3.sh` | Apptainer/Singularity `.sif` image for AlphaFold 3. |
| `AF3_DB_PATH` | `runner/scripts/af3.sh` | AF3 sequence database directory (e.g. `/databases/af3_DB`). |
| `AF3_ROOT` | `runner/scripts/af3.sh`, AF3 wrapper | Your AlphaFold 3 source clone (so `import run_alphafold` works). |
| `AF3_MODEL_DIR` | AF3 wrapper | Directory containing the AF3 model parameters (`af3.bin` / `af3.bin.zst`). Passed through to AF3 as `--model_dir`; required by any wrapper that doesn't want AF3 to fall back to its `$HOME/models` default. |
| `AF3_WRAPPER` | `runner/scripts/af3.sh` | Override the injection wrapper path. Defaults to the bundled copy in `third_party/alphafold3-patches/`. |
| `COPAL_AF3_HOST_CONFIDENCE_TYPES` | AF3 wrapper | Override path to the patched `confidence_types.py`. Defaults to the bundled copy. |
| `BOLTZ2_BIN` | `runner/scripts/boltz2.sh` | Optional — absolute path to a `boltz` / `run_boltz.sh` launcher. If unset, `boltz` must be on PATH. |
| `BOLTZ2_DIFFUSION_SAMPLES` | `runner/scripts/boltz2.sh` | Optional, default 5. |
| `RF3_INSTALL` | `runner/scripts/rf3.sh` | Path to your Foundry/RF3 install (must contain `.venv/bin/activate`). |
| `RF3_CKPT` | `runner/scripts/rf3.sh` | Path to RF3 weights. |
| `RF3_EXTRA_ARGS` | `runner/scripts/rf3.sh` | Optional extra `rf3 fold` flags. |
| `COPAL_AF3_SCRIPT` | `copal.runner.af3` | Override the AF3 wrapper script with your own copy. |
| `COPAL_BOLTZ2_SCRIPT` | `copal.runner.boltz2` | Same, for Boltz2. |
| `COPAL_RF3_SCRIPT` | `copal.runner.rf3` | Same, for RF3. |
| `COPAL_POCKET_REF_CIF` | `copal.utils.pocket` | Reference CIF for the optional pocket-residue annotation utility. The built-in default path is `data/reference/4y17.cif`, which is **not shipped in the repo** — set this variable to your own CIF if you use pocket annotation. |
| `COPAL_REPO_ROOT` | Hydra `paths` config (`config/paths/repo.yaml`) | Root the `copal analyze` / `full-ranking` / `sim-matrix` config paths resolve against. Defaults to the current working directory; **set it explicitly** when running these commands from anywhere other than the repo root (e.g. on a cluster), so `data/` and `outputs_official_workflow/` resolve correctly. |

A `.env` style template:

```bash
# Copy to .env and edit. Source it from your shell or sbatch launcher.
export AF3_CONTAINER=/software/alphafold3.sif
export AF3_DB_PATH=/databases/af3_DB
export AF3_ROOT=$HOME/src/alphafold3
export AF3_MODEL_DIR=/databases/alphafold  # dir containing af3.bin
export BOLTZ2_BIN=$(which boltz)           # or a run_boltz.sh path
export RF3_INSTALL=$HOME/src/foundry
export RF3_CKPT=$RF3_INSTALL/rf3_latest.pt
export COPAL_REPO_ROOT=$PWD                # so `copal analyze` resolves data/ + outputs/
```

## 3. AlphaFold 3

> **Non-container install is not officially supported upstream.** The AF3
> authors distribute AF3 as a Docker image; a native `uv pip install` of
> the source tree is *not* an officially supported path. In practice,
> two working options are:
>   1. Build an Apptainer/Singularity `.sif` from the upstream AF3
>      Dockerfile (the generic route), or
>   2. Reuse an institution-provided container (e.g. the IPD cluster
>      ships one at `/software/containers/af3.sif`).
>
> Either way, the cloned `AF3_ROOT` source tree is still used at runtime:
> the bundled injection wrapper (`third_party/alphafold3-patches/`)
> supplies `run_alphafold.py` with the host-patched `confidence_types.py`.
> If you prefer, you can instead apply the patch in-place inside
> `AF3_ROOT` (Route B below) and skip the runtime injection.

1. Clone the official source (no venv needed — this is just the source
   tree the container / wrapper reads from):

   ```bash
   git clone https://github.com/google-deepmind/alphafold3.git $AF3_ROOT
   ```

2. Obtain the model parameters from Google DeepMind (per their
   [Model Parameters Terms of Use](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md))
   and place them where AF3 expects.

3. Build / pull an AF3 container (Apptainer / Singularity / Docker — see
   the upstream README).

4. **DB size + MSA shortcut.** AF3's full sequence+template database is
   ~252 GB compressed / **~630 GB unzipped** (see upstream
   [installation docs](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md)).
   If you have pre-computed MSAs (the paper benchmark ships 281 of them
   in `data/benchmark_msas.tar.xz`) and your AF3 inputs set
   `"templates": []` — which `copal.runner.af3` does — you can skip the
   DB entirely by running AF3 with `--run_data_pipeline false`. The
   bundled `src/copal/runner/scripts/af3.sh` uses `true` by default;
   copy it and flip the flag (or edit in place) if you want the
   shortcut and haven't downloaded the DB. If you stick with
   `--run_data_pipeline true`, make sure `hmmsearch`/`jackhmmer` are on
   `PATH` inside whatever environment the wrapper runs under
   (e.g. `export PATH=/software/hmmer/latest/bin:$PATH`); AF3 shells
   out to them during the data pipeline.

5. **Apptainer bind mounts.** `apptainer exec` only auto-binds a small
   set of host paths (typically `$HOME`, `/tmp`, `$PWD`). If your
   `AF3_ROOT`, `AF3_MODEL_DIR`, `AF3_DB_PATH`, or the per-run
   input/output directories live outside those (e.g. `/databases/…`,
   `/net/scratch/…`), the container will see them as missing and AF3
   will crash with a `FileNotFoundError`. Your wrapper must pass
   `--bind <host>:<host>` for each such path. The template in
   `external/wrappers/af3.sh` shows the five canonical binds
   (`AF3_ROOT`, `AF3_MODEL_DIR`, `AF3_DB_PATH`, input dir, tmp output
   dir); the bundled `src/copal/runner/scripts/af3.sh` omits them
   intentionally to stay generic, so copy+edit is the expected path.

6. **(Optional but recommended for COPAL ensemble ranking):**
   apply the 7-decimal-precision patch shipped here. See
   [`third_party/alphafold3-patches/README.md`](third_party/alphafold3-patches/README.md)
   for the runtime-injection (no AF3 source edit) or in-place-edit route.

> **Container vs source install**: the bundled `af3.sh` wrapper assumes
> an Apptainer/Singularity image. If you `uv pip install` AF3 from
> source, copy the wrapper, remove the `apptainer exec --nv "$AF3_CONTAINER"`
> prefix, and point `COPAL_AF3_SCRIPT` at your copy.

> **License**: AlphaFold 3 source is CC-BY-NC-SA 4.0 (non-commercial);
> the model weights have separate, restrictive terms of use. You are
> responsible for compliance.

### 3a. Native (no container) install

When neither a Docker daemon nor `apptainer build --fakeroot` is
available (e.g. you're on a shared cluster without subuid mappings or
the `docker` group), AF3 can be installed directly into a uv venv.
**Upstream AF3 does not officially support this path** — use it only if
the container routes are truly blocked — but in practice it works on
any Linux host with recent NVIDIA drivers and a working C++ toolchain.

**System prerequisites** (Ubuntu 22.04/24.04; adjust for other distros):

```bash
sudo apt install build-essential cmake ninja-build zstd
```

JAX ships its own bundled CUDA 12 libraries, so **no system CUDA
toolkit is required** — just an NVIDIA driver recent enough to expose
CUDA 12 user-space (driver ≥ 525 on Linux). No HMMER is needed as long
as you run with `--run_data_pipeline false` (COPAL's default path when
`msa_path` is set).

```bash
# 1. Clone + apply the 7-decimal ipTM patch in-place (same as step 1 in §3).
git clone https://github.com/google-deepmind/alphafold3.git $AF3_ROOT
cp third_party/alphafold3-patches/confidence_types.py \
   $AF3_ROOT/src/alphafold3/model/confidence_types.py

# 2. Build + install into a dedicated venv.
cd $AF3_ROOT
uv venv .venv --python 3.12 && source .venv/bin/activate
uv pip install -e .

# 3. Build AF3's chemical-component constants. The Dockerfile does this
# automatically (`RUN uv run build_data`); pip install does not, and AF3
# will crash at import time with a missing chemical_component_sets.pickle
# if you skip it.
build_data

# 4. Verify.
python -c "import alphafold3; import jax; print('jax devices:', jax.devices())"
```

> **Troubleshooting — broken pip-shim `cmake`.** If `cmake --version`
> errors with a Python `Traceback` (common symptom: you once ran
> `pip install cmake` into your user site), scikit-build-core will
> pick up that broken shim and `uv pip install -e .` will fail. Either
> `pip uninstall cmake` globally, or front-load the real system cmake
> for the build: `PATH=/usr/bin:$PATH uv pip install -e .`.

Then point COPAL at the matching wrapper (`external/wrappers/af3_native.sh`
in this repo), which `source`s `$AF3_ROOT/.venv/bin/activate` and calls
`run_alphafold.py` directly — no apptainer / docker:

```bash
export COPAL_AF3_SCRIPT=$(pwd)/external/wrappers/af3_native.sh
```

> **Performance caveat (tokamax autotune).** The `tokamax` wheel ships
> autotuning caches for H100 and B200, but **not for A100**. The
> official container bakes an A100 cache in during build; the native
> install does not. First AF3 inference on an A100 therefore spends a
> few extra minutes autotuning Pallas/Triton kernels — e.g. a 256-token
> complex measured ~1m30s in the IPD container vs ~4m50s native first
> run. The cache is per-process only (not persisted to `~/.cache/…`),
> but COPAL runs all predictions in a chunk in one process, so the
> penalty is paid once per `copal predict` call.

> **Numerical note.** Inference on CUDA-13 host libraries vs the
> CUDA-12.6 container can differ at the ~10⁻⁴ level on ipTM (e.g.
> `0.9702934` container vs `0.9701854` native on the same input). Fine
> for ranking; not bit-exact, so don't mix the two for reproducibility
> studies.

## 4. Generating MSAs for your own sequences

All three models accept pre-computed `.a3m` multiple-sequence alignments
(AF3 via `unpairedMsaPath`, Boltz2 via `msa`, RF3 via `msa_path`) — COPAL's
per-model input builders wire those through from the `msa_path` column of
your receptors CSV.

Use the ColabFold MMseqs2 cloud API (free, zero local setup, ~30 s per
sequence for typical protein lengths; rate-limited, so keep it civil):

```bash
uv pip install colabfold

# One FASTA with all your receptors; colabfold_search emits one .a3m per entry.
colabfold_search my_receptors.fasta out_msas/ --use-env 1 --db-load-mode 2
```

Then point each row's `msa_path` in your receptors CSV at the matching
`out_msas/<id>.a3m`.

The 281 MSAs shipped in `data/benchmark_msas.tar.xz` were produced with
this same pipeline. For larger runs (≳1000 sequences) or fully-offline /
bit-for-bit-reproducible operation, swap to [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)
(same `colabfold_search` CLI, against a local ~220 GB UniRef30 +
ColabFoldDB mirror).

## 5. Boltz2

Install Boltz2 into its own uv venv under `external/boltz/`:

```bash
uv venv external/boltz/.venv --python 3.12 && source external/boltz/.venv/bin/activate
uv pip install boltz
```

Then point COPAL at the venv's `boltz` entry point:

```bash
export BOLTZ2_BIN=<repo>/external/boltz/.venv/bin/boltz
```

Boltz2 weights download automatically on first run. License: MIT.

> **CUDA-version mismatch on recent Boltz releases.** Boltz 2.2.1's
> pypi wheel pulls `torch==2.11.0+cu130`, but its required
> `cuequivariance-ops-torch-cu12` kernel is built against CUDA 12 —
> loading it at runtime fails with
> `undefined symbol: cublasGemmGroupedBatchedEx, version libcublas.so.12`.
> Workaround: invoke `boltz predict --no_kernels` (documented upstream
> flag that bypasses the custom triangular-mult kernel; ~1.5–2x slower
> but numerically identical). The bundled `boltz2.sh` wrapper does not
> set this; copy it, add `--no_kernels`, and point
> `COPAL_BOLTZ2_SCRIPT` at your copy. An example is checked in at
> `external/wrappers/boltz2.sh`.

## 6. RosettaFold3 (Foundry)

```bash
git clone -b production https://github.com/RosettaCommons/foundry.git $RF3_INSTALL
cd $RF3_INSTALL
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[rf3]"
```

The `[rf3]` extras pull in `cuequivariance*`, which RF3 `fold` requires
at runtime — plain `uv pip install -e .` will install but blow up the
first time you call `rf3 fold`.

Download the RF3 checkpoint via the `foundry` CLI — this fetches the
exact ckpt matching the cloned branch (important: training-time
hyperparameters are baked into the ckpt, so a ckpt from a divergent
branch will fail to load with `TypeError: FabricTrainer.__init__() got
an unexpected keyword argument ...`):

```bash
foundry install rf3 --checkpoint-dir $RF3_INSTALL/checkpoints
ln -s checkpoints/rf3_foundry_01_24_latest_remapped.ckpt $RF3_INSTALL/rf3_latest.pt
export RF3_CKPT=$RF3_INSTALL/rf3_latest.pt
```

Run `foundry list-available` to see all registered checkpoints. If
`foundry` is not on PATH after `uv pip install -e ".[rf3]"`, run
`uv pip install -e .` once (plain, no extras) to register the
`foundry_cli` entry point, then retry.

License: BSD-3-Clause.

## 7. Shell wrappers (optional — you already have them)

Out of the box `copal predict` uses three bundled shell templates at
`src/copal/runner/scripts/{af3,boltz2,rf3}.sh`. They're documented,
generic, and assume an Apptainer/Singularity + `uv` install. If your
cluster needs something custom (Slurm module loads, a different launcher,
extra bind mounts), copy the relevant template to a writable location,
edit it, and point `copal` at it:

```bash
export COPAL_AF3_SCRIPT=$HOME/copal-scripts/af3.sh
export COPAL_BOLTZ2_SCRIPT=$HOME/copal-scripts/boltz2.sh
export COPAL_RF3_SCRIPT=$HOME/copal-scripts/rf3.sh
```

See `src/copal/runner/scripts/README.md` for the contract (input/output
directory shape) your wrapper must satisfy.

> **Bundled AF3 wrapper caveats.** `af3.sh` (a) omits `--model_dir`, so
> AF3 falls back to its `$HOME/models` default, and (b) hardcodes
> `--run_data_pipeline true`, which re-runs HMMer searches even when
> the input JSON already supplies `unpairedMsaPath`. For a typical COPAL
> run with precomputed MSAs, copy `af3.sh`, add
> `--model_dir "$AF3_MODEL_DIR"` and flip `--run_data_pipeline false`,
> then point `COPAL_AF3_SCRIPT` at your copy. An example is checked in
> at `external/wrappers/af3.sh` in this repo.

## 8. BLAST+ (optional)

Required only for `copal seqid-filter` and the RALPH BLAST baseline.

```bash
conda install -c bioconda blast        # easiest
# or:
sudo apt install ncbi-blast+
```

## 9. Extract the MSA archive (for benchmark reproduction)

```bash
mkdir -p data/benchmark_msas
tar -xJf data/benchmark_msas.tar.xz -C data/benchmark_msas/
```

Only needed if you want to re-run inference on the paper benchmark.
The `msa_path` column in `data/paper_luxrs.csv`
points at `data/benchmark_msas/<N>.a3m` (relative to the repo root),
so run `copal` commands from there.

## 10. Verify

```bash
uv run copal --help
uv run copal predict --help
uv run copal analyze --help
uv run copal full-ranking --help
uv run copal sim-matrix --help
uv run python -c "import copal, copal.pipeline, copal.analysis"
```

If those commands succeed, the local install is ready.

### Per-model venv smoke tests

Run each of these from the repo root — all three are cheap (no GPU
allocation, no real inference):

```bash
# AF3: patch applied + container + DB + model-params reachable + import works
diff -q third_party/alphafold3-patches/confidence_types.py \
        external/alphafold3/src/alphafold3/model/confidence_types.py   # silent = OK
ls -la "$AF3_CONTAINER"                                                # .sif exists
ls "$AF3_DB_PATH" | head                                               # FASTAs visible
ls "$AF3_MODEL_DIR"/af3.bin* 2>/dev/null                               # model params present
apptainer exec --bind "$AF3_ROOT:$AF3_ROOT" "$AF3_CONTAINER" \
    python -c "import alphafold3, sys; print('af3', alphafold3.__version__)"

# Boltz2
external/boltz/.venv/bin/boltz --help
external/boltz/.venv/bin/python -c "import boltz; print(boltz.__version__)"

# RF3
external/foundry/.venv/bin/rf3 --help
readlink external/foundry/rf3_latest.pt                                # resolves
external/foundry/.venv/bin/python -c "import cuequivariance_torch; print(cuequivariance_torch.__version__)"
```
