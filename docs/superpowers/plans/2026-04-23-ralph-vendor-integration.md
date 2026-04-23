# RALPH vendor integration — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Vendor the `ralph` Python package from `~/projects/ralph` into the COPAL publication repo at `src/ralph/`, add the user-facing Colab notebook at `applications/RALPH_blastp_inference.ipynb`, and make everything install in one `uv sync`.

**Architecture:** In-tree vendored snapshot. Six Python files copied byte-for-byte from `~/projects/ralph/ralph/` into `src/ralph/`. `pyproject.toml` grows by one package and two dependencies. A local-friendly notebook replaces the Colab setup cell. License/NOTICE updated. No changes to `src/copal/`; no changes to `baselines/ralph_blast/run_ralph_blast.py` (its top-level `from ralph.* import …` imports start resolving once `src/ralph/` exists).

**Tech Stack:** Python ≥3.10, `uv` + `hatchling`, Jupyter. Runtime deps already in place: `pandas`, `numpy`, `biotite`, `rdkit`. New deps: `matplotlib`, `ipywidgets`. System dep: `blastp` (already required for existing `copal seqid-filter`).

**Repo context:**
- This publication repo is **not a git repository** (`git rev-parse --is-inside-work-tree` fails). Skip every "commit" step below unless the engineer first `git init`s the repo at the user's request. Each task still has a logical stopping point for review.
- Reference spec: `docs/superpowers/specs/2026-04-23-ralph-vendor-integration-design.md`.
- Source of truth being vendored: `~/projects/ralph/` (git: `github.com/davinan/ralph.git`, currently at commit `92b226a` — record this in the README so future contributors know exactly what was copied).

---

## File structure overview

**Created:**
- `src/ralph/__init__.py`            (copy of `~/projects/ralph/ralph/__init__.py`, 28 lines)
- `src/ralph/alignment.py`           (copy, 90 lines; biotite BLASTDB + alignment plotting)
- `src/ralph/blastp.py`              (copy, 77 lines; `blastp` subprocess wrapper)
- `src/ralph/database.py`            (copy, 86 lines; scores/ligands loaders)
- `src/ralph/ranking.py`             (copy, 86 lines; RRF-style `aggregate_rankings`)
- `src/ralph/visualization.py`       (copy, 441 lines; interactive widgets + clustered rank plot)
- `src/ralph/LICENSE`                (new MIT LICENSE, copyright to RALPH authors)
- `applications/RALPH_blastp_inference.ipynb` (adapted from `~/projects/ralph/colab/RALPH_blastp_inference.ipynb`)

**Modified:**
- `pyproject.toml`                   (add package to wheel targets; add 2 deps)
- `NOTICE`                           (append one stanza about `src/ralph/`)
- `baselines/ralph_blast/README.md`  (drop "clone + pip install ralph" section)

**Untouched:**
- `src/copal/**`
- `baselines/ralph_blast/run_ralph_blast.py`
- `baselines/ralph_blast/run_blast_analysis.sh`
- `INSTALLATION.md`, `README.md` (root)

---

## Task 1: Vendor the ralph source tree

**Files:**
- Create: `src/ralph/__init__.py`
- Create: `src/ralph/alignment.py`
- Create: `src/ralph/blastp.py`
- Create: `src/ralph/database.py`
- Create: `src/ralph/ranking.py`
- Create: `src/ralph/visualization.py`

- [ ] **Step 1: Create the destination directory**

From the publication repo root:

```bash
mkdir -p src/ralph
```

- [ ] **Step 2: Copy all six ralph source files verbatim**

```bash
cp ~/projects/ralph/ralph/__init__.py     src/ralph/__init__.py
cp ~/projects/ralph/ralph/alignment.py    src/ralph/alignment.py
cp ~/projects/ralph/ralph/blastp.py       src/ralph/blastp.py
cp ~/projects/ralph/ralph/database.py     src/ralph/database.py
cp ~/projects/ralph/ralph/ranking.py      src/ralph/ranking.py
cp ~/projects/ralph/ralph/visualization.py src/ralph/visualization.py
```

- [ ] **Step 3: Verify byte-for-byte equality**

Run:

```bash
for f in __init__.py alignment.py blastp.py database.py ranking.py visualization.py; do
  diff -q ~/projects/ralph/ralph/$f src/ralph/$f || echo "MISMATCH: $f"
done
```

Expected: no output (no diffs, no MISMATCH lines).

- [ ] **Step 4: Confirm the line-count inventory matches the spec**

Run:

```bash
wc -l src/ralph/*.py
```

Expected (must match exactly):

```
    90 src/ralph/alignment.py
    77 src/ralph/blastp.py
    86 src/ralph/database.py
    28 src/ralph/__init__.py
    86 src/ralph/ranking.py
   441 src/ralph/visualization.py
   808 total
```

- [ ] **Step 5: Record the upstream commit**

Write the upstream SHA into a plain-text provenance note that future readers of `src/ralph/` will find. Append this block to `src/ralph/__init__.py` *above* the existing module docstring is **not** acceptable (we want the file byte-for-byte identical). Instead, create a new file:

Create: `src/ralph/VENDOR.md`

```markdown
# Vendored RALPH source

This directory is a byte-for-byte snapshot of
<https://github.com/davinan/ralph> at commit `92b226a`, imported into this
publication repository on 2026-04-23.

Do not modify files in this directory as a substitute for upstream fixes.
To re-sync with a newer upstream commit, re-copy all six `.py` files and
update the commit hash on this line.
```

- [ ] **Step 6: Checkpoint — list the new tree and confirm**

```bash
ls -la src/ralph/
```

Expected: `__init__.py`, `alignment.py`, `blastp.py`, `database.py`, `ranking.py`, `visualization.py`, `VENDOR.md` (and the `LICENSE` added in Task 2 not yet present).

No commit step — repo is not a git repository.

---

## Task 2: Add MIT LICENSE for the vendored subtree, update root NOTICE

**Files:**
- Create: `src/ralph/LICENSE`
- Modify: `NOTICE` (append one stanza)

- [ ] **Step 1: Write `src/ralph/LICENSE`**

Create the file with this exact content (same template as the root `LICENSE`, just rewording the copyright line to attribute RALPH):

```
MIT License

Copyright (c) 2026 The RALPH authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

- [ ] **Step 2: Append the RALPH stanza to `NOTICE`**

The current `NOTICE` has a "Third-party components" section listing the `third_party/alphafold3-patches/confidence_types.py` entry. Add a second bullet immediately under that one (before the "External models invoked at runtime" section). After editing, the "Third-party components" section must read:

```
Third-party components
----------------------

`third_party/alphafold3-patches/confidence_types.py`
    Derived from AlphaFold 3 source code by Google DeepMind. Distributed
    under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0
    International License (CC-BY-NC-SA 4.0). See
    `third_party/alphafold3-patches/LICENSE` for the full license text and
    `third_party/alphafold3-patches/README.md` for provenance and usage.

`src/ralph/`
    Vendored snapshot of the RALPH package (https://github.com/davinan/ralph),
    redistributed under the MIT License. See `src/ralph/LICENSE` for the
    full license text and `src/ralph/VENDOR.md` for the exact upstream
    commit.
```

The rest of `NOTICE` (copyright line at top, "External models invoked at runtime" section) is unchanged.

- [ ] **Step 3: Verify both files**

```bash
head -1 src/ralph/LICENSE
grep -c "src/ralph/" NOTICE
```

Expected:
- `head -1 src/ralph/LICENSE` → `MIT License`
- `grep -c "src/ralph/" NOTICE` → `1` (or higher)

---

## Task 3: Wire the new package into `pyproject.toml`

**Files:**
- Modify: `pyproject.toml`

Context: current pyproject has one package (`src/copal`) under `[tool.hatch.build.targets.wheel]` and a flat `dependencies` list. Two small edits.

- [ ] **Step 1: Add `src/ralph` to the wheel packages list**

Find the existing block:

```toml
[tool.hatch.build.targets.wheel]
packages = ["src/copal"]
```

Replace with:

```toml
[tool.hatch.build.targets.wheel]
packages = ["src/copal", "src/ralph"]
```

Leave the `[tool.hatch.build.targets.wheel.force-include]` block that follows it unchanged.

- [ ] **Step 2: Add the two new runtime dependencies**

Find the existing `dependencies = [ ... ]` list. It currently ends with:

```toml
    "openpyxl>=3.1",
]
```

Insert two new entries before the closing `]`, preserving alphabetical order would normally suggest placing `matplotlib` after `gemmi` and `ipywidgets` after... but the existing list is not alphabetized (it's grouped by role), so **just append** both entries at the end for minimal diff. After editing the list must end like:

```toml
    "openpyxl>=3.1",
    "matplotlib>=3.7",
    "ipywidgets>=7.6",
]
```

- [ ] **Step 3: Verify syntax by re-parsing with uv**

```bash
uv lock --check 2>&1 | head -40
```

A clean run will either accept the existing lockfile (if the new deps happen to already be satisfied by the current environment) or announce that the lock is out of date. Either is fine at this step — the important thing is that the TOML parses. If you see a parse error, re-read the file and fix.

---

## Task 4: Re-sync the environment and smoke-test imports

**Files:**
- Modify: `uv.lock` (regenerated)

- [ ] **Step 1: Regenerate the lockfile and sync**

```bash
uv sync
```

Expected: completes without errors; installs `matplotlib`, `ipywidgets`, and any of their transitive dependencies that were not already in the environment.

- [ ] **Step 2: Smoke-test `import ralph`**

```bash
uv run python -c "import ralph; print(ralph.__version__)"
```

Expected: prints `0.1.0`.

- [ ] **Step 3: Smoke-test the specific symbols the baseline imports**

```bash
uv run python -c "from ralph.blastp import blastp_query; from ralph.ranking import aggregate_rankings; from ralph.database import load_scores; print('ok')"
```

Expected: prints `ok`. (These are the three imports at the top of `baselines/ralph_blast/run_ralph_blast.py`.)

- [ ] **Step 4: Smoke-test that `copal` still imports**

```bash
uv run python -c "import copal; import copal.cli; import copal.analysis.pipeline; print('ok')"
```

Expected: prints `ok`. (Guards against accidental breakage of the existing package.)

- [ ] **Step 5: Smoke-test the baseline script's `--help`**

```bash
uv run python baselines/ralph_blast/run_ralph_blast.py --help
```

Expected: prints the argparse help text documented in `baselines/ralph_blast/README.md` (flags `--benchmark-csv`, `--ralph-database`, etc.). This proves the `from ralph.* import …` statements at the top of the script resolve.

If any of Steps 2–5 fail, do not proceed to Task 5. The failure is either a missing file (re-check Task 1 step 3's diff output) or a pyproject misconfiguration (re-check Task 3).

---

## Task 5: Adapt the Colab notebook to `applications/`

**Files:**
- Create: `applications/RALPH_blastp_inference.ipynb`

Context: the upstream Colab notebook has 4 cells. Cells 2 and 3 (blastp_query invocation, interactive widget) are reused verbatim. Cell 1 is rewritten entirely. A new markdown cell 0 is added. The empty trailing cell is dropped.

- [ ] **Step 1: Create the target directory**

```bash
mkdir -p applications
```

- [ ] **Step 2: Start from a copy of the upstream notebook**

```bash
cp ~/projects/ralph/colab/RALPH_blastp_inference.ipynb applications/RALPH_blastp_inference.ipynb
```

- [ ] **Step 3: Edit cell 0 (currently the query-sequence cell)**

Leave this cell **unchanged**. It defines `query_sequence` and `top_k_BLASTP_hits`. This was called "cell-0" in the upstream notebook.

- [ ] **Step 4: Replace cell 1 (the setup cell) entirely**

The upstream cell 1 runs `pip install -e /content/ralph` / `pip install git+…`, `apt-get install ncbi-blast+`, and `gdown` to fetch the database. None of that applies in-repo. Replace the **entire cell 1 source** with:

```python
#@title Environment check and database path

import os
import shutil

# --- ralph availability ---
# ralph ships in this repository under src/ralph/ and is installed by
# `uv sync` from the repo root. If this import fails, run `uv sync`.
from ralph import (
    blastp_query,
    aggregate_rankings,
    clustered_rank_plot,
    load_scores,
    load_ligands,
    interactive_blastp_results,
    interactive_alignment_viewer,
)
import pandas as pd
import numpy as np

# --- BLAST+ availability ---
if shutil.which("blastp") is None:
    raise RuntimeError(
        "blastp not on PATH. Install BLAST+ and re-run. "
        "See INSTALLATION.md in the repo root."
    )

# --- Database path ---
# Default: the baseline's working directory relative to this notebook.
# Override with the RALPH_DATABASE environment variable if your database
# lives elsewhere.
database_path = os.environ.get(
    "RALPH_DATABASE",
    os.path.abspath(os.path.join(os.getcwd(), "..", "baselines", "ralph_blast", "ralph_database")),
)
if not os.path.isdir(database_path):
    raise FileNotFoundError(
        f"RALPH database directory not found at {database_path}. "
        f"See baselines/ralph_blast/README.md for build instructions, or "
        f"set RALPH_DATABASE to an existing database directory."
    )

df = load_scores(database_path=database_path)
ligand_df = load_ligands(database_path=database_path)
```

- [ ] **Step 5: Edit cell 2 (BLASTP query cell)**

Change the single line:

```python
database_path = "/content/ralph_database"
```

to:

```python
# database_path is already defined in the setup cell above.
```

The rest of the cell (the `db_path = os.path.join(...)`, `hits = blastp_query(...)`, `hits = hits.sort_values(...)`, `scores = df[...]`) is unchanged.

- [ ] **Step 6: Leave cell 3 (interactive widget cell) unchanged**

Confirm the cell still contains:

```python
interactive_blastp_results(
    hits_df=hits,
    scores_df=df,
    query_sequence=query_sequence,
    ligand_df=ligand_df
)
```

- [ ] **Step 7: Remove the trailing empty cell 4**

The upstream notebook has an empty final cell. Delete it.

- [ ] **Step 8: Insert a new markdown cell at position 0 (before the query-sequence cell)**

With this source:

```markdown
# RALPH BLASTP inference

This notebook is the user-facing entry point for the RALPH BLAST baseline.
Given a LuxR-homolog query sequence, it runs BLASTP against the RALPH
database and visualizes the top-ranked ligands using the ENSEMBLE method.

## Before running

1. From the repo root, run `uv sync` to install the `ralph` package and
   its dependencies (ships in-tree under `src/ralph/`).
2. Install BLAST+ (`blastp` on PATH). See `INSTALLATION.md`.
3. Build the RALPH database. See `baselines/ralph_blast/README.md`. By
   default this notebook looks for the database at
   `../baselines/ralph_blast/ralph_database` relative to this notebook,
   or at the path given by the `RALPH_DATABASE` environment variable.

## Running

Change the `query_sequence` variable in the next cell, then run all cells.
```

- [ ] **Step 9: Validate the notebook JSON**

```bash
uv run python -c "import json; nb = json.load(open('applications/RALPH_blastp_inference.ipynb')); cells = nb['cells']; print(f'{len(cells)} cells:'); [print(f'  [{i}] {c[\"cell_type\"]}') for i, c in enumerate(cells)]"
```

Expected output:

```
5 cells:
  [0] markdown
  [1] code
  [2] code
  [3] code
  [4] code
```

- [ ] **Step 10: Start Jupyter and run-all**

```bash
uv run jupyter nbconvert --to notebook --execute applications/RALPH_blastp_inference.ipynb --output /tmp/ralph-smoke.ipynb
```

Notes:
- If no RALPH database is available locally, this will fail at cell 1 with the "RALPH database directory not found" error — that's the *expected* failure mode, and it demonstrates cell 1's guard works.
- If a database *is* available (e.g., set `RALPH_DATABASE` first), the cell should run end-to-end. The interactive widget cell will render without showing the widget in the `nbconvert` output — that's fine; it only works in a live kernel.

If the failure mode is anything other than (a) success or (b) the expected "database not found" `FileNotFoundError`, stop and debug.

---

## Task 6: Update docs — baseline README

**Files:**
- Modify: `baselines/ralph_blast/README.md`

- [ ] **Step 1: Locate the "Dependencies" section**

The current file has a numbered Dependencies list with three entries: (1) BLAST+, (2) the `ralph` Python package + manual install instructions, (3) the `copal` package. Entry (2) is the one to change.

- [ ] **Step 2: Replace entry (2)**

Find:

```markdown
2. **The `ralph` Python package** (provides `ralph.blastp.blastp_query`,
   `ralph.ranking.aggregate_rankings`, `ralph.database.load_scores`).
   Install from source:

   ```bash
   git clone <ralph-repo-url> /path/to/ralph
   uv pip install -e /path/to/ralph
   ```

   Or add it to a shared environment.
```

Replace with:

```markdown
2. **The `ralph` Python package** — ships in-tree under
   [`src/ralph/`](../../src/ralph/) and is installed by `uv sync` from the
   repo root. No separate install needed. Provides
   `ralph.blastp.blastp_query`, `ralph.ranking.aggregate_rankings`, and
   `ralph.database.load_scores`.
```

- [ ] **Step 3: Verify the rest of the README still reads coherently**

```bash
grep -n -i "ralph" baselines/ralph_blast/README.md
```

Check that every remaining mention of `ralph` refers to either the *database* (`ralph_database/`, `ralph_scores.csv`, `ralph.fasta`), the *package functions* (`ralph.database.load_scores`), or the in-tree *package location* (`src/ralph/`). There should be no remaining instruction to `git clone` or `pip install` ralph.

- [ ] **Step 4: Spot-check rendered links**

In the replacement, `../../src/ralph/` is a relative link from `baselines/ralph_blast/README.md` to the repo root's `src/ralph/`. Verify the path resolves:

```bash
( cd baselines/ralph_blast && ls ../../src/ralph/ )
```

Expected: lists the six `.py` files, `LICENSE`, `VENDOR.md`.

---

## Task 7: Final end-to-end verification

**Files:** none (verification only).

- [ ] **Step 1: Fresh install simulation**

```bash
rm -rf .venv && uv sync
```

Expected: clean install; no missing deps; no reference to `~/projects/ralph`.

- [ ] **Step 2: Re-run all smoke tests from Task 4**

```bash
uv run python -c "import ralph; print(ralph.__version__)"
uv run python -c "from ralph.blastp import blastp_query; from ralph.ranking import aggregate_rankings; from ralph.database import load_scores; print('ok')"
uv run python -c "import copal; print('ok')"
uv run python baselines/ralph_blast/run_ralph_blast.py --help | head -5
```

All four should succeed as in Task 4.

- [ ] **Step 3: If a RALPH database is available locally, end-to-end baseline run**

This is the paper-level verification. Only feasible on a machine that has
a built `ralph_database/` with `ralph_scores.csv`, `ralph.fasta`, and the
`makeblastdb`-produced `.phr`/`.pin`/`.psq` files. If the database is not
available, skip this step and flag it to the user.

```bash
cd baselines/ralph_blast
./run_blast_analysis.sh
```

Expected: produces `outputs/results/ranking_performance/ranking_performance_conf{1,2,3}.csv`, `outputs/results/full_rankings/all_splits.csv`, and `outputs/results/full_rankings/RALPH.xlsx`. This proves the vendor integration has not changed the baseline's numerical output.

- [ ] **Step 4: Report status to the user**

Summarize in one message: what changed, which verifications passed, and whether Step 3 (end-to-end baseline) was actually run or skipped for lack of a local database.

---

## Self-review notes

- **Spec coverage:** every section of the spec maps to a task — §1 layout → Task 1; §2 packaging → Task 3; §3 vendored source → Task 1; §4 notebook → Task 5; §5 licensing → Task 2; §6 docs → Task 6; §7 untouched files (verified by Task 4 step 4 / Task 7 step 2).
- **Placeholders:** no TBDs, no "implement later" — every code/config block is a complete replacement.
- **Type consistency:** package name `ralph` used consistently; symbol names `blastp_query`, `aggregate_rankings`, `load_scores`, `load_ligands`, `interactive_blastp_results`, `interactive_alignment_viewer`, `clustered_rank_plot` match `~/projects/ralph/ralph/__init__.py` exactly.
- **TDD note:** the spec explicitly scopes out new tests for ralph. Verification is done via smoke imports (Task 4) and end-to-end baseline run (Task 7 step 3) rather than unit tests. This matches the current testing posture of `baselines/ralph_blast/` (which has none).
- **Git note:** the repo is not a git repository; "commit" steps are omitted. If the user later runs `git init`, each task's stopping point is a natural commit boundary.
