# RALPH vendor integration — design

**Date:** 2026-04-23
**Status:** approved (pending user review of this document)
**Scope:** integrate `~/projects/ralph` into the COPAL publication repo
(`/home/dnan/projects/luxr/src/publication_repo`) so that the existing
`baselines/ralph_blast/` baseline works out of a single `uv sync`, and ship
the Colab notebook at `~/projects/ralph/colab/RALPH_blastp_inference.ipynb`
as an in-tree application.

## Motivation

The repo already has `baselines/ralph_blast/run_ralph_blast.py` which imports
`ralph.blastp`, `ralph.ranking`, and `ralph.database` at module load, but
`ralph` is not declared anywhere in `pyproject.toml`. The current README
tells users to `git clone` ralph and `uv pip install -e` it separately. For
a publication repo this is brittle — the paper's numbers depend on a
specific ralph commit, so the repo should pin exactly that commit in-tree.

RALPH is ~30KB of pure Python with no model weights. Vendor scale is
appropriate; `external/` (sized for AF3/Boltz2/RF3, multi-GB each) is not.

## Non-goals

- Refactor of `ralph/alignment.py` (biotite) vs `src/copal/analysis/compute_seqid.py` (Biopython). Different purposes, stay separate.
- Contributing a LICENSE back to the upstream `github.com/davinan/ralph` repo. Out of scope here; should be done upstream independently.
- Changes to the ralph API surface.
- New tests for ralph. The vendored copy is a snapshot; `baselines/ralph_blast/` has no tests today, and the integration does not introduce them.

## Code-duplication audit

Compared `~/projects/ralph/ralph/*.py` against `src/copal/**/*.py`:

| Symbol | ralph location | copal location | Conclusion |
|---|---|---|---|
| `blastp_query` (1-vs-DB BLASTP) | `ralph/blastp.py` | — | unique to ralph |
| `blastp` subprocess (N×N seqid) | — | `src/copal/analysis/compute_seqid.py` | different purpose; keep both |
| `aggregate_rankings` (RRF over `ensemble_*` cols of scores table) | `ralph/ranking.py` | — | unique to ralph |
| `make_ensemble_ranking` (fuses per-model ranks AF3/Boltz2/RF3) | — | `src/copal/analysis/pipeline.py` | different inputs/outputs; shared `ensemble_rank`/`ensemble_score` names are an intentional contract |
| `load_ligands` (ralph_database/ligands_58.csv + RDKit) | `ralph/database.py` | — | — |
| `load_ligands_csv` (generic copal ligand CSV) | — | `src/copal/io/schema.py` | name collision only; namespacing (`ralph.*` vs `copal.io.*`) resolves it |
| `plot_sequence_alignment` (biotite BLOSUM62) | `ralph/alignment.py` | — | — |
| Biopython `substitution_matrices.load("BLOSUM62")` | — | `src/copal/analysis/compute_seqid.py` | two alignment libs in one repo, acceptable |

**No literal code duplication.** Namespacing under `src/ralph/` vs
`src/copal/` keeps the name-collision surface invisible to consumers.

## Design

### 1. Directory layout (delta)

```
src/ralph/                             # NEW — vendored pure-Python package
├── __init__.py                        # copied verbatim from ~/projects/ralph/ralph/
├── blastp.py                          # ditto
├── ranking.py                         # ditto
├── database.py                        # ditto
├── alignment.py                       # ditto
├── visualization.py                   # ditto
└── LICENSE                            # NEW — MIT, credits RALPH authors

applications/                          # NEW — top-level user-facing notebooks
└── RALPH_blastp_inference.ipynb       # adapted from ~/projects/ralph/colab/

docs/superpowers/specs/
└── 2026-04-23-ralph-vendor-integration-design.md    # this file
```

### 2. Packaging (`pyproject.toml`)

- `[tool.hatch.build.targets.wheel]` → `packages = ["src/copal", "src/ralph"]`
  (both packages ship from the one wheel).
- Add to `[project] dependencies`:
  - `matplotlib>=3.7` — used by `ralph.alignment` / `ralph.visualization`; not currently in copal deps.
  - `ipywidgets>=7.6` — used by `ralph.visualization` interactive widgets.
- Keep unchanged: `pandas`, `numpy`, `biotite`, `rdkit` — already present.
- No new `[project.scripts]` entry (ralph has no CLI).
- No new optional-dependencies extras required. (If notebook use turns out to feel heavy on fresh installs, a future `[notebook] = ["jupyter", "ipywidgets"]` extra is a drop-in — not part of this spec.)

### 3. Vendored source

Copy all six files from `~/projects/ralph/ralph/` to `src/ralph/` **byte-for-byte**. Internal imports (`from ralph.alignment import …`) continue to resolve because the package name is preserved.

Verify no modifications are needed: `ralph.database.load_ralph_database` currently hard-codes `/content/ralph_database` as a default and prints a "please download" hint. This is Colab-flavored but harmless — the baseline and the new notebook both pass `database_path` explicitly, so the default is never hit. Leave as-is for the vendor snapshot.

### 4. Notebook (`applications/RALPH_blastp_inference.ipynb`)

Adapt cell 1 of the upstream notebook. Cells 2–3 (blastp_query, interactive widget) are unchanged.

Cell 1 changes:
- **Remove** the `pip install -e /content/ralph` / `pip install git+…` branch. Replace with a one-line assertion `import ralph; print(ralph.__version__)`.
- **Remove** the `apt-get install ncbi-blast+` line. Replace with:
  ```python
  import shutil
  if shutil.which("blastp") is None:
      raise RuntimeError("blastp not on PATH — see INSTALLATION.md for BLAST+ install")
  ```
- **Remove** `gdown` + the hard-coded Google Drive file ID. Replace with:
  ```python
  import os
  database_path = os.environ.get(
      "RALPH_DATABASE",
      "../baselines/ralph_blast/ralph_database",
  )
  ```
- Add a markdown cell above cell 1 explaining: (a) run this notebook from `applications/`, (b) where the RALPH database comes from (Zenodo DOI placeholder, matching the rest of the repo's data-hosting conventions), (c) how to override the path with `RALPH_DATABASE`.

### 5. Licensing

- Add `src/ralph/LICENSE` — MIT text, copyright attributed to the RALPH authors (Davi Nakajima An et al.), same format as the repo root `LICENSE`.
- Append one line to the root `NOTICE` stating: `src/ralph/` vendored from `github.com/davinan/ralph` under MIT.
- The repo root `LICENSE` does not need changes — `src/ralph/LICENSE` governs the vendored subtree.

### 6. Docs updates

- `baselines/ralph_blast/README.md`:
  - Delete the "git clone ralph / `uv pip install -e /path/to/ralph`" section (lines ~18–27).
  - Keep: BLAST+ install instructions, RALPH database build, invocation example.
  - Add a one-liner: "`ralph` ships in-tree under `src/ralph/`; no separate install required."
- `INSTALLATION.md`: no structural changes. The existing "BLAST+ required only for `copal seqid-filter` and the RALPH BLAST baseline" line (:382) remains accurate.
- `README.md` (root): no mandatory change. Optional: one-line pointer under a quickstart to the new `applications/RALPH_blastp_inference.ipynb`.

### 7. What does NOT change

- `baselines/ralph_blast/run_ralph_blast.py` — `from ralph.blastp import …` imports resolve once `src/ralph/` exists in the package.
- `baselines/ralph_blast/run_blast_analysis.sh` — unchanged.
- `src/copal/**` — untouched.
- `external/` directory — ralph does not go here.
- The RALPH database build process (still `makeblastdb` from `ralph.fasta` under `baselines/ralph_blast/ralph_database/`).

## Risks / open questions

1. **Upstream ralph divergence.** After this vendor, edits made to `~/projects/ralph/` do not automatically flow into the publication repo. Acceptable for a frozen-at-submission publication repo; worth a note in `baselines/ralph_blast/README.md` so future contributors know to re-sync intentionally if ralph is still evolving.
2. **Notebook-execution CI.** No existing notebook-execution tests in the repo. This spec does not add any; the notebook is a user-facing explorer, not a test artifact. If that changes, a `nbval` or `papermill` smoke test would be the right addition — out of scope here.
3. **`ralph.database.load_ralph_database` Colab defaults.** Kept as-is (see §3). If a user in a non-Colab environment calls it without `database_path`, they'll get a misleading "/content/..." printout. Low risk; the documented usage (baseline script + new notebook) never hits this code path. Fix is a trivial follow-up if it bites.

## Verification

After implementation, the integration is verified by:

1. `uv sync` from the repo root completes without network access to any external ralph repo.
2. `python -c "import ralph; print(ralph.__version__)"` prints `0.1.0`.
3. `baselines/ralph_blast/run_blast_analysis.sh` runs end-to-end against a built `ralph_database/` and produces the documented outputs (`ranking_performance_conf{1,2,3}.csv`, `full_rankings/all_splits.csv`, `full_rankings/RALPH.xlsx`).
4. `applications/RALPH_blastp_inference.ipynb` executes top-to-bottom in Jupyter with a `RALPH_DATABASE` env var set, and the interactive widget renders.

Step (3) is the hard verification — it was the original reason ralph had to be importable at all.
