# RALPH BLAST baseline

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/davinan/COPAL/blob/main/applications/RALPH_blastp_inference.ipynb)

A homology-search baseline / inference tool: BLAST a query LuxR against the
RALPH reference database and read off the cognate-ligand rankings of the top
hits. The interactive entry point is
[`applications/RALPH_blastp_inference.ipynb`](../../applications/RALPH_blastp_inference.ipynb)
— click the badge above to run it in Google Colab (a fresh runtime in your own
Google account), or run it locally.

## Dependencies

1. **BLAST+** (`blastp`, `makeblastdb`).
   - Conda: `conda install -c bioconda blast`
   - Debian/Ubuntu: `sudo apt install ncbi-blast+`
   - Or the official static binaries from
     <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>.

   Verify with `which blastp makeblastdb`.

2. **The `ralph` + `copal` packages** — ship in-tree; `uv sync` (or
   `pip install -e .`) from the repo root installs them (RDKit + ipywidgets
   included).

## The database (ships in-repo, compressed)

`ralph_database/` contains the database as **compressed sources** (git-committed):

```
ralph_database/
├── ralph_scores.csv.xz        # wide rankings: luxr_name, source, organism,
│                              #   true_ahl_id, mpKT, sequence, {model}_1..58
│                              #   (AHL id at each rank) for model in
│                              #   ensemble/af3/rf3/boltz2/boltz2affinity
├── ralph_model_scores.csv.xz  # per-pair scores by (luxr_name, ahl):
│                              #   copal, af3, rf3, boltz2, boltz2affinity
├── ralph_db.faa.xz            # BLAST FASTA (headers = luxr_name)
└── ligands_58.csv             # AHLid, abbreviation, smi
```

It covers **10,496 RALPH homologs + 153 RefAHL benchmark receptors** (10,649
sequences). RefAHL rows carry the experimentally-known cognate (`true_ahl_id`).

### Setup (one command)

```bash
bash baselines/ralph_blast/setup_db.sh
```

This decompresses the `.xz` sources to plain `.csv`/`.faa` and runs
`makeblastdb` to build the BLAST index. It is idempotent (safe to re-run) and
the materialized files are git-ignored. The notebook auto-runs it if the index
is missing.

### Rebuilding the shipped sources (maintainers only)

The compressed sources are produced by `build_ralph_db.py` from the published
supplementary data (RALPH rankings/scores workbook + RefAHL benchmark rankings +
`ralph_v3.csv` sequences); those inputs are external to the repo. End users do
**not** run this — they use the shipped `.xz` + `setup_db.sh`.

```bash
python build_ralph_db.py \
    --xlsx <Supp_Table6_RALPH_rankings.xlsx> \
    --ralph-v3 <ralph_v3.csv> \
    --bench-rankings <benchmark_data_all_rankings.csv> \
    --bench-ensemble <benchmark_data_all_ensemble_rankings.csv>
```

(run from the repo root; `--ligands` / `--paper-luxrs` default to `data/`.)

## Running the CLI baseline (benchmark)

```bash
python run_ralph_blast.py \
    --benchmark-csv ../../data/paper_luxrs.csv \
    --ralph-database ./ralph_database \
    --ligands-path ../../data/curated_ligands_annotated.csv \
    --output-root ./outputs/results \
    --min-pident 30 --max-hits 50
```

> **Note.** For a leakage-controlled *benchmark*, point `--ralph-database` at a
> leakage-filtered, RALPH-only database (see the paper). The shipped
> `ralph_database/` is the *inference* database: full, unfiltered, and includes
> the RefAHL receptors so a query can hit either — a different use case.

Outputs land under `outputs/results/` (`ranking_performance/`, `full_rankings/`).
