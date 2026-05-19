# RALPH BLAST baseline

A homology-search baseline: for each query LuxR, BLAST against the RALPH
reference database and use the cognate-ligand annotations of the top BLAST
hits to produce a ranked ligand list.

## Dependencies

1. **BLAST+** (`blastp`, `makeblastdb`).

   - Conda: `conda install -c bioconda blast`
   - Debian/Ubuntu: `sudo apt install ncbi-blast+`
   - Or download official binaries from
     <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>.

   Verify with `which blastp`.

2. **The `ralph` Python package** — ships in-tree under
   [`src/ralph/`](../../src/ralph/) and is installed by `uv sync` from the
   repo root. No separate install needed. Provides
   `ralph.blastp.blastp_query`, `ralph.ranking.aggregate_rankings`, and
   `ralph.database.load_scores`.

3. **The `copal` package** (this repository) — `uv sync` from the repo root.

## Building the RALPH database

The RALPH database is a directory containing:

```
ralph_database/
├── ralph_scores.csv     # per-LuxR ranked AHL scores (training labels)
├── ralph.fasta          # reference LuxR sequences
├── ralph.phr            # ─┐
├── ralph.pin            #  │ BLAST DB files produced by makeblastdb
├── ralph.psq            # ─┘
└── (optional .pdb/.pog/.pos/.pot/.ptf/.pto extras)
```

To build it from a source FASTA + scores CSV:

```bash
mkdir -p baselines/ralph_blast/ralph_database
cd baselines/ralph_blast/ralph_database

# 1. Place ralph_scores.csv here (per-LuxR annotated AHL scores).
# 2. Place ralph.fasta here (one record per RALPH LuxR sequence).

makeblastdb \
    -in ralph.fasta \
    -dbtype prot \
    -out ralph
```

Where `ralph_scores.csv` is expected to have at least the columns
`luxr_name` and per-AHL score columns (consult `ralph.database.load_scores`
for the exact schema).

> If you don't have a RALPH source FASTA + scores file, contact the COPAL
> authors — the data accompanying the paper is hosted externally (see the
> repository root README.md for the DOI).

## Running

```bash
cd baselines/ralph_blast
./run_blast_analysis.sh
```

Or directly:

```bash
python run_ralph_blast.py \
    --benchmark-csv ../../data/paper_luxrs.csv \
    --ralph-database ./ralph_database \
    --ligands-path ../../data/curated_ligands_annotated.csv \
    --output-root ./outputs/results \
    --workflow-dir ./workflow_ralph \
    --min-pident 30 \
    --max-hits 50
```

Outputs land under `outputs/results/`:

- `ranking_performance/ranking_performance_conf{1,2,3}.csv`
- `full_rankings/all_splits.csv`
- `full_rankings/RALPH.xlsx`

Tune `--min-pident` to control which hits are kept and how the baseline
behaves on RefAHL splits.
