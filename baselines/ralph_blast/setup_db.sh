#!/usr/bin/env bash
# Materialize the RALPH-BLAST inference database:
#   1. decompress the shipped .xz sources to plain .csv / .faa
#   2. build the BLAST protein index with makeblastdb
#
# Idempotent and safe to re-run: existing outputs are skipped. The decompressed
# files and BLAST index are git-ignored (only the .xz sources + ligands_58.csv
# are committed). Run once from anywhere; it resolves paths relative to itself.
set -euo pipefail

DB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/ralph_database" && pwd)"
cd "$DB_DIR"

if ! command -v makeblastdb >/dev/null 2>&1; then
    echo "ERROR: makeblastdb not found on PATH. Install BLAST+ first" >&2
    echo "  (conda install -c bioconda blast  |  apt install ncbi-blast+)." >&2
    echo "  See INSTALLATION.md." >&2
    exit 1
fi

for stem in ralph_scores ralph_model_scores; do
    if [ ! -f "$stem.csv" ]; then
        echo "Decompressing $stem.csv.xz ..."
        xz -dk "$stem.csv.xz"
    fi
done

if [ ! -f ralph_db.faa ]; then
    echo "Decompressing ralph_db.faa.xz ..."
    xz -dk ralph_db.faa.xz
fi

if [ ! -f ralph.pin ]; then
    echo "Building BLAST index (makeblastdb) ..."
    makeblastdb -in ralph_db.faa -dbtype prot -out ralph
else
    echo "BLAST index already present (ralph.pin) - skipping makeblastdb."
fi

echo "RALPH database ready at: $DB_DIR"
