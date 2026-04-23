# `data/`

Small input files used by the analysis pipeline and shipped with the
repository for full reproducibility of the paper benchmark.

| File | Size | Purpose |
| ---- | ---- | ------- |
| `benchmark_data_all_corrected.csv` | ~200 KB | Per-receptor labels: `refahl_id`, `seq`, `true_ahl_id`, `true_ahl_abbrev`, `split` (`RefAHL`, `wetlabexp`, `solo`, …), `refahl_confidence_rating`, `msa_path`. Used by `copal analyze --labels` and as a receptors CSV for `copal predict`. `msa_path` values are repo-root-relative (`data/benchmark_msas/*.a3m`); extract `benchmark_msas.tar.xz` before running inference (see below). |
| `benchmark_msas.tar.xz` | 31 MB | ColabFold-style `.a3m` MSAs for all 281 benchmark receptors (one file per `refahl_id`, keyed by the number in the `msa_path` column). Extract with `tar -xJf data/benchmark_msas.tar.xz -C data/benchmark_msas/`. |
| `benchmark_ligands.csv` | 4 KB | Benchmark ligand subset used at inference time (`AHLid`, `smiles`, …). |
| `Final_confirmed_ligands_58_annotated.csv` | 5 KB | Annotated AHL ligand metadata (58 entries; `AHLid`, `abbreviation`, `Formula`, `smi`, `subcategory`). Used by `copal analyze --ligands`. |
| `curated_conf1_no_leakage.csv` | 2 KB | RefAHL conf=1 receptors retained after sequence-identity leakage filter. |
| `curated_conf2_no_leakage.csv` | 2 KB | Conf=2 equivalent. |
| `curated_conf3_no_leakage.csv` | 1 KB | Conf=3 equivalent. |

## Extracting the MSA archive

Before running `copal predict` against the benchmark, unpack the
bundled `.tar.xz`:

```bash
mkdir -p data/benchmark_msas
tar -xJf data/benchmark_msas.tar.xz -C data/benchmark_msas/
```

The extracted directory is `.gitignore`d. `msa_path` values in
`benchmark_data_all_corrected.csv` are repo-root-relative
(`data/benchmark_msas/<N>.a3m`), so run `copal` commands from the repo
root, or rewrite the column to absolute paths if you relocate the
MSAs: `df['msa_path'] = '/abs/path/' + df['msa_path'].str.replace('data/benchmark_msas/', '')`.

## Larger artifacts hosted externally

Files that are too large to commit live on Zenodo (DOIs filled in at
publication time):

| Artifact | Hosted | Notes |
| -------- | ------ | ----- |
| `predictions.csv` (~27 MB) | Zenodo: `<DOI_ANALYSIS>` | Pre-collected per-prediction ipTM/PAE table. Drop into this directory to run `copal analyze --input-csv data/predictions.csv`. |
| Raw model predictions (CIF + JSON, ~GBs) | Zenodo: `<DOI_PREDICTIONS>` | Optional: only needed if you want to re-run the analysis collection phase from scratch. |

## Optional reference assets

`reference/4y17.cif` (added on demand) is the reference *E. coli* SdiA /
3OC8 crystal structure used by `copal.utils.pocket` for pocket-residue
identification. Download from RCSB <https://www.rcsb.org/structure/4Y17>
and place at `data/reference/4y17.cif`, or set `COPAL_POCKET_REF_CIF` to
point at your local copy.
