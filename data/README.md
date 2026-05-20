# `data/`

Small input files used by the analysis pipeline and shipped with the
repository for full reproducibility of the paper benchmark.

| File | Size | Purpose |
| ---- | ---- | ------- |
| `paper_luxrs.csv` | ~115 KB | Per-receptor labels: `refahl_id`, `seq`, `true_ahl_id`, `true_ahl_abbrev`, `split` (`RefAHL`, `wetlabexp`, `solo`, …), `refahl_confidence_rating`, `msa_path`. Referenced as `labels_path` by the Hydra `set` configs; also usable as a receptors CSV for `copal predict`. `msa_path` values are repo-root-relative (`data/benchmark_msas/*.a3m`); extract `benchmark_msas.tar.xz` before running inference (see below). |
| `benchmark_msas.tar.xz` | 31 MB | ColabFold-style `.a3m` MSAs for all 281 benchmark receptors (one file per `refahl_id`, keyed by the number in the `msa_path` column). Extract with `tar -xJf data/benchmark_msas.tar.xz -C data/benchmark_msas/`. |
| `curated_ligands_annotated.csv` | 5 KB | Annotated AHL ligand metadata (58 entries; `AHLid`, `abbreviation`, `Formula`, `smi`, `subcategory`). Referenced as `ligands_path` by the `set` configs (`copal analyze`), and also valid as a `copal predict` ligands CSV (the `smi`/`AHLid` columns are auto-aliased to `smiles`/`ligand_id`). |
| `curated_luxR_RefAHL_no_leakage.csv` | 5 KB | RefAHL receptors (all confidence tiers) retained after the sequence-identity leakage filter (`refahl_id`, `cognate_ligand`, `max_identity_to_any_pdb`, `closest_pdb`, `closest_pdb_ligand`, `reason`). Referenced by the `no_leakage_csvs` filter of the `*_no_leakage` set configs; per-confidence narrowing is done by the `confidences` filter axis. |

## Extracting the MSA archive

Before running `copal predict` against the benchmark, unpack the
bundled `.tar.xz`:

```bash
mkdir -p data/benchmark_msas
tar -xJf data/benchmark_msas.tar.xz -C data/benchmark_msas/
```

The extracted directory is `.gitignore`d. `msa_path` values in
`paper_luxrs.csv` are repo-root-relative (`data/benchmark_msas/<N>.a3m`),
so run `copal` commands from the repo root, or rewrite the column to
absolute paths if you relocate the MSAs:
`df['msa_path'] = '/abs/path/' + df['msa_path'].str.replace('data/benchmark_msas/', '')`.

## Larger artifacts hosted externally

Files that are too large to commit live on an external data host (host
and DOIs TBD, filled in at publication time):

| Artifact | Hosted | Notes |
| -------- | ------ | ----- |
| `predictions.csv` (~27 MB) | `<DOI_ANALYSIS>` (host TBD) | Pre-collected per-prediction ipTM/PAE/affinity table (`batch_idx, luxr, ahl, luxr_ahl, pae_interface, min_pae_interface, iptm, seed, model, affinity_*`). Download and place in this directory; not shipped in-repo (git-ignored). |
| Raw model predictions (CIF + JSON, ~GBs) | `<DOI_PREDICTIONS>` (host TBD) | Optional: only needed if you want to re-run the analysis collection phase from scratch. |

## Optional reference assets

`reference/4y17.cif` (added on demand) is the reference *E. coli* SdiA /
3OC8 crystal structure used by `copal.utils.pocket` for pocket-residue
identification. Download from RCSB <https://www.rcsb.org/structure/4Y17>
and place at `data/reference/4y17.cif`, or set `COPAL_POCKET_REF_CIF` to
point at your local copy.
