# Data dictionary

Schemas of every CSV that COPAL reads or writes. All file names assume the
default working-directory layout from
[`pipeline_overview.md`](pipeline_overview.md).

## Inputs

### `receptors.csv` (read by `copal predict`, `copal analyze`)

| Column | Type | Notes |
| ------ | ---- | ----- |
| `receptor_id` | str | Unique. |
| `seq` | str | Single-chain amino-acid sequence. |
| `msa_path` | str (optional) | Path to `.a3m`. |
| `chain_count` | int (optional) | Default 1. |

See [`input_formats.md`](input_formats.md) for the full spec.

### `ligands.csv` (read by `copal predict`, `copal analyze`)

| Column | Type | Notes |
| ------ | ---- | ----- |
| `ligand_id` | str | Unique. |
| `smiles` | str | SMILES. |
| `abbreviation`, `formula`, `subcategory` | optional | Pass-through. |

### Labels CSV (e.g. `data/benchmark_data_all_corrected.csv`)

Used only when running `copal analyze` in benchmark mode.

| Column | Notes |
| ------ | ----- |
| `refahl_id` | Receptor ID (mapped to `receptor_id`). |
| `seq` | Sequence. |
| `true_ahl_id` | Ground-truth AHL ID (e.g. `AHL04`). |
| `true_ahl_abbrev` | Ground-truth abbreviation (e.g. `3OC10`). |
| `split` | Cohort: `RefAHL`, `wetlabexp`, `solo`, … |
| `refahl_confidence_rating` | 1 / 2 / 3 — label confidence tier. |
| `crystal_pdb` | PDB ID, if any (informational). |
| `aa_length` | Sequence length. |
| `strain_or_metagenome`, `pubmed_ids`, `notes` | Free-text. |

### `Final_confirmed_ligands_58_annotated.csv`

Per-AHL annotation table consumed by the analysis generate phase.

| Column | Notes |
| ------ | ----- |
| `AHLid` | E.g. `AHL01`. |
| `abbreviation` | E.g. `3OC4`. |
| `Formula`, `Exact mas` | Pass-through. |
| `smi` | SMILES (`SMILES_COL` constant in analysis). |
| `subcategory` | Acyl-chain length category (`Short`, `Medium`, `Long`). |

### `curated_conf{1,2,3}_no_leakage.csv`

Per-confidence subset of RefAHL receptors that pass the leakage filter
(<40% identity to any reference PDB structure with the same cognate
ligand).

| Column | Notes |
| ------ | ----- |
| `refahl_id` | Surviving receptor ID. |
| `cognate_ligand` | Ground-truth ligand abbreviation. |
| `max_identity_to_any_pdb` | Highest sequence identity observed. |
| `closest_pdb`, `closest_pdb_ligand` | Closest reference structure. |
| `reason` | `<40% identity` or `≥40% but different ligand`. |

## Intermediate outputs

### `scheduled_ensemble.csv`

Built by `copal predict`; one row per (receptor, ligand) pair.

| Column | Notes |
| ------ | ----- |
| `target_name` | `<receptor_id>_<ligand_id>`. |
| `protein_name`, `ligand_name` | Source IDs. |
| `batch_id` | Integer batch index (controls SLURM array ID). |
| `is_done_<model>` | True once the model's prediction artifacts exist. |

### `predictions.csv`

Output of analysis Phase 1 (collect). One row per (receptor, ligand,
model, batch) prediction. The most important columns:

| Column | Notes |
| ------ | ----- |
| `protein_name` | Receptor ID. |
| `ligand_name` | Ligand ID. |
| `model` | `af3`, `boltz2`, `boltz2affinity`, or `rf3`. |
| `batch_idx` | Per-prediction sample index (0–4 typically). |
| `iptm` | Interface pTM. |
| `min_pae_interface`, `pae_interface` | PAE-derived interface scores. |
| `affinity_pred_value*`, `affinity_probability_binary*` | Boltz2 affinity head outputs. |

### `per_model_rankings.csv`

One row per (receptor, ligand, model). Columns: `luxr`, `model`, `ahl`,
`rank`, `iptm`.

### `ensemble_rankings.csv`

One row per (receptor, ligand). Columns: `luxr_ahl`, `ensemble_rank`,
`ensemble_score`.

### `ablation_table.csv`

One row per ablation cell (model subset × split × confidence). Columns
include the ablation key, top-K accuracy, Kendall τ, etc.

### `done_pairs.json`

`{ "RECEPTOR_LIGAND": ["af3", "boltz2", ...], ... }` — pairs marked
already complete; can be passed to `--done-lookup` on subsequent
`copal predict` runs to skip re-checking the filesystem.

## Analysis outputs (`results/`)

### `ranking_performance/`

`ranking_performance_conf<K>.csv`: top-K accuracy, Kendall τ, Spearman ρ
per model + ensemble.

### `full_rankings/`

`all_splits.csv` — per-receptor ranked ligand list.
`<MODEL>.xlsx` — Excel workbooks per model (one sheet per receptor with
predicted vs. true ranking).

### `ablation/`

Filtered, re-formatted ablation tables for figure consumption.

### `chemical_similarity/` (requires RDKit)

`chemical_similarity_conf<K>.csv` — MCS / Morgan-fingerprint similarity
between predicted top-1 and ground-truth ligand, per RefAHL confidence
stratum.

### `acyl_moiety/`

Acyl-chain-length match diagnostics.

### `mpkt/`

Per-receptor ranked tables formatted for the mpKT chemoinformatic
follow-up analyses.

