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

### Labels CSV (`data/paper_luxrs.csv`)

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

### `curated_ligands_annotated.csv`

Per-AHL annotation table (58 entries). Consumed by `copal analyze` as the
set's `ligands_path`, and also valid as a `copal predict` ligands CSV.

| Column | Notes |
| ------ | ----- |
| `AHLid` | E.g. `AHL01`. |
| `abbreviation` | E.g. `3OC4`. |
| `Formula`, `Exact mas` | Pass-through. |
| `smi` | SMILES (`SMILES_COL` constant in analysis; aliased to `smiles` by `copal predict`). |
| `subcategory` | Acyl-chain length category (`Short`, `Medium`, `Long`, `Atypical`). |

### `curated_luxR_RefAHL_no_leakage.csv`

RefAHL receptors (all confidence tiers) that pass the leakage filter
(<40% identity to any reference PDB structure with the same cognate
ligand). Referenced by the `no_leakage_csvs` filter of the `*_no_leakage`
set configs; per-confidence narrowing is done by the `confidences` filter.

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

The pre-collected raw prediction table — one row per (receptor, ligand,
model, batch) prediction. Produced by the **collect** step (the deprecated
`python -m copal.analysis.pipeline` entrypoint; see [`../REPRODUCE.md`](../REPRODUCE.md)
Path B) or downloaded from the data archive. The most important columns:

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

One row per (receptor, ligand, model) — the postprocess output and the
source of truth for all `copal analyze` artifacts. Columns: `rank`, `luxr`,
`luxr_ahl`, `ahl`, `model`, `iptm`, `affinity_probability_binary` (the last
populated only for `boltz2affinity` rows; `iptm` is blank for those).

### `ensemble_rankings.csv`

One row per (receptor, ligand). Columns: `luxr`, `ahl`, `luxr_ahl`,
`ensemble_rank`, `ensemble_score`. The aggregation method (default `sum`;
`rrf` etc. optional) is set by `analysis.aggregation`.

### `ablation_table.csv`

Workflow-level ablation table written by postprocess: one row per
model-combination, with the receptor pool's top-K success under each
ensemble subset. The per-set, figure-ready version is re-emitted under
`results/<set.name>/ablation/` (see [`artifacts.md`](artifacts.md)).

### `done_pairs.json`

`{ "RECEPTOR_LIGAND": ["af3", "boltz2", ...], ... }` — pairs marked
already complete; can be passed to `--done-lookup` on subsequent
`copal predict` runs to skip re-checking the filesystem.

## Analysis outputs (`results/<set.name>/`)

`copal analyze` emits seven artifact groups; `copal full-ranking` adds the
stylized `.xlsx` workbooks. Each group's files, columns, and interpretation
are documented in **[`artifacts.md`](artifacts.md)** — the authoritative
reference for analysis outputs.

