# Input CSV formats

COPAL accepts two simple CSV files. Both are validated by
`copal.io.schema` on load (helpful errors on missing or duplicate IDs).

## `receptors.csv`

| Column | Required | Type | Description |
| ------ | -------- | ---- | ----------- |
| `receptor_id` | ✅ | string | Unique identifier (used in filenames). Stripped of whitespace. Must be unique across rows. |
| `seq` | ✅ | string | Single-chain amino-acid sequence (one-letter codes). Empty rows are dropped on load. |
| `msa_path` | optional | string | Path to a pre-built MSA in `.a3m` format. If omitted, the model runner will build one at inference time (requires MSA tooling per `INSTALLATION.md`). |
| `chain_count` | optional | integer | Number of identical chains for homo-oligomers. Default 1. |
| `template_paths` | optional | string | Semicolon-separated list of template CIF paths (AF3 only; pass-through). |

### Backward-compatible aliases

If your CSV uses one of these column names, COPAL renames it on load:

| Found name | Renamed to |
| ---------- | ---------- |
| `LuxR_name`, `refahl_id`, `target` | `receptor_id` |
| `sequence` | `seq` |

### Example

```csv
receptor_id,seq,chain_count
LasR,MALVDGFLELERSSGKL...,2
SdiA,MQDKDFFSWRRTMLLR...,1
```

## `ligands.csv`

| Column | Required | Type | Description |
| ------ | -------- | ---- | ----------- |
| `ligand_id` | ✅ | string | Unique identifier (used in filenames). |
| `smiles` | ✅ | string | SMILES string. |
| `inchi`, `name`, `formula`, `mol_weight`, … | optional | any | Pass-through metadata; ignored by inference. |

### Backward-compatible aliases

| Found name | Renamed to |
| ---------- | ---------- |
| `AHLid`, `ahl_id`, `abbreviation` | `ligand_id` |
| `Isomeric SMILES`, `smi`, `SMILES` | `smiles` |

### Example

```csv
ligand_id,smiles,name
3OC8,CCCCCC(=O)CC(=O)N[C@H]1CCOC1=O,3-oxo-C8 HSL
C4,CCCC(=O)N[C@H]1CCOC1=O,butyryl HSL
```

## Labels CSV (used only by `copal analyze` in benchmark mode)

This is a *receptors-like* CSV with extra annotation columns that enable
performance metrics:

| Column | Used for |
| ------ | -------- |
| `refahl_id` (or `receptor_id`) | Match against rankings. |
| `seq` | Sanity-check against `receptors.csv`. |
| `true_ahl_abbrev` | Ground-truth ligand abbreviation. |
| `true_ahl_id` | Ground-truth ligand ID. |
| `split` | Cohort label (`RefAHL`, `wetlabexp`, `solo`, …). |
| `refahl_confidence_rating` | 1 / 2 / 3 confidence tier of the label. |

The benchmark file shipped at `data/paper_luxrs.csv` illustrates the full
schema.

## Output schemas

For the column meanings of every CSV produced by `copal analyze` and
`copal full-ranking`, see [`artifacts.md`](artifacts.md). For input and
intermediate file schemas, see [`data_dictionary.md`](data_dictionary.md).
