# COPAL

**Co-folding ranking of protein–ligand pairs across AlphaFold 3, Boltz2, and RosettaFold3.**

COPAL takes a CSV of protein receptors and a CSV of small-molecule ligands,
runs each receptor × ligand pair through a configurable subset of
{AlphaFold 3, Boltz2, RosettaFold3} on a SLURM cluster, and ranks the
ligands per receptor using an ensemble of confidence-derived scores.

This repository accompanies the paper *“[paper title TBA]”* (Nakajima An et al.).

```
                  receptors.csv         ligands.csv
                       │                     │
                       ▼                     ▼
                       └──────► copal predict ◄──────┐
                                  │                  │
              ┌───────────────────┼───────────────────┐
              ▼                   ▼                   ▼
            AF3 JSON          Boltz2 YAML         RF3 JSON
              │                   │                   │
              └─── SLURM array job (per chunk) ──────┘
                                  │
                                  ▼
                       inputs/, outs/ (predictions)
                                  │
                                  ▼
                            copal analyze
                                  │
                       ┌──────────┴──────────┐
                       ▼                     ▼
              prediction-only        + benchmark mode
              (rankings.csv …)       (top-K, Kendall τ,
                                      ablation, etc.)
```

## Quickstart

Install the package with [`uv`](https://docs.astral.sh/uv/):

```bash
uv sync
```

Run the toy example (2 receptors × 3 ligands):

```bash
cd examples/quickstart
./run.sh                # writes per-model inputs and a SLURM array script
```

To predict on **your own** data, just supply two CSVs:

```bash
# Cluster (SLURM) — default
uv run copal predict \
    --receptors my_proteins.csv \
    --ligands   my_ligands.csv \
    --out-dir   ./run \
    --models    af3,boltz2,rf3

# Workstation (no SLURM)
uv run copal predict \
    --receptors my_proteins.csv \
    --ligands   my_ligands.csv \
    --out-dir   ./run \
    --models    af3,boltz2,rf3 \
    --executor  local
```

CSV schemas are documented in [`docs/input_formats.md`](docs/input_formats.md).
The minimum is `receptor_id,seq` and `ligand_id,smiles`.

After the SLURM array completes:

```bash
uv run copal analyze \
    --workflow-dir ./run/workflow \
    --output-root  ./run/results \
    --labels       data/benchmark_data_all_corrected.csv \
    --ligands      data/Final_confirmed_ligands_58_annotated.csv \
    --basepath     ./run/store/all
```

(Run `copal organize --pred-dir ./run/outs --dst-dir ./run/store` first
to build the hierarchical store that `analyze --basepath` expects.)

`--labels` is optional — without it, `analyze` produces ranked predictions
per receptor; with it, it additionally computes top-K accuracy, Kendall τ,
ablation tables, and chemical-similarity analyses.

## What you need installed

COPAL itself is pure Python (managed via `uv sync`). The structure
prediction models are external — the repository ships **documentation +
patches**, not the model source:

- **AlphaFold 3** — clone <https://github.com/google-deepmind/alphafold3>,
  apply the optional `third_party/alphafold3-patches/` patch (see
  [`INSTALLATION.md`](INSTALLATION.md)). CC-BY-NC-SA 4.0; non-commercial only.
- **Boltz2** — `pip install boltz` or clone <https://github.com/jwohlwend/boltz>.
  MIT licensed.
- **RosettaFold3 / Foundry** — clone the `production` branch of
  <https://github.com/RosettaCommons/foundry>. BSD-3-Clause.
- **BLAST+** — only required for the RALPH BLAST baseline and the
  `copal seqid-filter` leakage check.

Full setup (containers, weights, env vars) is in
[`INSTALLATION.md`](INSTALLATION.md).

## Reproducing the paper

Two reproduction paths in [`REPRODUCE.md`](REPRODUCE.md):

1. **Analysis-only** (no GPUs needed) — download a ~27 MB pre-collected
   ipTM CSV from a public data archive, run `copal analyze`, generate all
   paper tables.
2. **Full inference** — re-run AF3/Boltz2/RF3 on the benchmark, then
   `copal analyze`.

| Artifact | Hosted |
| -------- | ------ |
| `predictions.csv` | `<DOI_ANALYSIS>` |
| Raw model predictions (CIF + JSON) | `<DOI_PREDICTIONS>` |
| MSAs used at inference | `<DOI_MSAS>` |

Data host and DOIs will be filled in at publication time (TBD).

## Repository layout

```
copal/
├── src/copal/             # importable Python package
│   ├── pipeline.py        # schedule + emit AF3/Boltz2/RF3 inputs and SLURM script
│   ├── runner/            # in-tree AF3/Boltz2/RF3 managers + shell templates
│   ├── analysis/          # collect / postprocess / generate / leakage filter
│   ├── io/                # receptor & ligand CSV schema + validation
│   ├── prediction_path_managers.py
│   ├── organizer.py
│   └── cli.py             # `copal` entry point (Typer)
├── data/                  # small inputs shipped with the repo
├── examples/              # quickstart + benchmark pointers
├── baselines/ralph_blast/ # BLAST homology baseline
├── third_party/alphafold3-patches/  # 7-decimal AF3 confidence patch (CC-BY-NC-SA 4.0)
├── docs/                  # input formats, pipeline overview, data dictionary
├── pyproject.toml
└── README.md (this file)
```

## License

The COPAL source code is licensed under the **MIT License** (see
[`LICENSE`](LICENSE)). The AlphaFold 3 derivative file under
`third_party/alphafold3-patches/confidence_types.py` is **CC-BY-NC-SA 4.0**
(see that directory's `LICENSE` and `NOTICE`). External models invoked at
runtime carry their own licenses — see [`NOTICE`](NOTICE).

## Citation

```bibtex
@article{copal2026,
  title   = {<TBA>},
  author  = {Nakajima An, Davi and ...},
  journal = {<TBA>},
  year    = {2026},
  doi     = {<TBA>}
}
```
