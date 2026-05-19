# Reproducing the paper benchmark

This is the canonical command sequence for the published COPAL benchmark.
It assumes you already have the postprocessed ranking CSVs
(`per_model_rankings.csv` + `ensemble_rankings.csv`) in
`outputs_official_workflow/` — either from the data archive or rebuilt via
[`../REPRODUCE.md`](../REPRODUCE.md) Path A/B.

All three commands here are **Hydra-driven**: arguments are `key=value`
config overrides, not `--flags`.

## Environment

Set `COPAL_REPO_ROOT` so the Hydra `paths` config resolves `data/` and
`outputs_official_workflow/` correctly when you run from anywhere other
than the repo root:

```bash
export COPAL_REPO_ROOT=/path/to/COPAL          # the repo root
```

## The benchmark set: `final_n96`

The published set is `set=final_n96` (`src/copal/config/set/final_n96.yaml`):
RefAHL split, all confidence tiers {1,2,3}, intersected with the union of
the three curated no-leakage CSVs, with replicate (`_R1`, `_R2`) and litter
(`LONG`, `SHORT`, `MED`) suffixes dropped — **and** `Vib_001c` excluded (a
known problematic entry held out of the benchmark). After all filters,
exactly **N=96** LuxRs remain.

## Step 1 — precompute the ligand-similarity matrix

`chemical_similarity` compares predicted ligands to the true ligand by MCS
similarity. Precompute the pairwise matrix once with the `mcs_penalized`
method (see [`artifacts.md`](artifacts.md) for why this method — it handles
bond-order differences smoothly, unlike `mcs_order` / `mcs_any`):

```bash
uv run copal sim-matrix sim_matrix.method=mcs_penalized
# → data/chem_sim_matrices/mcs_penalized.csv  (58×58 symmetric)
```

## Step 2 — stylized ranking workbooks

```bash
uv run copal full-ranking set=final_n96
# → outputs_official_workflow/results/final_n96/full_rankings/
#     {AF3,Boltz2,Boltz2Affinity,RF3,ENSEMBLE}.xlsx + full_rankings.csv
```

## Step 3 — the analysis artifacts

```bash
uv run copal analyze \
    set=final_n96 \
    analysis.chem_sim_matrix=data/chem_sim_matrices/mcs_penalized.csv
# → outputs_official_workflow/results/final_n96/
#     ranking_performance/  mpkt/  acyl_moiety/  chemical_similarity/
#     score_vs_rank/  full_rankings/  ablation/
```

Every artifact group is documented in [`artifacts.md`](artifacts.md).

## Aggregation: SUM (default) vs RRF

The ensemble combines the four per-model rankings. The default aggregation
is **`sum`** (sum of per-model ranks) — the commands above use it
implicitly, no override needed.

**Reciprocal Rank Fusion (RRF)** is the main alternative. To produce the
RRF variant *alongside* the default run (a distinct `set.name` keeps the
outputs separate):

```bash
uv run copal full-ranking set=final_n96 set.name=final_n96_rrf \
    analysis.aggregation=rrf

uv run copal analyze set=final_n96 set.name=final_n96_rrf \
    analysis.aggregation=rrf \
    analysis.chem_sim_matrix=data/chem_sim_matrices/mcs_penalized.csv
```

`analysis.aggregation` also accepts `mean`, `median`, `min`, `max`;
`analysis.rrf_k` (default 60) tunes the RRF constant. Aggregation affects
only ensemble-derived columns — the per-model artifacts (`score_vs_rank/`,
per-model `acyl_moiety_*`, the model columns of `ranking_performance`) are
identical regardless of the choice.

## Case studies

The paper's case-study figures use a hand-picked LuxR subset,
`set=paper_case_studies` (an `include_names` allowlist; see
`src/copal/config/set/paper_case_studies.yaml`). Run it the same way:

```bash
uv run copal full-ranking set=paper_case_studies
uv run copal analyze set=paper_case_studies \
    analysis.chem_sim_matrix=data/chem_sim_matrices/mcs_penalized.csv
```
