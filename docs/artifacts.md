# Analysis artifacts

`copal analyze` writes seven artifact groups under
`<paths.output_root>/<set.name>/` (default
`outputs_official_workflow/results/<set.name>/`). `copal full-ranking` adds
the stylized `.xlsx` workbooks to the `full_rankings/` group.

Throughout, **rank 1 = best** (the model's most-confident prediction), and
"true AHL" is the experimentally validated cognate ligand of a LuxR.

| Group | Files | What it answers |
| ----- | ----- | --------------- |
| `ranking_performance/` | `ranking_performance.csv` | Where does each model rank the true AHL? |
| `mpkt/` | `mpkt.csv`, `mpkt_correlation.csv`, `mpkt_quartile_top6.csv` | Do the models agree, and does agreement predict success? |
| `acyl_moiety/` | `acyl_moiety.csv` + 4 per-model | True-AHL rank stratified by acyl-chain class. |
| `chemical_similarity/` | `chemical_similarity.csv` | How chemically close are the top picks to the true AHL? |
| `score_vs_rank/` | 4 per-model + `correlation.csv` | Does a model's own confidence predict its success? |
| `full_rankings/` | `full_rankings.csv` + 5 `.xlsx` | The complete rank matrix; stylized workbooks. |
| `ablation/` | `ablation_table.csv`, `ablation_table_raw.csv` | Which model combination ensembles best? |

---

## `ranking_performance/ranking_performance.csv`

One row per LuxR. Columns: `LuxR_name`, then one per model — `AF3`,
`Boltz2`, `Boltz2Affinity`, `RF3`, `Ensemble`. Each cell is the **rank of
the true AHL** in that model's ranking for that LuxR (1 = the model put the
cognate ligand first). Blank when the LuxR has no curated `true_ahl_id`.

## `mpkt/`

**`mpkt.csv`** — one row per LuxR. Columns:
- `LuxR_name`
- `mpKT` — *mean pairwise Kendall-τ*: the average rank correlation across
  all pairs of the four per-model rankings. High mpKT = the models agree
  with each other on how to rank the ligands (a label-free signal).
- `Ensemble True AHL Ranking` — rank of the true AHL in the ensemble.

**`mpkt_correlation.csv`** — three rows (`spearman`, `pearson`,
`kendall_tau`); columns `test`, `statistic`, `pvalue`, `n`. Correlation
between `mpKT` and `Ensemble True AHL Ranking` — i.e. does cross-model
agreement predict a better ensemble ranking of the true AHL.

**`mpkt_quartile_top6.csv`** — LuxRs split into mpKT quartiles; columns
`mpkt_quartile`, `quartile_min/max/center`, `n`,
`top1/top3/top5/top6_success_rate` (fraction of LuxRs in that quartile
whose ensemble ranked the true AHL within top-K).

## `acyl_moiety/`

**`acyl_moiety.csv`** — one row per LuxR. Columns `LuxR_name`, `All`,
`Short`, `Medium`, `Long`, `Atypical`. `All` is the ensemble's true-AHL
rank; the four category columns repeat that value only in the column
matching the true AHL's acyl-chain subcategory (else blank) — so a column
aggregates ranks for LuxRs whose cognate AHL is in that class.

**`acyl_moiety_AF3.csv`**, **`_Boltz2.csv`**, **`_Boltz2Affinity.csv`**,
**`_RF3.csv`** — the same table computed from a single model's ranking
instead of the ensemble.

## `chemical_similarity/chemical_similarity.csv`

One row per LuxR. Columns `LuxR_name`, `Top 1`, `Top 3`, `Top 5`, `Top 6`.
Each `Top K` value is the **maximum MCS similarity** between the true AHL
and any of the ensemble's top-K predicted ligands — a near-miss credit: if
the model's top pick is the wrong ligand but chemically very close to the
right one, the score stays high.

Similarity is a maximum-common-substructure (MCS) Jaccard index over heavy
atoms. The method is chosen by `analysis.chem_sim_matrix`: if it points at
a precomputed matrix (see below), values are looked up; if null, they are
computed on the fly with `mcs_order`.

## `score_vs_rank/`

Per-model tables asking whether a model's *own confidence* tracks its
success — using a **label-free** predictor computable at inference time.

**`AF3.csv`**, **`Boltz2.csv`**, **`RF3.csv`** — columns `LuxR_name`,
`mean_top6_iptm`, `True AHL Ranking`. `mean_top6_iptm` is the mean ipTM
over the model's six top-ranked ligand predictions for that LuxR.

**`Boltz2Affinity.csv`** — same, but the score column is
`mean_top6_affinity_probability_binary` (Boltz2's affinity-head output, not
an ipTM).

**`correlation.csv`** — one row per model; columns `model`, `score_col`,
`n`, `spearman_r`, `spearman_p`, `pearson_r`, `pearson_p`, `kendall_tau`,
`kendall_p`. Correlation between the `mean_top6_*` score and `True AHL
Ranking`. A negative correlation means higher confidence ↔ better (lower-
numbered) true-AHL rank.

## `full_rankings/`

**`full_rankings.csv`** — the wide rank matrix. One row per (LuxR, model)
pair (model ∈ {AF3, Boltz2, Boltz2Affinity, RF3, Ensemble}); meta columns
`LuxR_name`, `split`, `confidence`, `model_name`, `true_ahl_id`,
`true_ahl_abbrev`, `true_ahl_subcategory`, followed by one column per AHL
ligand holding that ligand's rank.

**`AF3.xlsx`**, **`Boltz2.xlsx`**, **`Boltz2Affinity.xlsx`**, **`RF3.xlsx`**,
**`ENSEMBLE.xlsx`** — stylized per-model workbooks, written by
`copal full-ranking` (not by `copal analyze`). Three rows per LuxR (AHL ids
/ abbreviations / scores). Highlighting: the true-AHL rank cell is filled
**yellow**; the `mpKT` column is filled **orange** when mpKT < 0.3 and
**green** when ≥ 0.3. In `set=inference` runs there is no ground truth, so
the yellow fill is absent and the meta columns are blank.

## `ablation/`

**`ablation_table_raw.csv`** — one row per non-empty subset of the four
models; the receptor pool's top-K success when only that subset is
ensembled. Shows which models carry the ensemble.

**`ablation_table.csv`** — the raw table summarized by number of models
included (mean / std of top-5 and top-6 success rate, count).

---

## `data/chem_sim_matrices/<method>.csv`

Precomputed pairwise ligand chemical-similarity matrices, produced by
`copal sim-matrix sim_matrix.method=<method>`. Each file is an N×N
symmetric matrix (N = ligand pool size, 58 for the paper benchmark),
indexed by `AHLid` on both axes, with the diagonal equal to 1.0.

`copal analyze … analysis.chem_sim_matrix=data/chem_sim_matrices/<method>.csv`
makes `chemical_similarity.csv` an O(1) lookup into the matrix instead of an
on-the-fly MCS computation.

Three methods are available (defined in `src/copal/analysis/chem_sim.py`),
all MCS-Jaccard over heavy atoms with element-aware atom matching and
`ringMatchesRingOnly=True`:

| Method | Bond handling | Behaviour on the AHL pool |
| ------ | ------------- | ------------------------- |
| `mcs_order` | `BondCompare.CompareOrder` — bond orders must match | Over-penalizes single bond-order differences: a saturated vs mono-unsaturated chain fragments the MCS and scores ~0.5. |
| `mcs_any` | `BondCompare.CompareAny` — bond order ignored | Under-penalizes: collapses saturated vs unsaturated (and 3-oxo vs 3-hydroxy) pairs to 1.0. |
| `mcs_penalized` | `CompareAny` match, then −1 atom-equivalent per bond-order mismatch in the matched substructure | Smooth middle ground; recommended for the AHL pool. Used for the paper benchmark. |

See [`reproduce_paper.md`](reproduce_paper.md) for the canonical
similarity-matrix → analysis workflow.
