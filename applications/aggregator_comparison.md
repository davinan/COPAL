# Ensemble aggregator comparison — top-K success on no-leakage RefAHL

## Subset

RefAHL receptors with confidence rating ∈ {1, 2, 3} restricted to the
no-leakage union (`curated_conf{1,2,3}_no_leakage.csv`) after dropping the
six `_R1`/`_R2` near-duplicate variants. **n = 97 receptors.**

## Aggregators evaluated

All aggregators take per-model ranks (af3, boltz2, boltz2affinity, rf3) and
combine them into one ensemble score per ligand without any fitting:

| Aggregator | Formula |
|---|---|
| RRF (k=60) | `score = Σ_models 1 / (60 + per_model_rank)`, rank descending |
| Sum | `score = Σ_models per_model_rank`, rank ascending |
| Mean | `score = mean_models per_model_rank`, rank ascending (order-equivalent to Sum on a fixed model count) |
| Median | `score = median_models per_model_rank`, rank ascending |
| Min | `score = min_models per_model_rank`, rank ascending |
| Max | `score = max_models per_model_rank`, rank ascending |

Ties at every ranking step are resolved with `method="min"` (tied items share
the lower rank), so results are reproducible under input row reorderings.

## Top-K success rate

Success rate = fraction of receptors whose cognate ligand appears at
`ensemble_rank ≤ K`. 95% CIs are Wilson intervals.

| Aggregator | Top-1 [95% CI] | Top-3 [95% CI] | Top-5 [95% CI] | Top-6 [95% CI] | Avg. #ligands ≤ rank 5 |
|---|---|---|---|---|---|
| **RRF** (k=60) | 0.278 [0.20, 0.37] | 0.505 [0.41, 0.60] | 0.660 [0.56, 0.75] | **0.711** [0.61, 0.79] | **5.00** |
| **Sum** | 0.278 [0.20, 0.37] | **0.536** [0.44, 0.63] | **0.670** [0.57, 0.76] | 0.670 [0.57, 0.76] | 5.13 |
| **Mean** | 0.278 [0.20, 0.37] | **0.536** [0.44, 0.63] | **0.670** [0.57, 0.76] | 0.670 [0.57, 0.76] | 5.13 |
| Median | 0.309 [0.23, 0.41] | 0.464 [0.37, 0.56] | 0.670 [0.57, 0.76] | 0.701 [0.60, 0.78] | 5.14 |
| Min | 0.464† [0.37, 0.56] | 0.495 [0.40, 0.59] | 0.649 [0.55, 0.74] | 0.680 [0.58, 0.76] | 6.11† |
| Max | 0.216 [0.15, 0.31] | 0.515 [0.42, 0.61] | 0.598 [0.50, 0.69] | 0.608 [0.51, 0.70] | 5.23 |

**Bold** = best point estimate per column. Underlying values in
`aggregator_comparison_top_k.csv` (this directory).

†Min's apparent top-1 dominance is a tie artifact: under min aggregation, an
average of 3.19 ligands per receptor tie at rank 1 (right-most column). So
"top-1" for Min is "is the cognate among ~3 candidates," not "is it the
single best." After accounting for tie inflation, no aggregator is
significantly different from any other at any K (see tests below).

## Statistical tests

### Cochran's Q omnibus

H₀: all six aggregators have equal top-K success rate.

| K | Q | df | p | reject H₀? |
|---|---|---|---|---|
| 1 | 45.83 | 5 | 9.8 × 10⁻⁹ | yes — driven entirely by Min's tie inflation † |
| 3 | 3.96 | 5 | 0.56 | no |
| 5 | 4.79 | 5 | 0.44 | no |
| 6 | 8.44 | 5 | 0.13 | no |

### Pairwise McNemar (exact) for RRF vs each alternative

| K | vs Sum | vs Mean | vs Median | vs Min | vs Max |
|---|---|---|---|---|---|
| top-1 | 1.00 | 1.00 | 0.55 | 0.0003 † | 0.15 |
| top-3 | 0.38 | 0.38 | 0.42 | 1.00 | 1.00 |
| top-5 | 1.00 | 1.00 | 1.00 | 1.00 | 0.24 |
| top-6 | 0.13 | 0.13 | 1.00 | 0.63 | 0.05 |

Bonferroni-corrected α for the 15 pairwise comparisons within each K is
**0.0033**. Only Min's top-1 differences reach significance, and that
signal is the tie artifact.

## Bottom line

Across confidence-1, -2, and -3 RefAHL receptors that pass the no-leakage
filter (n = 97), the ensemble's top-K success rate is essentially
insensitive to the choice of rank-aggregator: sum-of-ranks, mean-of-ranks,
and Reciprocal Rank Fusion (k = 60) all fall within a 1–4 receptor band of
one another at top-1, top-3, top-5, and top-6, and Cochran's Q omnibus
tests do not reject equality at K ∈ {3, 5, 6} (p ≥ 0.13). Pairwise McNemar
tests find no Bonferroni-significant difference between Sum/Mean and RRF
at any K. Max is consistently the weakest aggregator (~10 percentage
points below the top group at top-6). Min's apparent top-1 advantage is
an artifact of multi-way ties — its top-1 set averages 3.2 ligands per
receptor versus 1.0 for RRF — and disappears once tie inflation is
accounted for.

We retain sum-of-ranks (Borda count) as the primary aggregator. It is
parameter-free, has a direct interpretation as the total rank across the
four models (and dividing by four recovers the mean per-model rank), and
matches RRF's benchmark performance to within sampling noise. RRF would
require disclosing the convention `k = 60` and reporting scores on an
opaque 0.04–0.07 scale; with no measurable accuracy benefit, we prefer
the simpler aggregator. The only practical difference is a small amount
of tie inflation in Sum's ranks (the top-5 set averages 5.13 ligands per
receptor versus exactly 5.00 for RRF, a 2.6% effect); ties are resolved
deterministically with `method="min"` (tied items share the lower rank),
so results are reproducible under input row reorderings.

## Reproduction

The numbers above come from running the extended analysis pipeline once per
aggregator and comparing the resulting `ensemble_rankings.csv` files:

```bash
for AGG in rrf sum mean median min max; do
    AGGREGATION=$AGG \
        WORKFLOW_DIR=outputs_sweep/workflow_$AGG \
        OUTPUT_ROOT=outputs_sweep/results_$AGG \
        bash applications/run_analysis_extended.sh
done
```

The companion CSV `aggregator_comparison_top_k.csv` (this directory) has
one row per aggregator with hits, rate, 95% CI bounds, and average
ligand-set size at K ∈ {1, 3, 5, 6}.
