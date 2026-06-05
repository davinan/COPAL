# The COPAL Ranking Procedure

## Notation

| Symbol | Meaning |
|---|---|
| $p$ | a LuxR receptor (protein) |
| $\mathcal{A}_p = \{a_1,\dots,a_N\}$ | its panel of $N$ candidate AHL ligands ($N=58$) |
| $\mathcal{M}$ | ensemble of $M$ structure models, $\{\text{AF3},\text{Boltz2},\text{Boltz2Affinity},\text{RF3}\}$, $M=4$ |
| $\mathcal{Q}_m(p,a)$ | predicted poses emitted by model $m$ for pair $(p,a)$ |
| $\mathrm{score}_m(\cdot)$ | raw confidence: ipTM (AF3/Boltz2/RF3), affinity\_probability\_binary (Boltz2Affinity) |
| $s_m(p,a)$ | pose-pooled per-model score |
| $r_m(p,a)$ | per-model rank of ligand $a$ for receptor $p$ |
| $R_{\text{COPAL}}(p,a)$ | final COPAL ensemble rank |

## Procedure

COPAL ranks the $N$ candidate ligands of each receptor $p$ independently, in three stages. First, multiple predictions per receptor–ligand pair are pooled by their maximum confidence:

$$s_m(p, a) = \max_{q \in \mathcal{Q}_m(p, a)} \mathrm{score}_m(q).$$

Within each receptor, every model ranks ligands by descending score, with tied scores assigned the most generous (minimum) rank:

$$r_m(p, a) = 1 + \bigl|\{\, a' \in \mathcal{A}_p : s_m(p, a') > s_m(p, a) \,\}\bigr|.$$

The per-model ranks are fused by a parameter-free rank-sum, and the COPAL rank is ascending in this aggregate (lower is better):

$$R_{\text{COPAL}}(p, a) = 1 + \Bigl|\Bigl\{\, a' \in \mathcal{A}_p : \sum_{m\in\mathcal{M}} r_m(p, a') < \sum_{m\in\mathcal{M}} r_m(p, a) \,\Bigr\}\Bigr|.$$

## Alternative fusions

As an alternative to rank-sum, COPAL supports a general fusion $S(p,a)$, with the final rank ascending in $S$ for the rank-reduction variants and descending for Reciprocal Rank Fusion:

$$S(p, a) = \begin{cases} \displaystyle\sum_{m \in \mathcal{M}} \frac{1}{k + r_m(p, a)}, & \text{RRF}\ (k = 60), \\[2ex] \displaystyle\phi_{m \in \mathcal{M}}\bigl(r_m(p, a)\bigr), & \phi \in \{\mathrm{sum},\mathrm{mean},\mathrm{median},\mathrm{min},\mathrm{max}\}, \end{cases}$$

where the paper default is $\phi = \mathrm{sum}$, recovering the COPAL rank above.
