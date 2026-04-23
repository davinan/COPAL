"""Ranking aggregation functions for consensus ligand prediction."""

import pandas as pd
import numpy as np
from typing import Optional


def aggregate_rankings(
    df_top: pd.DataFrame,
    weight_col: Optional[str] = None,
    weight_col_min: float = 0.01,
    rrf_k: float = 0.0,
    power: float = 1.0,
) -> pd.DataFrame:
    """Aggregates ranking columns across df_top into a final consensus ranking.
    
    This function performs Reciprocal Rank Fusion (RRF) style aggregation
    across multiple ranking columns (e.g., ensemble_1, ensemble_2, ...).
    
    Args:
        df_top: DataFrame containing ranking columns (columns starting with "ensemble_")
        weight_col: Optional column name to weight rankings by (e.g., "KT")
        weight_col_min: Minimum weight value (clips weights below this)
        rrf_k: RRF k parameter (bigger -> flatter; smaller -> top ranks dominate)
        power: Power parameter for position decay
        
    Returns:
        DataFrame with columns: name, score, count, avg_position, ranking
        Sorted by score (descending)
    """
    ranking_cols = [c for c in df_top.columns if c.startswith("ensemble_") and not c.endswith("std")]
    
    # Sort by numeric suffix if present
    def _suffix_num(c):
        try:
            return int(c.split("_")[-1])
        except Exception:
            return 10**9
    
    ranking_cols = sorted(ranking_cols, key=_suffix_num)

    if len(ranking_cols) == 0:
        raise ValueError("No ranking columns found (expected columns like ensemble_1, ensemble_2, ...).")

    # Position decay (how much rank_1 vs rank_N matters)
    npos = len(ranking_cols)
    positions = np.arange(1, npos + 1)  # 1..npos

    # Reciprocal Rank Fusion style: 1 / (k + r)^power
    pos_weight = 1.0 / np.power(rrf_k + positions, power)

    # Accumulate votes
    scores = {}
    counts = {}
    pos_sums = {}
    
    if weight_col is not None:
        w = np.clip(df_top[weight_col].values, a_min=weight_col_min, a_max=None)
    else:
        w = np.ones(len(df_top))
    
    for row_i, (_, row) in enumerate(df_top.iterrows()):
        row_w = w[row_i]
        for j, col in enumerate(ranking_cols):
            name = row[col]
            if pd.isna(name):
                continue
            name = str(name).strip()
            if not name:
                continue

            vote = row_w * float(pos_weight[j])
            scores[name] = scores.get(name, 0.0) + vote
            counts[name] = counts.get(name, 0) + 1
            pos_sums[name] = pos_sums.get(name, 0.0) + (j + 1)  # 1-indexed position

    out = pd.DataFrame({
        "name": list(scores.keys()),
        "score": [scores[k] for k in scores.keys()],
        "count": [counts[k] for k in scores.keys()],
        "avg_position": [pos_sums[k] / counts[k] for k in scores.keys()],
    })

    out = out.sort_values(["score", "count"], ascending=[False, False]).reset_index(drop=True)
    out['ranking'] = np.arange(len(out)) + 1
    return out
