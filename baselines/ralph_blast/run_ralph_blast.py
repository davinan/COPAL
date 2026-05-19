#!/usr/bin/env python3
"""
RALPH BLAST predictor: run BLASTP per LuxR, aggregate rankings, write
full_rankings and ranking_performance using the COPAL analysis machinery.

Setup
-----
The `ralph` Python package is required (see baselines/ralph_blast/README.md
for installation and database build instructions).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import numpy as np

# RALPH package — see baselines/ralph_blast/README.md for install instructions.
from ralph.blastp import blastp_query
from ralph.ranking import aggregate_rankings
from ralph.database import load_scores


def _normalize_luxr(refahl_id: str) -> str:
    """Same convention as analysis: stripped uppercase."""
    return str(refahl_id).strip().upper()


def run_blast_and_aggregate(
    seq: str,
    scores_df: pd.DataFrame,
    db_path: str,
    max_hits: int,
    min_pident: float,
    max_pident: float,
) -> pd.DataFrame | None:
    """
    Run BLASTP, filter hits by pident in [min_pident, max_pident), get score rows, aggregate_rankings.
    Returns aggregated DataFrame (name, ranking, ...) or None if no hits in range.
    """
    hits = blastp_query(query_seq=seq, db_path=db_path, max_hits=max_hits)
    if hits.empty:
        return None
    hits = hits.sort_values("bitscore", ascending=False)
    print(hits['pident'].describe())
    hits = hits[(hits["pident"] >= min_pident) & (hits["pident"] < max_pident)]
    if hits.empty:
        return None
    hit_ids = hits["sseqid"].tolist()
    score_rows = scores_df[scores_df["luxr_name"].isin(hit_ids)]
    if score_rows.empty:
        return None
    return aggregate_rankings(score_rows)

EXCLUDE_SUFFIXES = ("LONG", "SHORT", "MED")

def apply_global_suffix_exclusion(labels_df: pd.DataFrame) -> pd.DataFrame:
    """Drop label rows whose LuxR_name ends with EXCLUDE_SUFFIXES."""
    return labels_df[
        ~labels_df["LuxR_name"].str.endswith(EXCLUDE_SUFFIXES, na=False)
    ].copy()

def main() -> None:
    parser = argparse.ArgumentParser(description="RALPH BLAST predictor: BLASTP + aggregate -> full_rankings & ranking_performance")
    parser.add_argument("--benchmark-csv", type=Path, required=True, help="Path to paper_luxrs.csv (labels)")
    parser.add_argument("--ralph-database", type=Path, required=True, help="Path to RALPH database dir (contains ralph_scores.csv, ralph_db)")
    parser.add_argument("--ligands-path", type=Path, required=True, help="Path to ligands CSV (AHLid, abbreviation, subcategory)")
    parser.add_argument("--output-root", type=Path, required=True, help="Root for full_rankings/ and ranking_performance/")
    parser.add_argument("--workflow-dir", type=Path, default=None, help="Dir for rankings/ensemble CSVs (default: output_root/../workflow_ralph)")
    parser.add_argument("--min-pident", type=float, default=30.0, help="Min percent identity (default: 30)")
    parser.add_argument("--max-pident", type=float, default=100.0, help="Max percent identity, exclusive (default: 100)")
    parser.add_argument("--max-hits", type=int, default=50, help="BLASTP max_target_seqs (default: 50)")
    args = parser.parse_args()

    if args.workflow_dir is None:
        args.workflow_dir = args.output_root.parent / "workflow_ralph"
    args.workflow_dir = args.workflow_dir.resolve()
    args.output_root = args.output_root.resolve()
    db_path = str(args.ralph_database.resolve() / "ralph")
    if not (args.ralph_database / "ralph_scores.csv").exists():
        raise FileNotFoundError(f"RALPH scores not found at {args.ralph_database / 'ralph_scores.csv'}")

    # Load benchmark and RALPH scores
    benchmark = pd.read_csv(args.benchmark_csv)
    benchmark = benchmark[benchmark["seq"].notna() & (benchmark["seq"].astype(str).str.strip() != "")]
    scores_df = load_scores(database_path=str(args.ralph_database.resolve()))

    ranking_rows = []
    ensemble_rows = []

    for _, row in benchmark.iterrows():
        refahl_id = row["refahl_id"]
        seq = str(row["seq"]).strip()
        luxr = _normalize_luxr(refahl_id)
        agg = run_blast_and_aggregate(
            seq=seq,
            scores_df=scores_df,
            db_path=db_path,
            max_hits=args.max_hits,
            min_pident=args.min_pident,
            max_pident=args.max_pident,
        )
        if agg is None or agg.empty:
            continue
        for _, arow in agg.iterrows():
            ahl = str(arow["name"]).strip().upper()
            rank = int(arow["ranking"])
            ranking_rows.append({"luxr": luxr, "model": "ralph", "ahl": ahl, "rank": rank, "iptm": 1.0 / max(rank, 1)})
            ensemble_rows.append({
                "luxr_ahl": f"{luxr}_{ahl}",
                "ensemble_rank": rank,
                "ensemble_score": 1.0 / max(rank, 1),
            })

    rankings_df = pd.DataFrame(ranking_rows)
    ensemble_df = pd.DataFrame(ensemble_rows)
    if rankings_df.empty:
        raise RuntimeError("No rankings produced (no hits in identity range for any LuxR).")

    args.workflow_dir.mkdir(parents=True, exist_ok=True)
    rankings_path = args.workflow_dir / "per_model_rankings.csv"
    ensemble_path = args.workflow_dir / "ensemble_rankings.csv"
    rankings_df.to_csv(rankings_path, index=False)
    ensemble_df.to_csv(ensemble_path, index=False)
    print(f"Wrote {rankings_path} ({len(rankings_df)} rows)")
    print(f"Wrote {ensemble_path} ({len(ensemble_df)} rows)")

    # ----- Generate phase: reuse the COPAL analysis module, monkey-patch RALPH, write outputs -----
    from copal.analysis import pipeline as analysis

    analysis.MODEL_NAME_MAP["ralph"] = "RALPH"
    analysis.RANKING_MODELS_ORDER = list(analysis.RANKING_MODELS_ORDER) + ["RALPH"]

    labels_raw = pd.read_csv(args.benchmark_csv)
    ligands_raw = pd.read_csv(args.ligands_path)
    rankings_raw = pd.read_csv(rankings_path)
    ensemble_raw = pd.read_csv(ensemble_path)

    labels_norm, rankings_norm, ensemble_norm, ligands_norm = analysis.normalize_inputs(
        labels_raw, rankings_raw, ensemble_raw, ligands_raw
    )
    labels_norm = apply_global_suffix_exclusion(labels_norm)

    full_dir = args.output_root / "full_rankings"
    ranking_dir = args.output_root / "ranking_performance"
    full_dir.mkdir(parents=True, exist_ok=True)
    ranking_dir.mkdir(parents=True, exist_ok=True)

    # Ranking performance: RefAHL conf 1, 2, 3
    for conf in [1, 2, 3]:
        subset = labels_norm[
            (labels_norm["split"].str.upper() == "REFAHL") & (labels_norm["confidence"] == conf)
        ].copy()
        subset = subset.drop_duplicates(subset=["LuxR_name", "split", "confidence"])
        true_rank_lookup = analysis.build_true_rank_lookup(subset, rankings_norm, ensemble_norm)
        perf_df = analysis.build_ranking_performance(subset, true_rank_lookup)
        perf_path = ranking_dir / f"ranking_performance_conf{conf}.csv"
        analysis.write_results_csv(perf_df, perf_path)
        print(f"Wrote {perf_path}")

    # Full rankings: all_splits.csv
    full_rankings_df = analysis.build_full_rankings(
        labels_norm, rankings_norm, ensemble_norm, ligands_norm
    )
    full_path = full_dir / "all_splits.csv"
    analysis.write_results_csv(full_rankings_df, full_path)
    print(f"Wrote {full_path}")

    # mpKT placeholder for Excel (one row per LuxR in rankings)
    luxrs_in_rankings = rankings_norm["luxr"].unique()
    mpkt_df = pd.DataFrame({"LuxR_name": luxrs_in_rankings, "mpKT": np.nan})

    excel_path = full_dir / "RALPH.xlsx"
    analysis.build_ranked_excel(
        "ralph", rankings_norm, ensemble_norm, labels_norm, ligands_norm, mpkt_df, excel_path
    )
    print(f"Wrote {excel_path}")


if __name__ == "__main__":
    main()
