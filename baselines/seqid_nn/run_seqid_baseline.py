#!/usr/bin/env python3
"""Sequence-identity "nearest-leak" baseline.

A deliberately simple, sequence-only baseline for the COPAL AHL-ranking benchmark.

Idea
----
The benchmark is evaluated on a *leak-free* receptor set (receptors NOT too similar
to any training PDB). The receptors that *were* excluded ("leaked" examples) come with
a known cognate AHL. We use those leaked receptors as a sequence -> AHL database:

  for each leak-free benchmark receptor (query):
      1. find the closest leaked receptor by global-alignment %% identity (BLOSUM62)
      2. take that leaked receptor's cognate AHL as the "anchor"
      3. score every candidate AHL by its MCS bond-penalized chemical similarity
         to the anchor (precomputed ``mcs_penalized`` matrix) and rank descending.

Outputs the same benchmarking CSVs as the other COPAL models, by reusing the COPAL
analysis machinery exactly as the RALPH BLAST baseline does
(see ``baselines/ralph_blast/run_ralph_blast.py``).

The sequence-distance metric is the *same* global-alignment %% identity that originally
defined the no-leakage set (``src/copal/analysis/compute_seqid.py``), so the database /
query split is internally consistent.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

# Reuse the exact identity metric used to define the no-leakage set.
from copal.analysis.compute_seqid import compute_pairwise_identity_and_coverage

MODEL_KEY = "seqid"          # lowercase key in per_model_rankings.csv "model" column
MODEL_DISPLAY = "SeqID"      # display name in MODEL_NAME_MAP / RANKING_MODELS_ORDER
# Matches the official `refAHL_all_conf_no_leakage` set's exclude_suffixes.
EXCLUDE_SUFFIXES = ("LONG", "SHORT", "MED", "_R1", "_R2")


def _normalize_luxr(refahl_id: str) -> str:
    """Same convention as the analysis module: stripped uppercase."""
    return str(refahl_id).strip().upper()


def build_global_aligner() -> Align.PairwiseAligner:
    """BLOSUM62 global aligner, matching compute_seqid.py exactly."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def closest_leaked_anchor(
    aligner: Align.PairwiseAligner,
    query_seq: str,
    db: list[tuple[str, str, str]],
) -> tuple[str | None, str | None, float]:
    """Return (db_id, anchor_ahl, identity) for the highest-identity database row.

    ``db`` is a list of (db_id, db_seq, anchor_ahl) tuples, pre-sorted by db_id so the
    argmax is a deterministic tie-break (first id wins on equal identity).
    """
    best_ident = -1.0
    best_id: str | None = None
    best_ahl: str | None = None
    for db_id, db_seq, anchor_ahl in db:
        ident, _cov = compute_pairwise_identity_and_coverage(aligner, query_seq, db_seq)
        if ident > best_ident:
            best_ident = ident
            best_id = db_id
            best_ahl = anchor_ahl
    return best_id, best_ahl, best_ident


def select_database(
    labels: pd.DataFrame,
    no_leak_ids: set[str],
    pool_set: set[str],
    db_scope: str,
    query_ids: list[str] | None,
) -> pd.DataFrame:
    """Return the anchor-database rows (``_rid``, ``seq``, ``true_ahl_id``), sorted by ``_rid``.

    Benchmark mode (no ``query_ids``) keeps the original semantics: the database is the
    *leaked* receptors (NOT in the no-leakage set), optionally restricted to RefAHL split.

    Inference mode (explicit ``query_ids``) uses inclusive scopes — the no-leakage set is a
    valid anchor source because leakage-vs-the-RefAHL-benchmark is irrelevant for out-of-
    benchmark queries: ``refahl_no_leak`` = RefAHL receptors IN the no-leakage set;
    ``refahl_leak`` = RefAHL receptors NOT in the no-leakage set; ``refahl_all`` = all RefAHL
    receptors; ``all`` = every receptor. It always drops the query receptors so a receptor
    can never be its own anchor. In both modes, rows whose cognate AHL is not a real ligand
    in the pool are dropped (removes unknown/pending).
    """
    df = labels.copy()
    is_ref = df["split"].astype(str).str.upper() == "REFAHL"
    if query_ids:
        if db_scope == "refahl_no_leak":
            df = df[is_ref & df["_rid"].isin(no_leak_ids)]
        elif db_scope == "refahl_leak":
            df = df[is_ref & ~df["_rid"].isin(no_leak_ids)]
        elif db_scope == "refahl_all":
            df = df[is_ref]
        elif db_scope == "all":
            pass  # every receptor, any split
        else:
            raise SystemExit(
                f"--db-scope {db_scope!r} is invalid with --query-ids; use one of "
                "refahl_no_leak | refahl_leak | refahl_all | all."
            )
        df = df[~df["_rid"].isin(set(query_ids))]
    else:
        if db_scope not in ("refahl", "all"):
            raise SystemExit(
                f"--db-scope {db_scope!r} requires --query-ids; benchmark mode supports "
                "only refahl | all."
            )
        df = df[~df["_rid"].isin(no_leak_ids)]
        if db_scope == "refahl":
            df = df[is_ref]
    df = df[df["true_ahl_id"].astype(str).isin(pool_set)]
    return df.sort_values("_rid")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Sequence-identity nearest-leak baseline: closest leaked sequence "
        "-> its AHL -> MCS-penalized ranking; writes COPAL benchmarking CSVs.",
    )
    parser.add_argument("--labels", type=Path, required=True,
                        help="paper_luxrs.csv (refahl_id, seq, split, true_ahl_id, ...).")
    parser.add_argument("--no-leakage-csv", type=Path, required=True,
                        help="curated_luxR_RefAHL_no_leakage.csv (defines the leak-free set).")
    parser.add_argument("--ligands-path", type=Path, required=True,
                        help="curated_ligands_annotated.csv (AHLid, abbreviation, subcategory, smi).")
    parser.add_argument("--sim-matrix", type=Path, required=True,
                        help="Precomputed chem-sim matrix CSV (e.g. mcs_penalized.csv).")
    parser.add_argument("--output-root", type=Path, required=True,
                        help="Root for full_rankings/ and ranking_performance/.")
    parser.add_argument("--workflow-dir", type=Path, default=None,
                        help="Dir for per_model_rankings/ensemble CSVs "
                             "(default: output_root/../workflow_seqid).")
    parser.add_argument("--db-scope",
                        choices=["refahl", "all", "refahl_no_leak", "refahl_leak",
                                 "refahl_all"],
                        default="refahl",
                        help="Anchor-database scope. Benchmark mode (no --query-ids): 'refahl' = "
                             "RefAHL-split leaked receptors only (strict); 'all' = every "
                             "non-no-leakage receptor with a usable in-pool AHL. Inference mode "
                             "(with --query-ids/--query-csv) uses inclusive scopes instead: "
                             "'refahl_no_leak' = RefAHL receptors IN the no-leakage set; "
                             "'refahl_leak' = RefAHL receptors NOT in the no-leakage set; "
                             "'refahl_all' = all RefAHL receptors; 'all' = every receptor (any "
                             "split). Query receptors are always excluded from the database "
                             "(no self-matching). Default: refahl.")
    parser.add_argument("--query-ids", type=str, default=None,
                        help="Comma-separated refahl_id list to use as the query set "
                             "(case-insensitive). Switches to inference mode: the inclusive "
                             "--db-scope semantics apply and only per_model_rankings.csv / "
                             "ensemble_rankings.csv are written (the RefAHL conf{1,2,3} "
                             "performance phase is skipped). Add --write-xlsx to also emit "
                             "full_rankings/all_splits.csv and SeqID.xlsx.")
    parser.add_argument("--query-csv", type=Path, default=None,
                        help="CSV with a 'refahl_id' column whose ids form the query set "
                             "(case-insensitive). Merged with --query-ids if both are given.")
    parser.add_argument("--write-xlsx", action="store_true",
                        help="In inference mode (--query-ids/--query-csv), also run the "
                             "full-ranking generate phase to emit full_rankings/all_splits.csv "
                             "and SeqID.xlsx for the query receptors (the RefAHL conf{1,2,3} "
                             "performance tables are still skipped). No effect in benchmark mode, "
                             "which always writes them.")
    args = parser.parse_args()

    # ---- Resolve the (optional) explicit query set -------------------------
    explicit_query_ids: list[str] | None = None
    if args.query_ids or args.query_csv:
        raw_ids: list[str] = []
        if args.query_ids:
            raw_ids += args.query_ids.split(",")
        if args.query_csv:
            qcsv = pd.read_csv(args.query_csv)
            if "refahl_id" not in qcsv.columns:
                raise SystemExit(f"--query-csv {args.query_csv} has no 'refahl_id' column.")
            raw_ids += qcsv["refahl_id"].astype(str).tolist()
        # Normalize with the same convention as the labels (_normalize_luxr == strip().upper()),
        # de-duplicate, and keep input order.
        seen: set[str] = set()
        explicit_query_ids = []
        for raw in raw_ids:
            nid = _normalize_luxr(raw)
            if nid and nid not in seen:
                seen.add(nid)
                explicit_query_ids.append(nid)
        if not explicit_query_ids:
            raise SystemExit("--query-ids/--query-csv resolved to an empty id list.")

    if args.workflow_dir is None:
        # Inference mode: keep each scope's rankings under its own --output-root so sweeping
        # scopes never collides. Benchmark mode default is unchanged.
        args.workflow_dir = (
            args.output_root if explicit_query_ids
            else args.output_root.parent / "workflow_seqid"
        )
    args.workflow_dir = args.workflow_dir.resolve()
    args.output_root = args.output_root.resolve()

    # ---- Load inputs -------------------------------------------------------
    labels = pd.read_csv(args.labels)
    labels = labels[labels["seq"].notna() & (labels["seq"].astype(str).str.strip() != "")].copy()
    labels["_rid"] = labels["refahl_id"].apply(_normalize_luxr)

    sim = pd.read_csv(args.sim_matrix, index_col=0)
    sim.index = sim.index.map(str)
    sim.columns = sim.columns.map(str)
    pool_ahls = list(sim.index)  # the candidate ligand pool (matrix rows/cols)
    pool_set = set(pool_ahls)

    no_leak_ids = set(
        pd.read_csv(args.no_leakage_csv)["refahl_id"].apply(_normalize_luxr)
    )

    # ---- Database (anchor source) -----------------------------------------
    db_df = select_database(labels, no_leak_ids, pool_set, args.db_scope, explicit_query_ids)
    db = list(zip(db_df["_rid"], db_df["seq"].astype(str).str.strip(),
                  db_df["true_ahl_id"].astype(str), strict=True))
    if not db:
        raise RuntimeError(f"Empty anchor database for --db-scope {args.db_scope}.")

    # ---- Query set --------------------------------------------------------
    if explicit_query_ids:
        qset = set(explicit_query_ids)
        query_df = labels[labels["_rid"].isin(qset)].copy().sort_values("_rid")
        missing = [q for q in explicit_query_ids if q not in set(query_df["_rid"])]
        if missing:
            print(f"WARNING: {len(missing)} query id(s) absent from labels "
                  f"(missing/blank seq or not found): {', '.join(missing)}")
        if query_df.empty:
            raise RuntimeError("Empty query set (no labels matched --query-ids/--query-csv).")
    else:
        query_df = labels[labels["_rid"].isin(no_leak_ids)].copy().sort_values("_rid")
        if query_df.empty:
            raise RuntimeError("Empty query set (no labels matched the no-leakage ids).")

    print(f"db-scope={args.db_scope}: {len(db)} database receptors, "
          f"{len(query_df)} query receptors, {len(pool_ahls)} candidate AHLs.")

    aligner = build_global_aligner()

    ranking_rows: list[dict] = []
    ensemble_rows: list[dict] = []
    for i, (_, q) in enumerate(query_df.iterrows(), start=1):
        luxr = _normalize_luxr(q["refahl_id"])
        qseq = str(q["seq"]).strip()
        db_id, anchor, ident = closest_leaked_anchor(aligner, qseq, db)
        if anchor is None:
            continue
        scores = sim.loc[anchor]                              # similarity to every pool AHL
        ranks = scores.rank(ascending=False, method="min")   # 1 = most similar (anchor itself)
        for ahl in pool_ahls:
            ahl_u = _normalize_luxr(ahl)
            rank = int(ranks[ahl])
            score = float(scores[ahl])
            # Column set mirrors RALPH so analysis.normalize_inputs consumes it unchanged.
            ranking_rows.append(
                {"luxr": luxr, "model": MODEL_KEY, "ahl": ahl_u, "rank": rank, "iptm": score}
            )
            ensemble_rows.append(
                {"luxr_ahl": f"{luxr}_{ahl_u}", "ensemble_rank": rank, "ensemble_score": score}
            )
        if i % 10 == 0:
            print(f"  [{i}/{len(query_df)}] {luxr}: closest={db_id} "
                  f"(ident={ident:.3f}) anchor={anchor}")

    rankings_df = pd.DataFrame(ranking_rows)
    ensemble_df = pd.DataFrame(ensemble_rows)
    if rankings_df.empty:
        raise RuntimeError("No rankings produced.")

    args.workflow_dir.mkdir(parents=True, exist_ok=True)
    rankings_path = args.workflow_dir / "per_model_rankings.csv"
    ensemble_path = args.workflow_dir / "ensemble_rankings.csv"
    rankings_df.to_csv(rankings_path, index=False)
    ensemble_df.to_csv(ensemble_path, index=False)
    print(f"Wrote {rankings_path} ({len(rankings_df)} rows)")
    print(f"Wrote {ensemble_path} ({len(ensemble_df)} rows)")

    # ---- Inference mode: rankings only, skip the benchmark generate phase --
    # The conf{1,2,3} / full-ranking / xlsx phase below assumes the RefAHL no-leakage
    # benchmark (and its EXCLUDE_SUFFIXES, which would drop e.g. NJ3_R1). For an explicit
    # out-of-benchmark query list it does not apply, so we stop after the rankings.
    if explicit_query_ids and not args.write_xlsx:
        print(f"Inference mode (db-scope={args.db_scope}): wrote rankings for "
              f"{rankings_df['luxr'].nunique()} query receptors; skipping the RefAHL "
              "conf{1,2,3} / full-ranking / xlsx generate phase. (Add --write-xlsx for "
              "full_rankings/all_splits.csv + SeqID.xlsx.)")
        return

    # ---- Generate phase: reuse COPAL analysis (mirrors run_ralph_blast.py) --
    from copal.analysis import pipeline as analysis

    analysis.MODEL_NAME_MAP[MODEL_KEY] = MODEL_DISPLAY
    if MODEL_DISPLAY not in analysis.RANKING_MODELS_ORDER:
        analysis.RANKING_MODELS_ORDER = list(analysis.RANKING_MODELS_ORDER) + [MODEL_DISPLAY]

    labels_raw = pd.read_csv(args.labels)
    ligands_raw = pd.read_csv(args.ligands_path)
    rankings_raw = pd.read_csv(rankings_path)
    ensemble_raw = pd.read_csv(ensemble_path)

    labels_norm, rankings_norm, ensemble_norm, ligands_norm = analysis.normalize_inputs(
        labels_raw, rankings_raw, ensemble_raw, ligands_raw
    )
    if explicit_query_ids:
        # Inference mode: keep exactly the requested receptors. Do NOT apply
        # EXCLUDE_SUFFIXES (it would drop e.g. NJ3_R1, which the user asked for).
        labels_norm = labels_norm[
            labels_norm["LuxR_name"].isin(set(explicit_query_ids))
        ].copy()
    else:
        # Evaluate on the leak-free benchmark exactly: restrict to the no-leakage
        # receptors and apply the official suffix exclusions, matching
        # `set=refAHL_all_conf_no_leakage`. Without this, leaked RefAHL receptors we
        # never ranked would show up as NaN rows in the performance table.
        labels_norm = labels_norm[labels_norm["LuxR_name"].isin(no_leak_ids)].copy()
        labels_norm = labels_norm[
            ~labels_norm["LuxR_name"].str.endswith(EXCLUDE_SUFFIXES, na=False)
        ].copy()

    full_dir = args.output_root / "full_rankings"
    full_dir.mkdir(parents=True, exist_ok=True)

    # Ranking performance: RefAHL conf 1, 2, 3 — only meaningful for the benchmark set.
    if not explicit_query_ids:
        ranking_dir = args.output_root / "ranking_performance"
        ranking_dir.mkdir(parents=True, exist_ok=True)
        for conf in [1, 2, 3]:
            subset = labels_norm[
                (labels_norm["split"].str.upper() == "REFAHL")
                & (labels_norm["confidence"] == conf)
            ].copy()
            subset = subset.drop_duplicates(subset=["LuxR_name", "split", "confidence"])
            true_rank_lookup = analysis.build_true_rank_lookup(
                subset, rankings_norm, ensemble_norm
            )
            perf_df = analysis.build_ranking_performance(subset, true_rank_lookup)
            perf_path = ranking_dir / f"ranking_performance_conf{conf}.csv"
            analysis.write_results_csv(perf_df, perf_path)
            n_scored = perf_df[MODEL_DISPLAY].notna().sum() if MODEL_DISPLAY in perf_df else 0
            print(f"Wrote {perf_path} ({len(perf_df)} LuxRs, {n_scored} with a true rank)")

    # Full rankings (all splits present in the query set).
    full_rankings_df = analysis.build_full_rankings(
        labels_norm, rankings_norm, ensemble_norm, ligands_norm
    )
    full_path = full_dir / "all_splits.csv"
    analysis.write_results_csv(full_rankings_df, full_path)
    print(f"Wrote {full_path}")

    # Stylized Excel workbook (mpKT placeholder, as RALPH does).
    luxrs_in_rankings = rankings_norm["luxr"].unique()
    mpkt_df = pd.DataFrame({"LuxR_name": luxrs_in_rankings, "mpKT": np.nan})
    excel_path = full_dir / f"{MODEL_DISPLAY}.xlsx"
    analysis.build_ranked_excel(
        MODEL_KEY, rankings_norm, ensemble_norm, labels_norm, ligands_norm, mpkt_df, excel_path
    )
    print(f"Wrote {excel_path}")


if __name__ == "__main__":
    main()
