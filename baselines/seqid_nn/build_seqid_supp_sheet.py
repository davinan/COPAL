#!/usr/bin/env python3
"""Build a standalone Seq-ID ranking sheet matching An_Supp_Table3_rankings.xlsx.

For every receptor in the referenced supplementary workbook (all 138), this:
  1. pulls the receptor's metadata (cols A-G) straight from an existing model sheet
     so the Seq-ID sheet aligns row-for-row with AF3/Boltz2/RF3;
  2. runs the Seq-ID nearest-leak method (closest RefAHL-leaked sequence by
     global-alignment %% identity -> its cognate AHL as anchor -> rank all candidate
     AHLs by MCS bond-penalized similarity to the anchor);
  3. writes a standalone .xlsx whose single 'Seq-ID' sheet uses the exact supp-table
     3-rows-per-receptor layout (AHL ids / abbreviations / scores), green/orange mpKT
     fill, and a yellow fill on the true-AHL's rank position.

The leaked database is the RefAHL leaked set (db-scope=refahl, 50 receptors), identical
to baselines/seqid_nn/run_seqid_baseline.py benchmark mode, so the 103 leak-free
receptors reproduce that baseline exactly. Database members appear in the query set too,
so they self-match (rank-1 = own AHL), as requested.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import openpyxl
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from openpyxl.styles import Font, PatternFill

# Stable copal API (the same identity metric that defined the no-leakage set).
from copal.analysis.compute_seqid import compute_pairwise_identity_and_coverage

N_META = 8  # A..H: LuxR name, Category, conf, AHLpref, AHLpref abbrev, AHLpref category, mpKT, rank
N_RANKS = 58


def _norm(s) -> str:
    return str(s).strip().upper()


def build_global_aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    return a


def closest_leaked_anchor(aligner, qseq, db):
    """db: list of (db_id, db_seq, anchor_ahl) sorted by db_id. Returns (ident, db_id, anchor)."""
    best_ident, best_id, best_ahl = -1.0, None, None
    for db_id, db_seq, anchor_ahl in db:
        ident, _cov = compute_pairwise_identity_and_coverage(aligner, qseq, db_seq)
        if ident > best_ident:
            best_ident, best_id, best_ahl = ident, db_id, anchor_ahl
    return best_ident, best_id, best_ahl


def read_supp(supp_path: Path, source_sheet: str):
    """Return (receptors, ahl_abbrev). receptors: list of dicts with name + meta A..G.

    ahl_abbrev: {AHLid -> abbreviation} read from the supp sheet itself (row1 ids vs row2 abbrev).
    """
    wb = openpyxl.load_workbook(supp_path, read_only=True, data_only=True)
    ws = wb[source_sheet]
    rows = [r for r in ws.iter_rows(values_only=True)]
    wb.close()

    data = rows[1:]
    # Group into receptor blocks separated by blank (col-A-empty) rows.
    blocks, block = [], []
    for r in data:
        if r[0] in (None, ""):
            if block:
                blocks.append(block)
                block = []
            continue
        block.append(r)
    if block:
        blocks.append(block)

    receptors = []
    ahl_abbrev: dict[str, str] = {}
    for blk in blocks:
        head = blk[0]
        receptors.append({
            "name": head[0],
            "meta": list(head[0:7]),  # A..G (LuxR name .. mpKT)
        })
        if len(blk) >= 2:
            ids_row, abbr_row = blk[0], blk[1]
            for c in range(N_META, N_META + N_RANKS):
                aid, abbr = ids_row[c], abbr_row[c]
                if aid not in (None, "") and _norm(aid) not in ahl_abbrev:
                    ahl_abbrev[_norm(aid)] = abbr
    return receptors, ahl_abbrev


def build_leaked_db(labels: pd.DataFrame, no_leak_ids: set[str], pool_set: set[str]):
    """RefAHL-leaked database (db-scope=refahl), identical to run_seqid_baseline benchmark mode."""
    is_ref = labels["split"].astype(str).str.upper() == "REFAHL"
    df = labels[~labels["_rid"].isin(no_leak_ids) & is_ref]
    df = df[df["true_ahl_id"].astype(str).isin(pool_set)].sort_values("_rid")
    return list(zip(df["_rid"], df["seq"].astype(str).str.strip(),
                    df["true_ahl_id"].astype(str), strict=True))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--supp-xlsx", type=Path, required=True)
    ap.add_argument("--source-sheet", default="AF3")
    ap.add_argument("--labels", type=Path, required=True)
    ap.add_argument("--no-leakage-csv", type=Path, required=True)
    ap.add_argument("--sim-matrix", type=Path, required=True)
    ap.add_argument("--out-xlsx", type=Path, required=True)
    ap.add_argument("--sheet-name", default="Seq-ID")
    args = ap.parse_args()

    receptors, ahl_abbrev = read_supp(args.supp_xlsx, args.source_sheet)
    print(f"supp receptors: {len(receptors)} (sheet '{args.source_sheet}'); "
          f"abbrev map: {len(ahl_abbrev)} AHLs")

    sim = pd.read_csv(args.sim_matrix, index_col=0)
    sim.index = sim.index.map(_norm)
    sim.columns = sim.columns.map(_norm)
    pool_ahls = list(sim.index)
    pool_set = set(pool_ahls)
    assert len(pool_ahls) == N_RANKS, f"expected {N_RANKS} ligands, got {len(pool_ahls)}"

    labels = pd.read_csv(args.labels)
    labels = labels[labels["seq"].notna() & (labels["seq"].astype(str).str.strip() != "")].copy()
    labels["_rid"] = labels["refahl_id"].apply(_norm)
    seq_map = dict(zip(labels["_rid"], labels["seq"].astype(str).str.strip()))

    no_leak_ids = set(pd.read_csv(args.no_leakage_csv)["refahl_id"].apply(_norm))
    db = build_leaked_db(labels, no_leak_ids, pool_set)
    print(f"leaked DB (refahl scope): {len(db)} receptors")

    aligner = build_global_aligner()

    # ---- styling ----
    yellow = PatternFill(start_color="FFFF00", end_color="FFFF00", fill_type="solid")
    orange = PatternFill(start_color="FFA500", end_color="FFA500", fill_type="solid")
    green = PatternFill(start_color="92D050", end_color="92D050", fill_type="solid")

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = args.sheet_name
    header = [
        "LuxR name", "LuxR Category", "RefAHL confidence class", "AHLpref",
        "AHLpref abbrev", "AHLpref category", "mpKT", "AHLpref Seq-ID rank",
    ] + [f"Seq-ID rank_{k}" for k in range(1, N_RANKS + 1)]
    ws.append(header)
    for cell in ws[1]:
        cell.font = Font(bold=True)

    n_blank_seq = 0
    n_selfmatch = 0
    for rec in receptors:
        name = rec["name"]
        meta_ag = list(rec["meta"])  # A..G (7 values)
        ahlpref = _norm(meta_ag[3]) if meta_ag[3] not in (None, "") else None
        mpkt_val = meta_ag[6]
        rid = _norm(name)

        if rid in seq_map:
            _ident, db_id, anchor = closest_leaked_anchor(aligner, seq_map[rid], db)
            scores = sim.loc[anchor]
            # deterministic order: best score first, ties broken by AHLid
            sorted_ahls = sorted(pool_ahls, key=lambda a: (-float(scores[a]), a))
            sorted_scores = [float(scores[a]) for a in sorted_ahls]
            true_rank_pos = (sorted_ahls.index(ahlpref) + 1) if ahlpref in pool_set else None
            if db_id == rid:
                n_selfmatch += 1
        else:
            sorted_ahls = [None] * N_RANKS
            sorted_scores = [None] * N_RANKS
            true_rank_pos = None
            n_blank_seq += 1

        meta_vals = meta_ag + [true_rank_pos]  # A..H (8 values)
        row_ids = meta_vals + sorted_ahls
        row_abbr = meta_vals + [ahl_abbrev.get(a) if a is not None else None for a in sorted_ahls]
        row_score = meta_vals + sorted_scores
        for rv in (row_ids, row_abbr, row_score):
            ws.append(rv)

        cur = ws.max_row
        block_rows = (cur - 2, cur - 1, cur)
        # mpKT fill (col G = 7)
        if mpkt_val is not None and not (isinstance(mpkt_val, float) and math.isnan(mpkt_val)):
            try:
                fill = orange if float(mpkt_val) < 0.3 else green
                for r in block_rows:
                    ws.cell(row=r, column=7).fill = fill
            except (TypeError, ValueError):
                pass
        # yellow fill on the true-AHL rank position (col = N_META + pos)
        if true_rank_pos is not None:
            col_idx = N_META + true_rank_pos
            for r in block_rows:
                ws.cell(row=r, column=col_idx).fill = yellow

        ws.append([])  # blank separator

    ws.freeze_panes = "A2"
    args.out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    wb.save(args.out_xlsx)
    print(f"Wrote {args.out_xlsx}: {len(receptors)} receptors "
          f"({n_selfmatch} DB self-matches, {n_blank_seq} blank/no-seq)")


if __name__ == "__main__":
    main()
