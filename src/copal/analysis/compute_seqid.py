#!/usr/bin/env python3
"""Compute sequence identity, alignment coverage, and BLAST metrics between
a set of *reference* protein structures (e.g. PDB-deposited LuxR-family
receptors) and a *query* set (e.g. RefAHL-split benchmark sequences).

For each query, optionally exclude rows whose closest reference shares both
high sequence identity (>= --threshold) AND the same cognate ligand
("leakage" filter for benchmark splits).

Inputs
------
--queries CSV   columns required: refahl_id (or receptor_id), seq
                optional: true_ahl_abbrev (used for leakage filter)
--references CSV columns required: refahl_id (or receptor_id), seq
                 optional: cocrystal_ligand (used for leakage filter)

Outputs
-------
--out-table       wide table of pairwise identity / coverage / BLAST hits
--out-leakage     per-query "no-leakage" subset (only if leakage filter ran)

BLAST binary
------------
Set --blast-bin to a directory containing makeblastdb and blastp; defaults
to the parent dir of `which blastp`.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices


def _resolve_blast_bin(explicit: str | None) -> Path:
    if explicit:
        return Path(explicit)
    found = shutil.which("blastp")
    if found is None:
        raise RuntimeError(
            "blastp not on PATH; install BLAST+ or pass --blast-bin."
        )
    return Path(found).parent


def run(cmd, **kwargs):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kwargs)


def compute_pairwise_identity_and_coverage(aligner, seq1: str, seq2: str):
    alignment = aligner.align(seq1, seq2)[0]
    s1, s2 = alignment[0], alignment[1]
    identical = sum(1 for a, b in zip(s1, s2) if a == b and a != "-")
    aln_len = len(s1)
    aligned = sum(1 for a, b in zip(s1, s2) if a != "-" and b != "-")
    coverage = aligned / max(len(seq1), len(seq2))
    return identical / aln_len, coverage


def write_fasta(records, path: Path):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def run_blast(query_fasta: Path, subject_fasta: Path, workdir: Path, blast_bin: Path):
    db_path = workdir / "subject_db"
    run([
        str(blast_bin / "makeblastdb"),
        "-in", str(subject_fasta),
        "-dbtype", "prot",
        "-out", str(db_path),
    ])

    out_path = workdir / "blast_out.tsv"
    run([
        str(blast_bin / "blastp"),
        "-query", str(query_fasta),
        "-db", str(db_path),
        "-outfmt", "6",
        "-evalue", "1e-3",
        "-max_target_seqs", "100",
        "-out", str(out_path),
    ])

    hits: dict = {}
    if out_path.exists() and out_path.stat().st_size > 0:
        cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        df = pd.read_csv(out_path, sep="\t", header=None, names=cols)
        df = df.sort_values("bitscore", ascending=False).drop_duplicates(
            subset=["qseqid", "sseqid"], keep="first"
        )
        for _, row in df.iterrows():
            hits[(row["qseqid"], row["sseqid"])] = {
                "evalue": float(row["evalue"]),
                "bitscore": float(row["bitscore"]),
            }
    return hits


def _normalize_id_column(df: pd.DataFrame, candidates=("refahl_id", "receptor_id")) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"CSV missing identifier column (tried {candidates}).")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--queries", type=Path, required=True,
                    help="CSV of query sequences (refahl_id|receptor_id, seq).")
    ap.add_argument("--references", type=Path, required=True,
                    help="CSV of reference sequences (refahl_id|receptor_id, seq, optional cocrystal_ligand).")
    ap.add_argument("--out-table", type=Path, required=True,
                    help="Output wide CSV of pairwise metrics.")
    ap.add_argument("--out-leakage", type=Path, default=None,
                    help="If set, write a no-leakage subset CSV.")
    ap.add_argument("--threshold", type=float, default=0.40,
                    help="Identity threshold for leakage exclusion (default 0.40).")
    ap.add_argument("--ligand-col-query", default="true_ahl_abbrev",
                    help="Column in queries CSV holding the cognate ligand abbreviation.")
    ap.add_argument("--ligand-col-reference", default="cocrystal_ligand",
                    help="Column in references CSV holding the co-crystallized ligand.")
    ap.add_argument("--exclude-suffix-rx", default=r"_(short|med|long)$",
                    help="Regex; query IDs matching it are excluded (case-insensitive).")
    ap.add_argument("--blast-bin", default=None,
                    help="Directory containing blastp/makeblastdb (default: dir of `which blastp`).")
    args = ap.parse_args()

    blast_bin = _resolve_blast_bin(args.blast_bin)

    queries_df = pd.read_csv(args.queries)
    refs_df = pd.read_csv(args.references)
    q_id_col = _normalize_id_column(queries_df)
    r_id_col = _normalize_id_column(refs_df)

    if args.exclude_suffix_rx:
        suffix_rx = re.compile(args.exclude_suffix_rx, re.IGNORECASE)
        queries_df = queries_df[~queries_df[q_id_col].astype(str).str.contains(suffix_rx)].copy()

    queries_df = queries_df[queries_df["seq"].astype(str).str.len() > 0].copy()
    refs_df = refs_df[refs_df["seq"].astype(str).str.len() > 0].copy()

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        query_fa = tmp / "queries.fasta"
        subj_fa = tmp / "refs.fasta"
        write_fasta(zip(queries_df[q_id_col], queries_df["seq"]), query_fa)
        write_fasta(zip(refs_df[r_id_col], refs_df["seq"]), subj_fa)
        print("Running BLAST...")
        blast_hits = run_blast(query_fa, subj_fa, tmp, blast_bin)
        print(f"BLAST returned {len(blast_hits)} (query, subject) hits")

    print("Computing pairwise global alignments...")
    aln_metrics: dict = {}
    for qid, qseq in zip(queries_df[q_id_col], queries_df["seq"]):
        for rid, rseq in zip(refs_df[r_id_col], refs_df["seq"]):
            aln_metrics[(qid, rid)] = compute_pairwise_identity_and_coverage(aligner, qseq, rseq)

    print("Building wide identity table...")
    rows = []
    for _, q in queries_df.iterrows():
        qid = q[q_id_col]
        row: dict = {"query_id": qid}
        if args.ligand_col_query in queries_df.columns:
            row["query_ligand"] = q[args.ligand_col_query]
        for _, r in refs_df.iterrows():
            rid = r[r_id_col]
            ident, cov = aln_metrics[(qid, rid)]
            hit = blast_hits.get((qid, rid))
            row[f"{rid}_seq_id"] = round(ident, 4)
            row[f"{rid}_aln_coverage"] = round(cov, 4)
            row[f"{rid}_blast_evalue"] = hit["evalue"] if hit else float("inf")
            row[f"{rid}_blast_bitscore"] = hit["bitscore"] if hit else 0.0
            if args.ligand_col_reference in refs_df.columns and args.ligand_col_query in queries_df.columns:
                row[f"{rid}_same_ligand"] = (
                    str(q[args.ligand_col_query]) == str(r[args.ligand_col_reference])
                )
        rows.append(row)
    wide = pd.DataFrame(rows)
    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(args.out_table, index=False)
    print(f"Wrote {args.out_table} ({len(wide)} rows, {len(wide.columns)} cols)")

    if args.out_leakage is None:
        return

    if args.ligand_col_query not in queries_df.columns or args.ligand_col_reference not in refs_df.columns:
        raise SystemExit(
            f"--out-leakage requires '{args.ligand_col_query}' in queries and "
            f"'{args.ligand_col_reference}' in references."
        )

    print("Building leakage filter...")
    ref_records = list(refs_df[[r_id_col, args.ligand_col_reference]].itertuples(index=False, name=None))
    ref_id_set = {rid for rid, _ in ref_records}
    curated = []
    for _, q in queries_df.iterrows():
        qid = q[q_id_col]
        if qid in ref_id_set:
            continue
        max_ident = 0.0
        closest_id = None
        closest_lig = None
        leaked = False
        for rid, rlig in ref_records:
            ident, _ = aln_metrics[(qid, rid)]
            if ident > max_ident:
                max_ident = ident
                closest_id = rid
                closest_lig = rlig
            if ident >= args.threshold and str(q[args.ligand_col_query]) == str(rlig):
                leaked = True
        if leaked:
            continue
        curated.append({
            "query_id": qid,
            "cognate_ligand": q[args.ligand_col_query],
            "max_identity_to_any_reference": round(max_ident, 4),
            "closest_reference": closest_id,
            "closest_reference_ligand": closest_lig,
            "reason": (
                f"<{args.threshold:.2f} identity"
                if max_ident < args.threshold
                else f"≥{args.threshold:.2f} but different ligand"
            ),
        })
    args.out_leakage.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(curated).to_csv(args.out_leakage, index=False)
    print(f"Wrote {args.out_leakage} ({len(curated)} rows)")


if __name__ == "__main__":
    main()
