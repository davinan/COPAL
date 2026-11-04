#!/usr/bin/env python3
"""Build the shipped RALPH-BLAST inference database (compressed sources).

Authoritative inputs (external to the repo; see baselines/ralph_blast/README.md):

  RALPH rankings/scores/mpKT  -- Supp Table 6 workbook, 5 sheets
      (COPAL, AF3, RF3, Boltz2, Boltz2Affinity); per RALPH, 3 rows:
      row0 = AHL id at rank i, row1 = abbreviation, row2 = score at rank i.
  RALPH sequences             -- ralph_v3.csv (luxr_name -> sequence, mmseqs_entry
      which holds the source FASTA header / organism annotation).
  Ligands                     -- curated_ligands_annotated.csv (AHLid, abbreviation, smi).
  RefAHL rankings/scores      -- benchmark_data_all_rankings.csv +
      benchmark_data_all_ensemble_rankings.csv (long per-model iptm/rank),
      restricted to split=='RefAHL' via paper_luxrs.csv (refahl_id, seq,
      split, true_ahl_id).

Emits into ralph_database/ (git-committed, compressed):
  ralph_scores.csv.xz         -- wide: luxr_name, source, organism, sequence,
      true_ahl_id, mpKT, {model}_1..58 (AHL ids) for model in
      ensemble/af3/rf3/boltz2/boltz2affinity. RALPH + RefAHL rows.
  ralph_model_scores.csv.xz   -- per-pair COPAL scores:
      luxr_name, ahl, copal (ensemble sum-of-ranks).
  ralph_db.faa.xz             -- BLAST FASTA, headers = luxr_name (RALPH id /
      RefAHL id), so blastp sseqid matches the scores key.
  ligands_58.csv              -- AHLid, abbreviation, smi.

This is a provenance/build tool. End users do NOT run it (its inputs are the
external paper artifacts); they use the shipped compressed files + setup_db.sh.
"""
from __future__ import annotations

import argparse
import lzma
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

SHEET_TO_MODEL = {
    "COPAL": "ensemble",
    "AF3": "af3",
    "RF3": "rf3",
    "Boltz2": "boltz2",
    "Boltz2Affinity": "boltz2affinity",
}
MODELS = ["ensemble", "af3", "rf3", "boltz2", "boltz2affinity"]
RANK_I = list(range(1, 59))  # 58 ligands
RANK_COLS = [f"rank_{i}" for i in RANK_I]


def _norm(s) -> str:
    return str(s).strip().upper()


# ---------------------------------------------------------------- RALPH (xlsx)
def parse_supp_xlsx(xlsx_path: Path):
    """Return (wide_rank, long_score, mpkt) for RALPH from Supp Table 6.

    wide_rank:  DataFrame indexed by luxr_name, cols {model}_1..58 (AHL ids).
    long_score: DataFrame [luxr_name, ahl, <model>] one col per model score.
    mpkt:       Series luxr_name -> mpKT.
    """
    wide = None
    score_frames = []
    mpkt = None
    for sheet, model in SHEET_TO_MODEL.items():
        df = pd.read_excel(xlsx_path, sheet_name=sheet, header=0)
        df = df[df["RALPH_id"].notna()].copy()
        # 3 rows per RALPH: 0=AHLid, 1=abbrev, 2=score (row order preserved)
        df["ri"] = df.groupby("RALPH_id", sort=False).cumcount()
        idrows = df[df.ri == 0].set_index("RALPH_id")
        scorerows = df[df.ri == 2].set_index("RALPH_id")

        # rankings: {model}_i = AHL id at rank i
        ranks = idrows[RANK_COLS].copy()
        ranks.columns = [f"{model}_{i}" for i in RANK_I]
        wide = ranks if wide is None else wide.join(ranks, how="outer")

        if mpkt is None:
            mpkt = idrows["mpKT"].astype(float)

        # scores: pair each AHL id (idrow) with its score (scorerow), per RALPH
        ids = idrows[RANK_COLS]
        scs = scorerows[RANK_COLS].apply(pd.to_numeric, errors="coerce")
        long = pd.DataFrame(
            {
                "luxr_name": np.repeat(ids.index.values, 58),
                "ahl": ids.values.reshape(-1),
                model: scs.values.reshape(-1),
            }
        )
        score_frames.append(long.set_index(["luxr_name", "ahl"]))

    long_score = pd.concat(score_frames, axis=1).reset_index()
    return wide, long_score, mpkt


def build_ralph(xlsx_path: Path, v3_csv: Path):
    wide, long_score, mpkt = parse_supp_xlsx(xlsx_path)
    v3 = pd.read_csv(v3_csv, usecols=["luxr_name", "mmseqs_entry", "sequence"])
    v3 = v3.drop_duplicates("luxr_name").set_index("luxr_name")

    ralphs = wide.index
    missing = set(ralphs) - set(v3.index)
    assert not missing, f"{len(missing)} RALPHs lack a sequence in ralph_v3"

    meta = pd.DataFrame(index=ralphs)
    meta["source"] = "RALPH"
    meta["organism"] = v3.loc[ralphs, "mmseqs_entry"].astype(str).str.lstrip(">")
    meta["sequence"] = v3.loc[ralphs, "sequence"].astype(str).str.strip()
    meta["true_ahl_id"] = pd.NA
    meta["mpKT"] = mpkt.reindex(ralphs)

    wide_out = meta.join(wide)
    wide_out.index.name = "luxr_name"
    wide_out = wide_out.reset_index()
    long_score = long_score.rename(columns={"ensemble": "copal"})
    return wide_out, long_score


# ------------------------------------------------------------------- RefAHL
def build_refahl(bench_rank_csv: Path, bench_ens_csv: Path, paper_luxrs_csv: Path):
    pl = pd.read_csv(paper_luxrs_csv)
    pl["key"] = pl["refahl_id"].map(_norm)
    refahl = pl[pl["split"].astype(str).str.upper() == "REFAHL"].copy()
    keep = set(refahl["key"])

    br = pd.read_csv(bench_rank_csv)  # rank, luxr, luxr_ahl, ahl, model, iptm,...
    br["key"] = br["luxr"].map(_norm)
    br = br[br["key"].isin(keep)].copy()

    be = pd.read_csv(bench_ens_csv)  # luxr_ahl, ensemble_rank, ensemble_score
    be["luxr"] = be["luxr_ahl"].str.rsplit("_", n=1).str[0]
    be["ahl"] = be["luxr_ahl"].str.rsplit("_", n=1).str[1]
    be["key"] = be["luxr"].map(_norm)
    be = be[be["key"].isin(keep)].copy()

    # wide rankings per model: {model}_i = ahl at rank i
    model_map = {"af3": "af3", "rf3": "rf3", "boltz2": "boltz2",
                 "boltz2affinity": "boltz2affinity"}
    wide_parts = {}
    long_parts = {}
    for m_raw, m in model_map.items():
        sub = br[br["model"] == m_raw]
        piv = (
            sub.sort_values("rank")
            .groupby("key")["ahl"]
            .apply(list)
        )
        wide_parts[m] = pd.DataFrame(
            {f"{m}_{i}": piv.map(lambda v, i=i: v[i - 1] if len(v) >= i else np.nan)
             for i in RANK_I}
        )
        # boltz2affinity is scored by its affinity head (iptm is NaN for it);
        # the other models by ipTM. Both rank by their own score.
        score_src = "affinity_probability_binary" if m == "boltz2affinity" else "iptm"
        long_parts[m] = sub.set_index(["key", "ahl"])[score_src].rename(m)

    # ensemble from benchmark ensemble rankings
    epiv = be.sort_values("ensemble_rank").groupby("key")["ahl"].apply(list)
    wide_parts["ensemble"] = pd.DataFrame(
        {f"ensemble_{i}": epiv.map(lambda v, i=i: v[i - 1] if len(v) >= i else np.nan)
         for i in RANK_I}
    )
    long_parts["copal"] = be.set_index(["key", "ahl"])["ensemble_score"].rename("copal")

    wide = pd.concat(wide_parts.values(), axis=1)
    # mpKT across the 4 model rank vectors (label-free, same recipe as RALPH)
    from scipy.stats import kendalltau

    rk = br.pivot_table(index=["key", "ahl"], columns="model", values="rank")
    mcols = [c for c in ["af3", "rf3", "boltz2", "boltz2affinity"] if c in rk.columns]

    def _mpkt(g):
        taus = []
        for i in range(len(mcols)):
            for j in range(i + 1, len(mcols)):
                t, _ = kendalltau(g[mcols[i]], g[mcols[j]])
                taus.append(t)
        return float(np.mean(taus)) if taus else np.nan

    mpkt = rk.groupby(level="key").apply(_mpkt)

    meta_lu = refahl.drop_duplicates("key").set_index("key")
    keys = wide.index
    meta = pd.DataFrame(index=keys)
    meta["source"] = "RefAHL"
    meta["organism"] = meta_lu.loc[keys, "strain_or_metagenome"].astype(str) \
        if "strain_or_metagenome" in meta_lu.columns else ""
    meta["sequence"] = meta_lu.loc[keys, "seq"].astype(str).str.strip()
    meta["true_ahl_id"] = meta_lu.loc[keys, "true_ahl_id"].astype(str) \
        if "true_ahl_id" in meta_lu.columns else pd.NA
    meta["mpKT"] = mpkt.reindex(keys)

    wide_out = meta.join(wide).reset_index().rename(columns={"key": "luxr_name"})

    long = pd.concat(long_parts.values(), axis=1).reset_index()
    long = long.rename(columns={"key": "luxr_name"})
    return wide_out, long


# --------------------------------------------------------------------- checks
def assert_wide(wide: pd.DataFrame, ligand_ids: set, tag: str):
    for m in MODELS:
        cols = [f"{m}_{i}" for i in RANK_I]
        block = wide[cols]
        # every row: a permutation of the 58 ligand ids, no dups/nulls
        for _, row in block.iterrows():
            vals = list(row)
            assert not any(pd.isna(v) for v in vals), f"{tag}:{m}: null in row"
            assert len(set(vals)) == 58, f"{tag}:{m}: not 58 unique in a row"
            assert set(vals) <= ligand_ids, f"{tag}:{m}: unknown AHL id"
    print(f"  [{tag}] wide asserts pass ({len(wide)} receptors x {len(MODELS)} models)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    # External (paper-side) inputs: required, no machine-local defaults.
    ap.add_argument("--xlsx", type=Path, required=True,
                    help="Supp Table 6 RALPH rankings workbook (.xlsx)")
    ap.add_argument("--ralph-v3", type=Path, required=True,
                    help="ralph_v3.csv (luxr_name -> sequence, mmseqs_entry)")
    ap.add_argument("--bench-rankings", type=Path, required=True,
                    help="benchmark_data_all_rankings.csv (RefAHL per-model ranks/iptm)")
    ap.add_argument("--bench-ensemble", type=Path, required=True,
                    help="benchmark_data_all_ensemble_rankings.csv (RefAHL ensemble)")
    # Repo-relative inputs: sensible defaults (run from repo root).
    ap.add_argument("--ligands", type=Path,
                    default=Path("data/curated_ligands_annotated.csv"))
    ap.add_argument("--paper-luxrs", type=Path, default=Path("data/paper_luxrs.csv"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("baselines/ralph_blast/ralph_database"))
    ap.add_argument("--no-refahl", action="store_true")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ligands
    lig = pd.read_csv(args.ligands)
    lig = lig[["AHLid", "abbreviation", "smi"]].dropna(subset=["AHLid", "smi"])
    ligand_ids = set(lig["AHLid"])
    assert len(ligand_ids) == 58, f"expected 58 ligands, got {len(ligand_ids)}"

    print("Parsing RALPH Supp Table 6 ...")
    r_wide, r_long = build_ralph(args.xlsx, args.ralph_v3)
    assert_wide(r_wide.set_index("luxr_name"), ligand_ids, "RALPH")

    frames_wide = [r_wide]
    frames_long = [r_long]
    if not args.no_refahl:
        print("Building RefAHL (split==RefAHL) ...")
        f_wide, f_long = build_refahl(args.bench_rankings, args.bench_ensemble,
                                      args.paper_luxrs)
        assert_wide(f_wide.set_index("luxr_name"), ligand_ids, "RefAHL")
        # id disjointness
        overlap = set(r_wide["luxr_name"]) & set(f_wide["luxr_name"])
        assert not overlap, f"RALPH/RefAHL id collision: {list(overlap)[:5]}"
        frames_wide.append(f_wide)
        frames_long.append(f_long)

    wide = pd.concat(frames_wide, ignore_index=True)
    long = pd.concat(frames_long, ignore_index=True)[["luxr_name", "ahl", "copal"]]

    # global uniqueness + clean FASTA ids
    assert wide["luxr_name"].is_unique, "duplicate luxr_name in combined wide"
    assert not wide["luxr_name"].astype(str).str.contains(r"\s").any(), "whitespace in id"

    long["copal"] = pd.to_numeric(long["copal"], errors="raise")

    # column order for wide
    front = ["luxr_name", "source", "organism", "true_ahl_id", "mpKT", "sequence"]
    rank_cols = [f"{m}_{i}" for m in MODELS for i in RANK_I]
    wide = wide[front + rank_cols]

    # write compressed
    def write_xz(df, name):
        path = args.out_dir / name
        with lzma.open(path, "wt", preset=9) as fh:
            df.to_csv(fh, index=False)
        return path.stat().st_size

    s1 = write_xz(wide, "ralph_scores.csv.xz")
    s2 = write_xz(long, "ralph_model_scores.csv.xz")

    # FASTA (headers = luxr_name)
    faa = args.out_dir / "ralph_db.faa"
    with open(faa, "w") as fh:
        for _, row in wide.iterrows():
            fh.write(f">{row['luxr_name']}\n{row['sequence']}\n")
    with open(faa, "rb") as fin, lzma.open(str(faa) + ".xz", "wb", preset=9) as fout:
        shutil.copyfileobj(fin, fout)
    s3 = (args.out_dir / "ralph_db.faa.xz").stat().st_size
    faa.unlink()

    lig.to_csv(args.out_dir / "ligands_58.csv", index=False)

    total = (s1 + s2 + s3) / 1e6
    print(f"\nWrote to {args.out_dir}:")
    print(f"  ralph_scores.csv.xz       {s1/1e6:6.2f} MB  ({len(wide)} receptors)")
    print(f"  ralph_model_scores.csv.xz {s2/1e6:6.2f} MB  ({len(long)} pairs)")
    print(f"  ralph_db.faa.xz           {s3/1e6:6.2f} MB")
    print(f"  ligands_58.csv            ({len(lig)} ligands)")
    print(f"  --- total compressed: {total:.2f} MB ---")


if __name__ == "__main__":
    main()
