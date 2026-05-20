#!/usr/bin/env python3
"""
LuxR benchmark analysis pipeline.

Runs three phases in order:

  1. Collect (optional): Gather raw benchmark records from basepath (and optionally
     af3_basepath), optionally computing AF3 ipTM from JSON or using base-folder CSV
     columns. Writes predictions.csv to the workflow directory.
     Skipped if --input-csv is provided (use pre-collected CSV instead).

  2. Postprocess: Read the processed CSV and labels; produce per_model_rankings.csv,
     ensemble_rankings.csv, ablation_table.csv, and done_pairs.json in the workflow directory.

  3. Generate: Read workflow outputs plus labels and ligands; write analysis CSVs
     (ranking_performance, acyl_moiety, chemical_similarity, mpkt, full_rankings,
     ablation) under --output-root.

Main CLI options:
  --input-csv       If provided, skip collection and use this file as the processed
                    benchmark CSV for postprocessing.
  --af3-iptm-source For collection: "base" = use AF3 ipTM from base-folder CSV columns;
                    "json" = compute ipTM from *_confidences.json (default).
  --workflow-dir    Directory for processed CSV (when collecting) and postprocess
                    outputs (rankings, ensemble, ablation, done.json).
  --output-root     Directory for generate outputs (ranking_performance, etc.).
"""
from __future__ import annotations

import argparse
import json
import re
from itertools import combinations
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr
from tqdm import tqdm


# -----------------------------------------------------------------------------
# Constants (used in generate phase and for model naming)
# -----------------------------------------------------------------------------

MODEL_NAME_MAP = {
    "af3": "AF3",
    "boltz2": "Boltz2",
    "boltz2affinity": "Boltz2Affinity",
    "rf3": "RF3",
}
RANKING_MODELS_ORDER = ["AF3", "Boltz2", "Boltz2Affinity", "RF3", "Ensemble"]
EXCLUDE_SUFFIXES = ("LONG", "SHORT", "MED")
SMILES_COL = "smi"

PROCESSED_CSV_NAME = "predictions.csv"
RANKINGS_CSV_NAME = "per_model_rankings.csv"
ENSEMBLE_CSV_NAME = "ensemble_rankings.csv"
DONE_PAIRS_NAME = "done_pairs.json"

_RDKIT = None


def load_no_leakage_names(path: Path) -> set[str]:
    """Load the no-leakage CSV and return a set of uppercased LuxR names."""
    df = pd.read_csv(path)
    return set(df["refahl_id"].astype(str).str.upper().str.strip())


# =============================================================================
# Phase 1: Collection
# =============================================================================


def compute_af3_iptm(confidences_json_path: Path) -> float:
    """Compute interface pTM from an AF3 confidences JSON (PAE + token_chain_ids)."""
    with confidences_json_path.open() as fh:
        data = json.load(fh)
    pae = np.array(data["pae"])
    token_chain_ids = np.array(data["token_chain_ids"])
    num_tokens = len(token_chain_ids)

    clipped_num_res = max(num_tokens, 19)
    d0 = 1.24 * (clipped_num_res - 15) ** (1.0 / 3) - 1.8
    tm_adjusted_pae = 1.0 / (1.0 + np.square(pae) / np.square(d0))

    interface_mask = (token_chain_ids[:, None] != token_chain_ids[None, :]).astype(float)
    if interface_mask.sum() == 0:
        return float("nan")

    normed = interface_mask / (1e-8 + interface_mask.sum(axis=-1, keepdims=True))
    per_alignment = (tm_adjusted_pae * normed).sum(axis=-1)
    return float(per_alignment.max())


def load_af3_iptm_from_json(model_dir: Path) -> dict[int, float]:
    """Load ipTM per batch from *_summary_confidences.json files in model_dir."""
    batch_to_iptm: dict[int, float] = {}
    for jf in sorted(model_dir.glob("*_summary_confidences.json")):
        match = re.search(r"_sample-(\d+)_summary_confidences\.json$", jf.name)
        if match is None:
            continue
        batch_idx = int(match.group(1))
        with jf.open() as fh:
            data = json.load(fh)
        batch_to_iptm[batch_idx] = float(data["chain_pair_iptm"][0][1])
    return batch_to_iptm


def collect_all_from_basepath(
    path: Path, af3_basepath: Path | None = None
) -> dict[str, dict[str, list[Path]]]:
    """Map run directory to (luxr_name -> ligand -> list of model dirs). AF3 can be from af3_basepath."""
    confs = list(path.glob("confidences/*/*"))
    name_to_confs: dict[str, dict[str, list[Path]]] = {}
    for c in confs:
        # Uppercase before splitting: AF3 lowercases its output dir names
        # (metha_001a_ahl03), so a case-sensitive split on "_AHL" would miss.
        dirname = c.parent.name.upper()
        name = dirname.split("_AHL")[0].split("_HEHEAA")[0]
        lig = dirname.split("_")[-1]
        name_to_confs.setdefault(name, {}).setdefault(lig, [])

        if c.name == "af3" and af3_basepath is not None:
            af3_candidate = af3_basepath / path.name / c.relative_to(path)
            if af3_candidate.exists():
                name_to_confs[name][lig].append(af3_candidate)
                continue
        name_to_confs[name][lig].append(c)
    return name_to_confs


def records_from_model_dir(
    name: str,
    lig: str,
    model_dir: Path,
    af3_iptm_source: Literal["base", "json"] = "json",
) -> list[dict[str, Any]]:
    """
    Build benchmark records from one model directory (confidences.csv or .score + optional metrics).
    For AF3: if af3_iptm_source is "base", use CSV columns only; if "json", compute ipTM from JSON.
    """
    model_name = model_dir.name.lower()
    luxr_ahl = f"{name}_{lig}"

    conf_csv = model_dir / "confidences.csv"
    score_files = list(model_dir.glob("*.score"))
    if len(score_files) == 1:
        cdf = pd.read_csv(score_files[0])
    elif conf_csv.exists():
        cdf = pd.read_csv(conf_csv)
    else:
        return []

    if model_name == "rf3" and "chain_i_interface" in cdf.columns:
        metrics_path = Path(str(score_files[0]).replace(".score", "_metrics.csv"))
        metrics_df = pd.read_csv(metrics_path)
        batch_to_iptm = {
            k: metrics_df[f"iptm.iptm_protein_ligand_{k}"].values[0] for k in range(5)
        }
        cdf = cdf[(cdf.chain_i_interface == "A_1") & (cdf.chain_j_interface == "B_1")]
    elif model_name == "af3":
        if af3_iptm_source == "json":
            batch_to_iptm = load_af3_iptm_from_json(model_dir)
        else:
            batch_to_iptm = {}
    else:
        batch_to_iptm = {}

    records: list[dict[str, Any]] = []
    for _, row in cdf.iterrows():
        row_dict = row.to_dict()
        if "pae_interface" not in row_dict:
            continue
        min_pae = row_dict.get("min_pae_interface")
        if model_name == "af3" and "min_pae_interface_af3" in row_dict:
            min_pae = row_dict["min_pae_interface_af3"]
        row_iptm = row_dict.get("iptm")
        if row_dict.get("batch_idx") is not None:
            batch_idx = int(row_dict["batch_idx"])
            if batch_idx in batch_to_iptm:
                row_iptm = batch_to_iptm[batch_idx]
        rec: dict[str, Any] = {
            "batch_idx": row_dict.get("batch_idx"),
            "structure_path": None,
            "luxr": name,
            "ahl": lig,
            "luxr_ahl": luxr_ahl,
            "pae_interface": row_dict.get("pae_interface"),
            "min_pae_interface": min_pae,
            "iptm": row_iptm,
            "seed": row_dict.get("seed"),
            "model": model_name,
        }
        for key, value in row_dict.items():
            if "affinity" in key:
                rec[key] = value
        records.append(rec)

    if model_name == "boltz2" and "model" in cdf.columns:
        for _, row in cdf.iterrows():
            row_dict = row.to_dict()
            if row_dict.get("model", "").lower() == "boltz2affinity" and "affinity_probability_binary" in row_dict:
                affinity_rec: dict[str, Any] = {
                    "batch_idx": row_dict.get("batch_idx"),
                    "structure_path": None,
                    "luxr": name,
                    "ahl": lig,
                    "luxr_ahl": luxr_ahl,
                    "pae_interface": None,
                    "min_pae_interface": None,
                    "iptm": None,
                    "seed": None,
                    "model": "boltz2affinity",
                }
                for key, value in row_dict.items():
                    if "affinity" in key:
                        affinity_rec[key] = value
                records.append(affinity_rec)
    return records


def collect_records(
    basepath: Path,
    af3_basepath: Path | None = None,
    af3_iptm_source: Literal["base", "json"] = "json",
) -> pd.DataFrame:
    """Collect all benchmark records from basepath (and optional af3_basepath) into one DataFrame."""
    paths = [p for p in basepath.glob("*") if "_to_" not in p.name]
    dfs: list[pd.DataFrame] = []
    for path in tqdm(paths, desc="Collecting benchmark records"):
        if not path.is_dir():
            continue
        name_to_confs = collect_all_from_basepath(path, af3_basepath=af3_basepath)
        subdirs = {
            model_dir.parent.parent.parent.name
            for name in name_to_confs
            for lig in name_to_confs[name]
            for model_dir in name_to_confs[name][lig]
        }
        for subdir in subdirs:
            for luxr in name_to_confs:
                records: list[dict[str, Any]] = []
                for lig in name_to_confs[luxr]:
                    for model_dir in name_to_confs[luxr][lig]:
                        if model_dir.parent.parent.parent.name != subdir:
                            continue
                        records.extend(
                            records_from_model_dir(
                                luxr, lig, model_dir, af3_iptm_source=af3_iptm_source
                            )
                        )
                if not records:
                    continue
                df = pd.DataFrame(records)
                if "structure_path" in df.columns:
                    df["structure_path"] = df["structure_path"].apply(
                        lambda x: str(x) if x is not None else None
                    )
                dfs.append(df)
    if not dfs:
        raise RuntimeError("No benchmark records were collected.")
    return pd.concat(dfs, ignore_index=True)


# =============================================================================
# Phase 2: Postprocess
# =============================================================================


def make_ranking(
    df: pd.DataFrame,
    col: str,
    ascending: bool = True,
    pool: str = "max",
    expected_ligand_count: int = 58,
) -> pd.DataFrame:
    """Rank luxr_ahl by col (after grouping by luxr_ahl and aggregating with pool).

    Guards against partial inputs by requiring exactly `expected_ligand_count`
    unique (luxr, ligand) pairs per receptor — defaults to 58 (the paper
    benchmark) but can be set to whatever the ligands CSV provides.
    """
    df = df.dropna(subset=[col])
    df = (
        df.groupby("luxr_ahl")
        .agg({col: pool, **{k: "first" for k in df.columns if k not in [col, "luxr_ahl"]}})
        .reset_index()
    )
    if df.luxr_ahl.nunique() != expected_ligand_count:
        raise ValueError(
            f"DataFrame contains {df.luxr_ahl.nunique()} unique luxr_ahl pairs, "
            f"expected {expected_ligand_count}"
        )
    if col not in df.columns:
        raise ValueError(f"Column {col} not found in DataFrame")
    if df[col].isna().any():
        raise ValueError(f"Column {col} contains NaN values")
    df = df.sort_values(col, ascending=ascending)
    # Ties at the per-model score level (e.g. exact ipTM ties) get the same
    # rank rather than a row-order-dependent break. This makes downstream
    # aggregation reproducible under input row reorderings.
    df["rank"] = df[col].rank(ascending=ascending, method="min")
    return df


def compute_all_rankings(
    all_df: pd.DataFrame, expected_ligand_count: int = 58
) -> pd.DataFrame:
    """Produce per-model rankings (iptm or affinity_probability_binary for boltz2affinity).

    Only ranks receptors that have a complete set of ligand predictions
    (`expected_ligand_count` unique luxr_ahl pairs) — this is the
    "all predictions present" guard. Set to the number of ligands in
    your ligands CSV.
    """
    rankings: list[pd.DataFrame] = []
    for model, mdf in all_df.groupby("model"):
        for _, ldf in mdf.groupby("luxr"):
            if ldf.luxr_ahl.nunique() != expected_ligand_count:
                continue
            col = "iptm" if model != "boltz2affinity" else "affinity_probability_binary"
            ranked_ldf = make_ranking(
                ldf, col, ascending=False, pool="max",
                expected_ligand_count=expected_ligand_count,
            )
            cols = ["rank", "luxr", "luxr_ahl", "ahl", "model", col]
            rankings.append(ranked_ldf[cols])
    if not rankings:
        raise RuntimeError(
            f"No rankings were produced. Every receptor had fewer than "
            f"{expected_ligand_count} unique luxr_ahl pairs — either the inference "
            f"is incomplete, or pass --expected-ligand-count to match your "
            f"actual ligand set."
        )
    return pd.concat(rankings, ignore_index=True)


AGGREGATIONS = ("rrf", "sum", "mean", "median", "min", "max")
DEFAULT_AGGREGATION = "sum"
DEFAULT_RRF_K = 60


def make_ensemble_ranking(
    rankings_df: pd.DataFrame,
    models_to_include: list[str],
    aggregation: str = DEFAULT_AGGREGATION,
    rrf_k: int = DEFAULT_RRF_K,
) -> pd.DataFrame:
    """Aggregate per-model ranks into an ensemble rank per (luxr, luxr_ahl).

    Aggregations (no fitting; benchmark-only):
      * ``rrf``    — Reciprocal Rank Fusion: ``score = sum_models 1 / (k + rank)``
                     (Cormack et al., CIKM 2009). Robust to a single bad model
                     and concentrates resolution at the top of each per-model
                     ranking, matching the top-k success metric used downstream.
      * ``sum`` / ``mean`` / ``median`` / ``min`` / ``max`` — corresponding
                     reduction over per-model ranks. ``mean`` is the right
                     choice when comparing across model-subset ablations
                     because it normalizes for the number of models included.

    For RRF the final rank is ascending in ``-score`` (higher score = better).
    For rank-based aggregations the final rank is ascending in the aggregated
    rank (lower rank = better). In both cases ties at the final step share the
    better rank (``method="min"``), which is reproducible under row reorderings
    and gives the standard "generous" interpretation of top-k membership.
    """
    if aggregation not in AGGREGATIONS:
        raise ValueError(
            f"Unknown aggregation {aggregation!r}; expected one of {AGGREGATIONS}"
        )
    df = rankings_df[rankings_df.model.isin(models_to_include)].copy()

    if aggregation == "rrf":
        df["_score"] = 1.0 / (rrf_k + df["rank"].astype(float))
        agg = (
            df.groupby(["luxr", "luxr_ahl"], as_index=False)["_score"].sum()
            .rename(columns={"_score": "ensemble_score"})
        )
        agg["ensemble_rank"] = agg.groupby("luxr")["ensemble_score"].rank(
            ascending=False, method="min"
        )
    else:
        agg = (
            df.groupby(["luxr", "luxr_ahl"], as_index=False)["rank"]
            .agg(aggregation)
            .rename(columns={"rank": "ensemble_score"})
        )
        agg["ensemble_rank"] = agg.groupby("luxr")["ensemble_score"].rank(
            ascending=True, method="min"
        )

    return agg[["luxr_ahl", "ensemble_rank", "ensemble_score"]]


def get_true_ahl_ranking(
    labels_df: pd.DataFrame, ensemble_df: pd.DataFrame
) -> pd.DataFrame:
    """For each label row, attach the ensemble rank of the true AHL (luxr + true_ahl_id)."""
    recs = []
    for _, row in labels_df.iterrows():
        luxr_ahl = f"{row['luxr']}_{row['true_ahl_id']}"
        ranking = ensemble_df[ensemble_df.luxr_ahl == luxr_ahl]
        if ranking.empty:
            continue
        recs.append({"true_pred_ranking": ranking["ensemble_rank"].values[0], **row.to_dict()})
    return pd.DataFrame(recs)


def build_ablation_table(
    rankings_df: pd.DataFrame,
    ablation_df: pd.DataFrame,
    aggregation: str = DEFAULT_AGGREGATION,
    rrf_k: int = DEFAULT_RRF_K,
) -> pd.DataFrame:
    """Compute top-5/top-6 success rates for all subsets of models (af3, boltz2, boltz2affinity, rf3)."""
    all_models = ["af3", "boltz2", "boltz2affinity", "rf3"]
    model_combos: list[list[str]] = []
    for r in range(1, len(all_models) + 1):
        for combo in combinations(all_models, r):
            model_combos.append(list(combo))

    results = []
    for combo in model_combos:
        ens = make_ensemble_ranking(rankings_df, combo, aggregation=aggregation, rrf_k=rrf_k)
        true_ranks = get_true_ahl_ranking(ablation_df, ens)
        if true_ranks.empty:
            continue
        top5 = (true_ranks.true_pred_ranking <= 5).mean()
        top6 = (true_ranks.true_pred_ranking <= 6).mean()
        top3 = (true_ranks.true_pred_ranking <= 3).mean()
        top1 = (true_ranks.true_pred_ranking <= 1).mean()
        results.append(
            {
                "top 5 success rate": top5,
                "top 6 success rate": top6,
                "top 3 success rate": top3,
                "top 1 success rate": top1,
                "af3 included": "af3" in combo,
                "boltz2 included": "boltz2" in combo,
                "boltz2 affinity included": "boltz2affinity" in combo,
                "rf3 included": "rf3" in combo,
            }
        )
    if not results:
        # No labeled receptors had their cognate ligand in the scheduled set,
        # so success rates are undefined — emit an empty table with the
        # expected columns so downstream generate steps can degrade gracefully.
        return pd.DataFrame(columns=[
            "top 5 success rate", "top 6 success rate",
            "top 3 success rate", "top 1 success rate",
            "af3 included", "boltz2 included",
            "boltz2 affinity included", "rf3 included",
        ])
    return (
        pd.DataFrame(results)
        .sort_values("top 6 success rate", ascending=False)
        .reset_index(drop=True)
    )


def build_done_dict(all_df: pd.DataFrame) -> dict[str, list[str]]:
    """Build luxr_ahl -> list of completed model names for done_pairs.json."""
    done_dict: dict[str, list[str]] = {}
    if "affinity_probability_binary" in all_df.columns:
        affinity_pairs = set(
            all_df[
                (all_df.model == "boltz2affinity")
                & (all_df.affinity_probability_binary.notna())
            ]["luxr_ahl"]
            .str.upper()
            .unique()
        )
    else:
        affinity_pairs = set()

    for luxr_ahl, model in all_df[["luxr_ahl", "model"]].drop_duplicates().values:
        luxr_ahl_upper = str(luxr_ahl).upper()
        done_dict.setdefault(luxr_ahl_upper, [])
        if model == "boltz2affinity" and luxr_ahl_upper not in affinity_pairs:
            continue
        done_dict[luxr_ahl_upper].append(str(model))
    return done_dict


def load_labels_with_true(labels_path: Path) -> pd.DataFrame:
    """Load labels CSV and filter to rows with known true_ahl_id (not unknown/pending)."""
    labels = pd.read_csv(labels_path)
    if "luxr" not in labels.columns and "refahl_id" in labels.columns:
        labels["luxr"] = labels["refahl_id"].astype(str).str.upper().str.strip()
    labels_with_true = labels[
        (labels.true_ahl_id != "unknown") & (labels.true_ahl_id != "pending")
    ].copy()
    return labels_with_true


def run_postprocess(
    processed_csv: Path,
    workflow_dir: Path,
    labels_path: Path,
    skip_rankings: bool = False,
    skip_ensemble: bool = False,
    skip_ablation: bool = False,
    skip_done_json: bool = False,
    expected_ligand_count: int = 58,
    aggregation: str = DEFAULT_AGGREGATION,
    rrf_k: int = DEFAULT_RRF_K,
) -> None:
    """Read processed CSV, write rankings, ensemble, ablation table, and done JSON to workflow_dir."""
    workflow_dir.mkdir(parents=True, exist_ok=True)
    all_df = pd.read_csv(processed_csv)

    all_rankings: pd.DataFrame | None = None
    if not skip_rankings:
        all_rankings = compute_all_rankings(
            all_df, expected_ligand_count=expected_ligand_count
        )
        rankings_path = workflow_dir / RANKINGS_CSV_NAME
        all_rankings.to_csv(rankings_path, index=False)
        print(f"Wrote {rankings_path} ({len(all_rankings)} rows)")

    labels_with_true: pd.DataFrame | None = None
    if labels_path.exists():
        labels_with_true = load_labels_with_true(labels_path)
    else:
        print(
            f"Labels file not found at {labels_path}; ensemble/ablation metrics may be skipped."
        )

    if (not skip_ensemble) and all_rankings is not None:
        ensemble_rankings = make_ensemble_ranking(
            all_rankings,
            ["af3", "boltz2", "boltz2affinity", "rf3"],
            aggregation=aggregation,
            rrf_k=rrf_k,
        )
        ensemble_path = workflow_dir / ENSEMBLE_CSV_NAME
        ensemble_rankings.to_csv(ensemble_path, index=False)
        print(f"Wrote {ensemble_path} ({len(ensemble_rankings)} rows; aggregation={aggregation})")

        if (not skip_ablation) and labels_with_true is not None:
            ablation_df = labels_with_true.copy()
            ablation_df = ablation_df[ablation_df.split == "RefAHL"].copy()
            ablation_df["refahl_confidence_rating"] = ablation_df[
                "refahl_confidence_rating"
            ].fillna(4)
            ablation_df = ablation_df[
                ablation_df.refahl_confidence_rating == 1.0
            ]
            ablation_df = ablation_df[
                ~(
                    ablation_df.luxr.str.endswith("SHORT")
                    | ablation_df.luxr.str.endswith("MED")
                    | ablation_df.luxr.str.endswith("LONG")
                    | ablation_df.luxr.str.endswith("CORRECT")
                )
            ]
            ablation_table = build_ablation_table(
                all_rankings, ablation_df, aggregation=aggregation, rrf_k=rrf_k,
            )
            ablation_path = workflow_dir / "ablation_table.csv"
            ablation_table.to_csv(ablation_path, index=False)
            print(f"Wrote {ablation_path} ({len(ablation_table)} rows)")

    if not skip_done_json:
        done_dict = build_done_dict(all_df)
        done_path = workflow_dir / DONE_PAIRS_NAME
        with done_path.open("w") as fh:
            json.dump(done_dict, fh)
        print(f"Wrote {done_path} ({len(done_dict)} pairs)")


# =============================================================================
# Phase 3: Generate CSVs
# =============================================================================


def get_rdkit() -> tuple[Any, Any]:
    """Lazy import of RDKit (Chem, rdFMCS) for MCS similarity."""
    global _RDKIT
    if _RDKIT is not None:
        return _RDKIT
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
    except ImportError as exc:
        raise RuntimeError("RDKit is required for MCS outputs.") from exc
    _RDKIT = (Chem, rdFMCS)
    return _RDKIT


def get_openpyxl() -> tuple[Any, Any]:
    """Lazy import of openpyxl and PatternFill for Excel output with highlighting."""
    try:
        import openpyxl
        from openpyxl.styles import PatternFill
    except ImportError as exc:
        raise RuntimeError("openpyxl is required for Excel outputs.") from exc
    return openpyxl, PatternFill


def normalize_token(value: object) -> str:
    """Normalize a string token to stripped uppercase."""
    return str(value).strip().upper()


def normalize_optional_token(value: object) -> str | None:
    """Normalize token or return None if empty/unknown/pending."""
    if pd.isna(value):
        return None
    token = normalize_token(value)
    if token in {"", "NAN", "NONE", "UNKNOWN", "PENDING"}:
        return None
    return token


def mcs_similarity_smiles(
    smiles_a: str, smiles_b: str, timeout: int = 10
) -> float:
    """Compute MCS-based similarity between two SMILES (common atoms / (|A| + |B| - common))."""
    Chem, rdFMCS = get_rdkit()
    mol_a = Chem.MolFromSmiles(smiles_a) if isinstance(smiles_a, str) else None
    mol_b = Chem.MolFromSmiles(smiles_b) if isinstance(smiles_b, str) else None
    if mol_a is None or mol_b is None:
        return float("nan")
    if Chem.MolToSmiles(mol_a, isomericSmiles=True) == Chem.MolToSmiles(
        mol_b, isomericSmiles=True
    ):
        return 1.0
    mcs = rdFMCS.FindMCS(
        [mol_a, mol_b],
        timeout=timeout,
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
    )
    if mcs is None or mcs.numAtoms == 0:
        return 0.0
    common = float(mcs.numAtoms)
    a = float(mol_a.GetNumHeavyAtoms())
    b = float(mol_b.GetNumHeavyAtoms())
    denom = a + b - common
    if denom <= 0:
        return float("nan")
    return float(common / denom)


def compute_avg_agreement_all(
    group: pd.DataFrame, rank_cols: list[str]
) -> float:
    """Average Kendall tau over all pairs of rank columns in group."""
    n = len(rank_cols)
    if n < 2:
        return float("nan")
    agreements: list[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            tau, _ = kendalltau(group[rank_cols[i]], group[rank_cols[j]])
            agreements.append(float(tau))
    return float(np.mean(agreements)) if agreements else float("nan")


def load_inputs(
    input_dir: Path,
    ligands_path: Path,
    ablation_path: Path,
    labels_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load labels, rankings, ensemble, ligands, ablation from paths."""
    labels = pd.read_csv(labels_path)
    rankings = pd.read_csv(input_dir / RANKINGS_CSV_NAME)
    ensemble = pd.read_csv(input_dir / ENSEMBLE_CSV_NAME)
    ligands = pd.read_csv(ligands_path)
    print("Loading ablation table from", ablation_path)
    ablation = pd.read_csv(ablation_path)
    return labels, rankings, ensemble, ligands, ablation


def normalize_inputs(
    labels: pd.DataFrame,
    rankings: pd.DataFrame,
    ensemble: pd.DataFrame,
    ligands: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Normalize and standardize columns for labels, rankings, ensemble, ligands."""
    labels_norm = labels.copy()
    labels_norm["LuxR_name"] = labels_norm["refahl_id"].apply(normalize_token)
    labels_norm["split"] = labels_norm["split"].astype(str)
    labels_norm["confidence"] = pd.to_numeric(
        labels_norm["refahl_confidence_rating"], errors="coerce"
    )
    labels_norm["true_ahl_id"] = labels_norm["true_ahl_id"].apply(
        normalize_optional_token
    )
    labels_norm = labels_norm[
        ["LuxR_name", "split", "confidence", "true_ahl_id", "true_ahl_abbrev"]
    ].drop_duplicates()

    rankings_norm = rankings.copy()
    rankings_norm["luxr"] = rankings_norm["luxr"].apply(normalize_token)
    rankings_norm["ahl"] = rankings_norm["ahl"].apply(normalize_token)
    rankings_norm["model"] = (
        rankings_norm["model"].astype(str).str.strip().str.lower()
    )
    rankings_norm["rank"] = pd.to_numeric(rankings_norm["rank"], errors="coerce")
    rankings_norm["luxr_ahl"] = rankings_norm["luxr"] + "_" + rankings_norm["ahl"]

    ensemble_norm = ensemble.copy()
    ensemble_norm["luxr_ahl"] = ensemble_norm["luxr_ahl"].apply(normalize_token)
    split = ensemble_norm["luxr_ahl"].str.rsplit("_", n=1)
    ensemble_norm["luxr"] = split.str[0]
    ensemble_norm["ahl"] = split.str[1]
    ensemble_norm["ensemble_rank"] = pd.to_numeric(
        ensemble_norm["ensemble_rank"], errors="coerce"
    )

    ligands_norm = ligands.copy()
    ligands_norm["AHLid"] = ligands_norm["AHLid"].apply(normalize_token)
    ligands_norm["abbreviation"] = ligands_norm["abbreviation"].where(
        ligands_norm["abbreviation"].notna(), np.nan
    )
    ligands_norm["subcategory"] = ligands_norm["subcategory"].where(
        ligands_norm["subcategory"].notna(), np.nan
    )
    ligands_norm[SMILES_COL] = ligands_norm[SMILES_COL].where(
        ligands_norm[SMILES_COL].notna(), np.nan
    )

    return labels_norm, rankings_norm, ensemble_norm, ligands_norm


def apply_global_suffix_exclusion(labels_df: pd.DataFrame) -> pd.DataFrame:
    """Drop label rows whose LuxR_name ends with EXCLUDE_SUFFIXES."""
    return labels_df[
        ~labels_df["LuxR_name"].str.endswith(EXCLUDE_SUFFIXES, na=False)
    ].copy()


def build_true_rank_lookup(
    labels_df: pd.DataFrame,
    rankings_df: pd.DataFrame,
    ensemble_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build LuxR_name x model_name -> true_rank (rank of true AHL) lookup."""
    valid = labels_df[labels_df["true_ahl_id"].notna()][
        ["LuxR_name", "true_ahl_id"]
    ].drop_duplicates().copy()
    valid["luxr_ahl"] = valid["LuxR_name"] + "_" + valid["true_ahl_id"]

    single_true = rankings_df.merge(
        valid[["LuxR_name", "luxr_ahl"]], on="luxr_ahl", how="inner"
    )[["LuxR_name", "model", "rank"]].copy()
    single_true["model_name"] = single_true["model"].map(MODEL_NAME_MAP)
    single_true = single_true[single_true["model_name"].notna()][
        ["LuxR_name", "model_name", "rank"]
    ]
    single_true = single_true.rename(columns={"rank": "true_rank"})

    ensemble_true = ensemble_df.merge(
        valid[["LuxR_name", "luxr_ahl"]], on="luxr_ahl", how="inner"
    )[["LuxR_name", "ensemble_rank"]].copy()
    ensemble_true["model_name"] = "Ensemble"
    ensemble_true = ensemble_true.rename(columns={"ensemble_rank": "true_rank"})

    return pd.concat([single_true, ensemble_true], ignore_index=True)


def build_ranking_performance(
    subset_labels: pd.DataFrame,
    true_rank_lookup: pd.DataFrame,
) -> pd.DataFrame:
    """One row per LuxR_name with columns AF3, Boltz2, ... Ensemble = true rank of true AHL."""
    out = subset_labels[["LuxR_name"]].drop_duplicates().copy()
    for model_name in RANKING_MODELS_ORDER:
        model_part = true_rank_lookup[
            true_rank_lookup["model_name"] == model_name
        ][["LuxR_name", "true_rank"]].copy()
        model_part = model_part.rename(columns={"true_rank": model_name})
        out = out.merge(model_part, on="LuxR_name", how="left")
    return out.sort_values("LuxR_name").reset_index(drop=True)


def build_acyl_moiety(
    subset_labels: pd.DataFrame,
    ligands_df: pd.DataFrame,
    true_rank_lookup: pd.DataFrame,
    rank_model_name: str = "Ensemble",
) -> pd.DataFrame:
    """Per-LuxR model true rank broken down by subcategory (Short, Medium, Long, Atypical)."""
    model_true = true_rank_lookup[
        true_rank_lookup["model_name"] == rank_model_name
    ][["LuxR_name", "true_rank"]].rename(columns={"true_rank": "All"})

    meta = subset_labels.merge(
        ligands_df[["AHLid", "subcategory"]],
        left_on="true_ahl_id",
        right_on="AHLid",
        how="left",
    )[["LuxR_name", "subcategory"]]

    out = subset_labels[["LuxR_name"]].drop_duplicates().merge(
        model_true, on="LuxR_name", how="left"
    )
    out = out.merge(meta, on="LuxR_name", how="left")

    for sub in ["Short", "Medium", "Long", "Atypical"]:
        out[sub] = np.where(out["subcategory"] == sub, out["All"], np.nan)
    out = out[["LuxR_name", "All", "Short", "Medium", "Long", "Atypical"]]
    return out.sort_values("LuxR_name").reset_index(drop=True)


def build_chemical_similarity(
    subset_labels: pd.DataFrame,
    ensemble_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    sim_table: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Max chemical similarity to true AHL in top-1/3/5/6 ensemble predictions per LuxR.

    If ``sim_table`` is provided (square AHLid-indexed DataFrame, e.g. produced
    by ``copal sim-matrix``), scores come from O(1) lookup. Otherwise falls
    back to on-the-fly ``mcs_similarity_smiles`` (the historical baseline
    using bondCompare=CompareOrder).
    """
    true_by_luxr = subset_labels[["LuxR_name", "true_ahl_id"]].drop_duplicates()

    top1 = ensemble_df[ensemble_df["ensemble_rank"] <= 1][["luxr", "ahl"]].copy()
    top3 = ensemble_df[ensemble_df["ensemble_rank"] <= 3][["luxr", "ahl"]].copy()
    top5 = ensemble_df[ensemble_df["ensemble_rank"] <= 5][["luxr", "ahl"]].copy()
    top6 = ensemble_df[ensemble_df["ensemble_rank"] <= 6][["luxr", "ahl"]].copy()

    if sim_table is not None:
        def sim_score(a_id: str, b_id: str) -> float:
            try:
                return float(sim_table.loc[a_id, b_id])
            except KeyError:
                return float("nan")
    else:
        ahl_to_smiles = dict(zip(ligands_df["AHLid"], ligands_df[SMILES_COL]))

        def sim_score(a_id: str, b_id: str) -> float:
            sa = ahl_to_smiles.get(a_id)
            sb = ahl_to_smiles.get(b_id)
            if not isinstance(sa, str) or not isinstance(sb, str):
                return float("nan")
            return mcs_similarity_smiles(sa, sb)

    def max_mcs_for_luxr(
        luxr: str, true_ahl: str | None, topk_df: pd.DataFrame
    ) -> float:
        if true_ahl is None:
            return float("nan")
        scores: list[float] = []
        for pred in topk_df[topk_df["luxr"] == luxr]["ahl"].tolist():
            score = sim_score(normalize_token(pred), true_ahl)
            if pd.notna(score):
                scores.append(float(score))
        return float(max(scores)) if scores else float("nan")

    recs = []
    for row in true_by_luxr.itertuples(index=False):
        recs.append(
            {
                "LuxR_name": row.LuxR_name,
                "Top 1": max_mcs_for_luxr(row.LuxR_name, row.true_ahl_id, top1),
                "Top 3": max_mcs_for_luxr(row.LuxR_name, row.true_ahl_id, top3),
                "Top 5": max_mcs_for_luxr(row.LuxR_name, row.true_ahl_id, top5),
                "Top 6": max_mcs_for_luxr(row.LuxR_name, row.true_ahl_id, top6),
            }
        )
    return pd.DataFrame(
        recs, columns=["LuxR_name", "Top 1", "Top 3", "Top 5", "Top 6"]
    ).sort_values("LuxR_name").reset_index(drop=True)


def build_mpkt(
    subset_labels: pd.DataFrame,
    rankings_df: pd.DataFrame,
    ensemble_true_rank: pd.DataFrame,
) -> pd.DataFrame:
    """Mean pairwise Kendall tau across models and ensemble true AHL rank per LuxR."""
    rank_pivot = (
        rankings_df.pivot_table(
            index=["luxr", "ahl"], columns="model", values="rank", aggfunc="first"
        )
        .reset_index()
        .copy()
    )
    model_rank_cols = [m for m in MODEL_NAME_MAP if m in rank_pivot.columns]
    recs = []
    for luxr, group in rank_pivot.groupby("luxr"):
        mpkt = compute_avg_agreement_all(group, model_rank_cols)
        recs.append({"LuxR_name": luxr, "mpKT": mpkt})
    mpkt_df = pd.DataFrame(recs)
    mpkt_df = subset_labels[["LuxR_name"]].drop_duplicates().merge(
        mpkt_df, on="LuxR_name", how="left"
    )
    out = mpkt_df.merge(ensemble_true_rank, on="LuxR_name", how="left")
    out = out.rename(columns={"true_rank": "Ensemble True AHL Ranking"})
    return out[
        ["LuxR_name", "mpKT", "Ensemble True AHL Ranking"]
    ].sort_values("LuxR_name").reset_index(drop=True)


def build_full_rankings(
    labels_df: pd.DataFrame,
    rankings_df: pd.DataFrame,
    ensemble_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
) -> pd.DataFrame:
    """Full pivot: LuxR x model x AHL rank, with metadata (split, confidence, true AHL info)."""
    ligand_order = ligands_df["AHLid"].dropna().astype(str).tolist()

    single = rankings_df[["luxr", "model", "ahl", "rank"]].copy()
    single["model_name"] = single["model"].map(MODEL_NAME_MAP)
    single = single[single["model_name"].notna()][
        ["luxr", "model_name", "ahl", "rank"]
    ]

    ensemble_part = ensemble_df[["luxr", "ahl", "ensemble_rank"]].copy()
    ensemble_part["model_name"] = "Ensemble"
    ensemble_part = ensemble_part.rename(columns={"ensemble_rank": "rank"})[
        ["luxr", "model_name", "ahl", "rank"]
    ]

    stacked = pd.concat([single, ensemble_part], ignore_index=True)
    pivot = (
        stacked.pivot_table(
            index=["luxr", "model_name"], columns="ahl", values="rank", aggfunc="first"
        )
        .reset_index()
        .rename(columns={"luxr": "LuxR_name"})
    )

    meta = labels_df.merge(
        ligands_df[["AHLid", "abbreviation", "subcategory"]],
        left_on="true_ahl_id",
        right_on="AHLid",
        how="left",
    ).rename(
        columns={
            "abbreviation": "true_ahl_abbrev",
            "subcategory": "true_ahl_subcategory",
        }
    )
    meta = meta[
        [
            "LuxR_name",
            "split",
            "confidence",
            "true_ahl_id",
            "true_ahl_abbrev",
            "true_ahl_subcategory",
        ]
    ].drop_duplicates()

    out = meta.merge(pivot, on="LuxR_name", how="inner")
    for lig in ligand_order:
        if lig not in out.columns:
            out[lig] = np.nan

    model_order_map = {name: i for i, name in enumerate(RANKING_MODELS_ORDER)}
    out["__model_order"] = out["model_name"].map(model_order_map).fillna(999)
    out = out.sort_values(
        ["LuxR_name", "split", "confidence", "__model_order"]
    ).drop(columns="__model_order")

    ordered_cols = [
        "LuxR_name",
        "split",
        "confidence",
        "model_name",
        "true_ahl_id",
        "true_ahl_abbrev",
        "true_ahl_subcategory",
    ] + ligand_order
    return out[ordered_cols].reset_index(drop=True)


def build_ranked_excel(
    model_key: str | None,
    rankings_df: pd.DataFrame,
    ensemble_df: pd.DataFrame,
    labels_df: pd.DataFrame,
    ligands_df: pd.DataFrame,
    mpkt_df: pd.DataFrame,
    out_path: Path,
) -> None:
    """
    Write a per-model ranked Excel file with 3 rows per LuxR:
      Row 1: rank_K = AHLid at rank K
      Row 2: rank_K = abbreviation of AHLid at rank K
      Row 3: rank_K = score of AHLid at rank K
    Cells in the rank_K column matching the true AHL are highlighted yellow.

    model_key: lowercase model name ("af3", "boltz2", etc.) or None for Ensemble.
    """
    import math
    openpyxl, PatternFill = get_openpyxl()
    yellow_fill = PatternFill(start_color="FFFF00", end_color="FFFF00", fill_type="solid")
    mpkt_orange_fill = PatternFill(start_color="FFA500", end_color="FFA500", fill_type="solid")
    mpkt_green_fill = PatternFill(start_color="92D050", end_color="92D050", fill_type="solid")

    ahl_to_abbrev = dict(zip(ligands_df["AHLid"], ligands_df["abbreviation"]))

    # Build per-LuxR sorted AHL list: [(ahl_id, score), ...]
    if model_key is None:
        # Ensemble
        score_col = "ensemble_score"
        rank_col = "ensemble_rank"
        model_df = ensemble_df.copy()
        model_df["luxr"] = model_df["luxr"].apply(normalize_token)
        model_df["ahl"] = model_df["ahl"].apply(normalize_token)
    else:
        score_col = "iptm" if model_key != "boltz2affinity" else "affinity_probability_binary"
        rank_col = "rank"
        model_df = rankings_df[rankings_df["model"] == model_key].copy()

    # mpKT lookup: LuxR_name -> mpKT
    mpkt_lookup = dict(zip(mpkt_df["LuxR_name"], mpkt_df["mpKT"]))

    # Labels: one row per LuxR (already de-duped upstream but ensure here)
    labels_meta = labels_df[
        ["LuxR_name", "split", "confidence", "true_ahl_id", "true_ahl_abbrev"]
    ].drop_duplicates(subset=["LuxR_name"]).copy()

    # Merge true_ahl_subcategory from ligands
    labels_meta = labels_meta.merge(
        ligands_df[["AHLid", "subcategory"]],
        left_on="true_ahl_id",
        right_on="AHLid",
        how="left",
    ).rename(columns={"subcategory": "true_ahl_subcategory"})
    labels_meta = labels_meta.sort_values("LuxR_name").reset_index(drop=True)

    n_ranks = 58
    meta_cols = [
        "LuxR_name", "split", "confidence", "true_ahl_id",
        "true_ahl_abbrev", "true_ahl_subcategory", "mpKT",
        "true_ahl_id_ranking",
    ]
    rank_cols = [f"rank_{k}" for k in range(1, n_ranks + 1)]
    all_cols = meta_cols + rank_cols

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(all_cols)

    for _, lrow in labels_meta.iterrows():
        luxr = lrow["LuxR_name"]
        true_ahl = lrow["true_ahl_id"]
        mpkt_val = mpkt_lookup.get(luxr, float("nan"))
        # Get sorted AHLs for this LuxR
        ldf = model_df[model_df["luxr"] == luxr].sort_values(rank_col)
        sorted_ahls = ldf["ahl"].tolist()
        sorted_scores = ldf[score_col].tolist()

        # Pad or trim to n_ranks
        sorted_ahls = sorted_ahls[:n_ranks]
        sorted_scores = sorted_scores[:n_ranks]
        while len(sorted_ahls) < n_ranks:
            sorted_ahls.append(None)
            sorted_scores.append(None)

        # Find which rank position contains the true AHL (1-indexed)
        true_rank_pos = None
        if true_ahl is not None:
            try:
                true_rank_pos = sorted_ahls.index(true_ahl) + 1  # 1-indexed
            except ValueError:
                true_rank_pos = None

        meta_vals = [
            luxr,
            lrow["split"],
            lrow["confidence"],
            true_ahl,
            lrow["true_ahl_abbrev"],
            lrow.get("true_ahl_subcategory", None),
            mpkt_val,
            true_rank_pos,
        ]

        # Row 1: AHLids
        row1 = meta_vals + sorted_ahls
        # Row 2: abbreviations
        row2 = meta_vals + [
            ahl_to_abbrev.get(a, "") if a is not None else None for a in sorted_ahls
        ]
        # Row 3: scores
        row3 = meta_vals + sorted_scores

        for row_vals in [row1, row2, row3]:
            ws.append(row_vals)

        current_row = ws.max_row
        # Color mpKT cells for all 3 rows: orange if < 0.3, green if >= 0.3
        mpkt_col_idx = meta_cols.index("mpKT") + 1
        if mpkt_val is not None and not (isinstance(mpkt_val, float) and math.isnan(mpkt_val)):
            fill = mpkt_orange_fill if mpkt_val < 0.3 else mpkt_green_fill
            for r in [current_row - 2, current_row - 1, current_row]:
                ws.cell(row=r, column=mpkt_col_idx).fill = fill

        if true_rank_pos is not None:
            # Column index in sheet: len(meta_cols) + true_rank_pos (1-indexed)
            col_idx = len(meta_cols) + true_rank_pos
            for r in [current_row - 2, current_row - 1, current_row]:
                ws.cell(row=r, column=col_idx).fill = yellow_fill

        # Empty separator row between LuxRs
        ws.append([])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)
    print(f"Wrote {out_path}")


def build_mpkt_correlation_test(mpkt_df: pd.DataFrame, out_path: Path) -> None:
    """Compute Spearman, Pearson, and Kendall correlations between mpKT and ensemble ranking."""
    df = mpkt_df.dropna(subset=["mpKT", "Ensemble True AHL Ranking"]).copy()
    if len(df) < 3:
        return
    x = df["mpKT"].values
    y = df["Ensemble True AHL Ranking"].values
    sp_r, sp_p = spearmanr(x, y)
    pe_r, pe_p = pearsonr(x, y)
    kt_tau, kt_p = kendalltau(x, y)
    results = pd.DataFrame([
        {"test": "spearman", "statistic": sp_r, "pvalue": sp_p, "n": len(df)},
        {"test": "pearson", "statistic": pe_r, "pvalue": pe_p, "n": len(df)},
        {"test": "kendall_tau", "statistic": kt_tau, "pvalue": kt_p, "n": len(df)},
    ])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


def build_mpkt_quartile_top6_success(mpkt_df: pd.DataFrame) -> pd.DataFrame:
    """Per mpKT quartile: n, top3/top5/top6 success rates."""
    cols = [
        "mpkt_quartile", "quartile_min", "quartile_max", "quartile_center", "n",
        "top1_success_rate", "top3_success_rate", "top5_success_rate", "top6_success_rate",
    ]
    df = mpkt_df.dropna(subset=["mpKT", "Ensemble True AHL Ranking"]).copy()
    if df.empty:
        return pd.DataFrame(columns=cols)
    labels = ["Q1", "Q2", "Q3", "Q4"]
    try:
        df["mpkt_quartile"] = pd.qcut(
            df["mpKT"], q=4, labels=labels, duplicates="drop"
        )
    except ValueError:
        return pd.DataFrame(columns=cols)

    recs = []
    for quartile in df["mpkt_quartile"].dropna().unique():
        qdf = df[df["mpkt_quartile"] == quartile]
        recs.append(
            {
                "mpkt_quartile": str(quartile),
                "quartile_min": qdf["mpKT"].min(),
                "quartile_max": qdf["mpKT"].max(),
                "quartile_center": (qdf["mpKT"].min() + qdf["mpKT"].max()) / 2,
                "n": int(len(qdf)),
                "top1_success_rate": float(
                    (qdf["Ensemble True AHL Ranking"] <= 1).mean()
                ),
                "top3_success_rate": float(
                    (qdf["Ensemble True AHL Ranking"] <= 3).mean()
                ),
                "top5_success_rate": float(
                    (qdf["Ensemble True AHL Ranking"] <= 5).mean()
                ),
                "top6_success_rate": float(
                    (qdf["Ensemble True AHL Ranking"] <= 6).mean()
                ),
            }
        )
    return pd.DataFrame(recs).sort_values("mpkt_quartile").reset_index(drop=True)


SCORE_COL_BY_MODEL: dict[str, str] = {
    "AF3": "iptm",
    "Boltz2": "iptm",
    "Boltz2Affinity": "affinity_probability_binary",
    "RF3": "iptm",
}


def build_score_vs_true_rank(
    subset_labels: pd.DataFrame,
    rankings_df: pd.DataFrame,
    model_name: str,
) -> pd.DataFrame:
    """Per-LuxR (mean model score over the top-6 AHLs, rank of true AHL) for one model.

    X column is the average of the model's confidence over its top-6 ranked AHL
    predictions for the LuxR — a label-free "peak confidence" signal,
    computable at inference. Column name is ``mean_top6_iptm`` for
    AF3/Boltz2/RF3 and ``mean_top6_affinity_probability_binary`` for
    Boltz2Affinity. Y column ``True AHL Ranking`` is the position of the
    cognate AHL in this model's ranking (NaN when no ``true_ahl_id``).
    """
    score_col = SCORE_COL_BY_MODEL[model_name]
    mean_col = f"mean_top6_{score_col}"
    inv = {v: k for k, v in MODEL_NAME_MAP.items()}
    model_key = inv[model_name]

    model_rows = rankings_df[rankings_df["model"] == model_key]

    means = (
        model_rows[model_rows["rank"] <= 6]
        .groupby("luxr")[score_col]
        .mean()
        .rename(mean_col)
        .reset_index()
        .rename(columns={"luxr": "LuxR_name"})
    )

    valid = subset_labels[subset_labels["true_ahl_id"].notna()][
        ["LuxR_name", "true_ahl_id"]
    ].drop_duplicates().copy()
    valid["luxr_ahl"] = valid["LuxR_name"] + "_" + valid["true_ahl_id"]
    true_rank = (
        model_rows[["luxr_ahl", "rank"]]
        .merge(valid[["LuxR_name", "luxr_ahl"]], on="luxr_ahl", how="inner")
        [["LuxR_name", "rank"]]
        .rename(columns={"rank": "True AHL Ranking"})
    )

    all_luxrs = subset_labels[["LuxR_name"]].drop_duplicates()
    out = (
        all_luxrs
        .merge(means, on="LuxR_name", how="left")
        .merge(true_rank, on="LuxR_name", how="left")
    )
    return out[["LuxR_name", mean_col, "True AHL Ranking"]].sort_values(
        "LuxR_name"
    ).reset_index(drop=True)


def build_score_vs_rank_correlation(
    tables: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Per-model Spearman / Pearson / Kendall between mean-top6 score and True AHL Ranking."""
    recs = []
    for model_name, df in tables.items():
        mean_col = f"mean_top6_{SCORE_COL_BY_MODEL[model_name]}"
        sub = df.dropna(subset=[mean_col, "True AHL Ranking"])
        n = len(sub)
        rec = {
            "model": model_name, "score_col": mean_col, "n": n,
            "spearman_r": float("nan"), "spearman_p": float("nan"),
            "pearson_r": float("nan"), "pearson_p": float("nan"),
            "kendall_tau": float("nan"), "kendall_p": float("nan"),
        }
        if n >= 3:
            x = sub[mean_col].astype(float).values
            y = sub["True AHL Ranking"].astype(float).values
            sp_r, sp_p = spearmanr(x, y)
            pe_r, pe_p = pearsonr(x, y)
            kt_tau, kt_p = kendalltau(x, y)
            rec.update({
                "spearman_r": float(sp_r), "spearman_p": float(sp_p),
                "pearson_r": float(pe_r), "pearson_p": float(pe_p),
                "kendall_tau": float(kt_tau), "kendall_p": float(kt_p),
            })
        recs.append(rec)
    return pd.DataFrame(recs)


def build_ablation_summary(ablation_table: pd.DataFrame) -> pd.DataFrame:
    """Summarize ablation by num_models: mean/std of top-5 and top-6 success rate, N."""
    df = ablation_table.copy()
    if df.empty:
        # No labeled receptors had their cognate in the scheduled ligand set;
        # `build_ablation_table` emitted an empty table. Carry the emptiness
        # forward so downstream writers don't hit KeyError on the missing
        # columns.
        return pd.DataFrame()
    df["num_models"] = (
        df["af3 included"].astype(int)
        + df["boltz2 included"].astype(int)
        + df["boltz2 affinity included"].astype(int)
        + df["rf3 included"].astype(int)
    )
    recs = []
    for num_models, mdf in df.groupby("num_models"):
        recs.append(
            {
                "num_models": num_models,
                "top 6 success rate": mdf["top 6 success rate"].mean(),
                "top 6 success rate std": mdf["top 6 success rate"].std(),
                "top 6 success rate min": mdf["top 6 success rate"].min(),
                "top 6 success rate max": mdf["top 6 success rate"].max(),
                "top 5 success rate": mdf["top 5 success rate"].mean(),
                "top 5 success rate std": mdf["top 5 success rate"].std(),
                "top 5 success rate min": mdf["top 5 success rate"].min(),
                "top 5 success rate max": mdf["top 5 success rate"].max(),
                "top 3 success rate":     mdf["top 3 success rate"].mean(),
                "top 3 success rate std": mdf["top 3 success rate"].std(),
                "top 3 success rate min": mdf["top 3 success rate"].min(),
                "top 3 success rate max": mdf["top 3 success rate"].max(),
                "top 1 success rate":     mdf["top 1 success rate"].mean(),
                "top 1 success rate std": mdf["top 1 success rate"].std(),
                "top 1 success rate min": mdf["top 1 success rate"].min(),
                "top 1 success rate max": mdf["top 1 success rate"].max(),
                "N": int(mdf.shape[0]),
            }
        )
    return pd.DataFrame(recs).sort_values("num_models").reset_index(drop=True)


def write_results_csv(df: pd.DataFrame, out_path: Path) -> None:
    """Write df to out_path, creating parent dirs as needed."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)


def _normalize_rankings(rankings_raw: pd.DataFrame) -> pd.DataFrame:
    """Apply the rankings-half of `normalize_inputs`."""
    out = rankings_raw.copy()
    out["luxr"] = out["luxr"].apply(normalize_token)
    out["ahl"] = out["ahl"].apply(normalize_token)
    out["model"] = out["model"].astype(str).str.strip().str.lower()
    out["rank"] = pd.to_numeric(out["rank"], errors="coerce")
    out["luxr_ahl"] = out["luxr"] + "_" + out["ahl"]
    return out


def _normalize_ensemble(ensemble_raw: pd.DataFrame) -> pd.DataFrame:
    """Apply the ensemble-half of `normalize_inputs`."""
    out = ensemble_raw.copy()
    out["luxr_ahl"] = out["luxr_ahl"].apply(normalize_token)
    split = out["luxr_ahl"].str.rsplit("_", n=1)
    out["luxr"] = split.str[0]
    out["ahl"] = split.str[1]
    out["ensemble_rank"] = pd.to_numeric(out["ensemble_rank"], errors="coerce")
    return out


def _derive_ensemble(
    rankings: pd.DataFrame,
    aggregation: str,
    rrf_k: int,
) -> pd.DataFrame:
    """Re-derive the ensemble ranking from per-model rankings with a custom aggregation."""
    models = sorted(set(rankings["model"]).intersection(MODEL_NAME_MAP.keys()))
    if not models:
        raise ValueError("rankings contains no recognized model names; cannot derive ensemble.")
    ensemble = make_ensemble_ranking(rankings, models, aggregation=aggregation, rrf_k=rrf_k)
    return _normalize_ensemble(ensemble)


def run_analysis_for_set(
    bset,
    workflow_dir: Path,
    output_root: Path,
    aggregation: str = DEFAULT_AGGREGATION,
    rrf_k: int = DEFAULT_RRF_K,
    per_model_split: tuple[str, ...] = ("AF3", "Boltz2", "Boltz2Affinity", "RF3"),
    skip_ablation: bool = False,
    chem_sim_matrix_path: Path | None = None,
) -> None:
    """Produce every analysis artifact for ONE resolved benchmark set.

    `bset` is a `copal.analysis.sets.BenchmarkSet` — already filtered by split,
    confidence, leakage, and suffix. No internal looping over conf or leakage.
    Outputs land under ``output_root / bset.name / <artifact_kind> / <file>.csv``.
    """
    import logging

    log = logging.getLogger(__name__)

    if bset.is_empty():
        log.error("Benchmark set %r is empty; nothing to write.", bset.name)
        raise SystemExit(2)

    workflow_dir = Path(workflow_dir)
    output_root = Path(output_root)

    rankings_csv = workflow_dir / RANKINGS_CSV_NAME
    ensemble_csv = workflow_dir / ENSEMBLE_CSV_NAME
    if not rankings_csv.exists():
        raise FileNotFoundError(f"Required rankings CSV missing: {rankings_csv}")
    if not ensemble_csv.exists():
        raise FileNotFoundError(f"Required ensemble CSV missing: {ensemble_csv}")

    rankings = _normalize_rankings(pd.read_csv(rankings_csv))
    ensemble = _normalize_ensemble(pd.read_csv(ensemble_csv))

    if aggregation != DEFAULT_AGGREGATION:
        log.info(
            "Aggregation %r differs from cached default %r; re-deriving ensemble in memory.",
            aggregation, DEFAULT_AGGREGATION,
        )
        ensemble = _derive_ensemble(rankings, aggregation, rrf_k)

    base = output_root / bset.name
    dirs = {
        "ranking_performance": base / "ranking_performance",
        "acyl_moiety": base / "acyl_moiety",
        "chemical_similarity": base / "chemical_similarity",
        "mpkt": base / "mpkt",
        "score_vs_rank": base / "score_vs_rank",
        "full_rankings": base / "full_rankings",
        "ablation": base / "ablation",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    labels = bset.labels
    ligands = bset.ligands

    lookup = build_true_rank_lookup(labels, rankings, ensemble)

    write_results_csv(
        build_ranking_performance(labels, lookup),
        dirs["ranking_performance"] / "ranking_performance.csv",
    )

    write_results_csv(
        build_acyl_moiety(labels, ligands, lookup),
        dirs["acyl_moiety"] / "acyl_moiety.csv",
    )
    for model_name in per_model_split:
        write_results_csv(
            build_acyl_moiety(labels, ligands, lookup, rank_model_name=model_name),
            dirs["acyl_moiety"] / f"acyl_moiety_{model_name}.csv",
        )

    sim_table = None
    if chem_sim_matrix_path is not None:
        sim_path = Path(chem_sim_matrix_path)
        if not sim_path.exists():
            raise FileNotFoundError(
                f"chem_sim_matrix path does not exist: {sim_path}. "
                "Generate one with `copal sim-matrix method=<m>`."
            )
        sim_table = pd.read_csv(sim_path, index_col=0)
        log.info("Loaded chem-sim matrix %s (shape=%s)", sim_path, sim_table.shape)

    write_results_csv(
        build_chemical_similarity(labels, ensemble, ligands, sim_table=sim_table),
        dirs["chemical_similarity"] / "chemical_similarity.csv",
    )

    ensemble_true = (
        lookup[lookup["model_name"] == "Ensemble"][["LuxR_name", "true_rank"]]
        .drop_duplicates()
    )
    mpkt_df = build_mpkt(labels, rankings, ensemble_true)
    write_results_csv(mpkt_df, dirs["mpkt"] / "mpkt.csv")

    valid_rows = mpkt_df.dropna(subset=["mpKT", "Ensemble True AHL Ranking"])
    if len(valid_rows) >= 4:
        quartile_df = build_mpkt_quartile_top6_success(mpkt_df)
        quartile_df.to_csv(dirs["mpkt"] / "mpkt_quartile_top6.csv", index=False)
    else:
        log.warning(
            "Set %r has %d rows with both mpKT and ensemble rank; skipping mpkt_quartile_top6.csv.",
            bset.name, len(valid_rows),
        )
    if len(valid_rows) >= 3:
        build_mpkt_correlation_test(mpkt_df, dirs["mpkt"] / "mpkt_correlation.csv")
    else:
        log.warning(
            "Set %r has %d rows for correlation; skipping mpkt_correlation.csv.",
            bset.name, len(valid_rows),
        )

    sr_tables: dict[str, pd.DataFrame] = {}
    for model_name in per_model_split:
        if model_name not in SCORE_COL_BY_MODEL:
            continue
        df = build_score_vs_true_rank(labels, rankings, model_name)
        write_results_csv(df, dirs["score_vs_rank"] / f"{model_name}.csv")
        sr_tables[model_name] = df
    if sr_tables:
        corr = build_score_vs_rank_correlation(sr_tables)
        corr.to_csv(dirs["score_vs_rank"] / "correlation.csv", index=False)

    write_results_csv(
        build_full_rankings(labels, rankings, ensemble, ligands),
        dirs["full_rankings"] / "full_rankings.csv",
    )

    available_models = sorted(set(rankings["model"]).intersection(MODEL_NAME_MAP.keys()))
    if skip_ablation:
        log.info("skip_ablation set; skipping ablation for set %r.", bset.name)
    elif len(available_models) < 2:
        log.warning(
            "Set %r has only %d model(s) in rankings (%r); skipping ablation.",
            bset.name, len(available_models), available_models,
        )
    else:
        abl_pool = labels[labels["true_ahl_id"].notna()].copy()
        abl_pool["luxr"] = abl_pool["LuxR_name"]
        abl_table = build_ablation_table(
            rankings, abl_pool, aggregation=aggregation, rrf_k=rrf_k,
        )
        abl_table.to_csv(dirs["ablation"] / "ablation_table_raw.csv", index=False)
        abl_summary = build_ablation_summary(abl_table)
        abl_summary.to_csv(dirs["ablation"] / "ablation_table.csv", index=False)

    log.info(
        "run_analysis_for_set complete: set=%r n_luxrs=%d out=%s",
        bset.name, bset.n_luxrs(), base,
    )


def run_full_ranking_for_set(
    bset,
    workflow_dir: Path,
    output_root: Path,
    aggregation: str = DEFAULT_AGGREGATION,
    rrf_k: int = DEFAULT_RRF_K,
    per_model_split: tuple[str, ...] = ("AF3", "Boltz2", "Boltz2Affinity", "RF3"),
) -> None:
    """Emit stylized per-model + ensemble ranking .xlsx files (and full_rankings.csv).

    Reads the canonical per_model/ensemble CSVs from ``workflow_dir`` and writes:
        output_root/<bset.name>/full_rankings/full_rankings.csv
        output_root/<bset.name>/full_rankings/{AF3,Boltz2,Boltz2Affinity,RF3,ENSEMBLE}.xlsx

    Works in benchmark mode (bset.labels / bset.ligands populated) and inference
    mode (both empty — synthesized from rankings; yellow-fill on true AHL is
    naturally skipped, mpKT-based orange/green coloring still renders).
    """
    import logging
    log = logging.getLogger(__name__)

    workflow_dir = Path(workflow_dir)
    output_root = Path(output_root)
    rankings_csv = workflow_dir / RANKINGS_CSV_NAME
    ensemble_csv = workflow_dir / ENSEMBLE_CSV_NAME
    if not rankings_csv.exists():
        raise FileNotFoundError(f"Required rankings CSV missing: {rankings_csv}")
    if not ensemble_csv.exists():
        raise FileNotFoundError(f"Required ensemble CSV missing: {ensemble_csv}")

    rankings = _normalize_rankings(pd.read_csv(rankings_csv))
    ensemble = _normalize_ensemble(pd.read_csv(ensemble_csv))
    if aggregation != DEFAULT_AGGREGATION:
        log.info(
            "Aggregation %r differs from cached default %r; re-deriving ensemble in memory.",
            aggregation, DEFAULT_AGGREGATION,
        )
        ensemble = _derive_ensemble(rankings, aggregation, rrf_k)

    if bset.labels.empty:
        log.info("Inference mode: synthesizing labels from rankings (no true_ahl_id).")
        labels = pd.DataFrame({
            "LuxR_name": sorted(rankings["luxr"].astype(str).unique()),
        })
        labels["split"] = pd.NA
        labels["confidence"] = pd.NA
        labels["true_ahl_id"] = pd.NA
        labels["true_ahl_abbrev"] = pd.NA
    else:
        labels = bset.labels

    if bset.ligands.empty:
        log.info("Inference mode: synthesizing ligands from rankings (abbreviation=AHLid).")
        ahls = sorted(rankings["ahl"].astype(str).unique())
        ligands = pd.DataFrame({
            "AHLid": ahls,
            "abbreviation": ahls,
            "subcategory": pd.NA,
        })
    else:
        ligands = bset.ligands

    if labels.empty:
        log.error("Set %r is empty after loading; nothing to write.", bset.name)
        raise SystemExit(2)

    full_dir = output_root / bset.name / "full_rankings"
    full_dir.mkdir(parents=True, exist_ok=True)

    # mpKT per LuxR (drives orange/green fill in the .xlsx). In inference mode
    # ensemble_true is empty -> "Ensemble True AHL Ranking" column is all NaN.
    lookup = build_true_rank_lookup(labels, rankings, ensemble)
    ensemble_true = (
        lookup[lookup["model_name"] == "Ensemble"][["LuxR_name", "true_rank"]]
        .drop_duplicates()
    )
    mpkt_df = build_mpkt(labels, rankings, ensemble_true)

    write_results_csv(
        build_full_rankings(labels, rankings, ensemble, ligands),
        full_dir / "full_rankings.csv",
    )

    excel_models: list[tuple[str, str | None]] = [
        *((m, m.lower()) for m in per_model_split),
        ("ENSEMBLE", None),
    ]
    for display_name, model_key in excel_models:
        build_ranked_excel(
            model_key, rankings, ensemble, labels, ligands, mpkt_df,
            full_dir / f"{display_name}.xlsx",
        )

    log.info(
        "run_full_ranking_for_set complete: set=%r n_luxrs=%d out=%s",
        bset.name, int(labels["LuxR_name"].nunique()), full_dir,
    )


def run_generate(
    workflow_dir: Path,
    output_root: Path,
    labels_path: Path,
    ligands_path: Path,
    ablation_path: Path | None = None,
    no_leakage_path: Path | None = None,
    aggregation: str = DEFAULT_AGGREGATION,
    rrf_k: int = DEFAULT_RRF_K,
) -> None:
    """
    Read rankings, ensemble, ablation, labels, ligands from workflow_dir (and ablation_path);
    write analysis CSVs under output_root.
    """
    if ablation_path is None:
        ablation_path = workflow_dir / "ablation_table.csv"

    no_leakage_names: set[str] | None = None
    if no_leakage_path is not None and no_leakage_path.exists():
        no_leakage_names = load_no_leakage_names(no_leakage_path)
        print(f"Loaded {len(no_leakage_names)} no-leakage LuxR names from {no_leakage_path}")

    labels_raw, rankings_raw, ensemble_raw, ligands_raw, ablation_raw = load_inputs(
        workflow_dir,
        ligands_path,
        ablation_path,
        labels_path,
    )
    labels, rankings, ensemble, ligands = normalize_inputs(
        labels_raw, rankings_raw, ensemble_raw, ligands_raw
    )
    labels = apply_global_suffix_exclusion(labels)

    ranking_dir = output_root / "ranking_performance"
    acyl_dir = output_root / "acyl_moiety"
    mcs_dir = output_root / "chemical_similarity"
    mpkt_dir = output_root / "mpkt"
    full_dir = output_root / "full_rankings"
    ablation_dir = output_root / "ablation"
    for d in [ranking_dir, acyl_dir, mcs_dir, mpkt_dir, full_dir, ablation_dir]:
        d.mkdir(parents=True, exist_ok=True)

    for conf in [1, 2, 3]:
        subset = labels[
            (labels["split"].str.upper() == "REFAHL") & (labels["confidence"] == conf)
        ].copy()
        subset = subset.drop_duplicates(
            subset=["LuxR_name", "split", "confidence"]
        )
        true_rank_lookup = build_true_rank_lookup(subset, rankings, ensemble)

        perf_df = build_ranking_performance(subset, true_rank_lookup)
        perf_path = ranking_dir / f"ranking_performance_conf{conf}.csv"
        write_results_csv(perf_df, perf_path)

        chem_df = build_acyl_moiety(subset, ligands, true_rank_lookup)
        chem_path = acyl_dir / f"acyl_moiety_conf{conf}.csv"
        write_results_csv(chem_df, chem_path)
        if conf == 1:
            for model_name in ["AF3", "Boltz2", "Boltz2Affinity", "RF3"]:
                model_chem_df = build_acyl_moiety(
                    subset,
                    ligands,
                    true_rank_lookup,
                    rank_model_name=model_name,
                )
                model_chem_path = acyl_dir / f"acyl_moiety_conf1_{model_name}.csv"
                write_results_csv(model_chem_df, model_chem_path)

        mcs_df = build_chemical_similarity(subset, ensemble, ligands)
        mcs_path = mcs_dir / f"chemical_similarity_conf{conf}.csv"
        write_results_csv(mcs_df, mcs_path)

        # No-leakage subset for conf=1
        if conf == 1 and no_leakage_names is not None:
            nl_subset = subset[subset["LuxR_name"].isin(no_leakage_names)].copy()
            nl_true_rank_lookup = build_true_rank_lookup(nl_subset, rankings, ensemble)

            nl_perf_df = build_ranking_performance(nl_subset, nl_true_rank_lookup)
            write_results_csv(
                nl_perf_df,
                ranking_dir / "ranking_performance_conf1_no_leakage.csv",
            )

            nl_chem_df = build_acyl_moiety(nl_subset, ligands, nl_true_rank_lookup)
            write_results_csv(
                nl_chem_df,
                acyl_dir / "acyl_moiety_conf1_no_leakage.csv",
            )
            for model_name in ["AF3", "Boltz2", "Boltz2Affinity", "RF3"]:
                nl_model_chem_df = build_acyl_moiety(
                    nl_subset, ligands, nl_true_rank_lookup,
                    rank_model_name=model_name,
                )
                write_results_csv(
                    nl_model_chem_df,
                    acyl_dir / f"acyl_moiety_conf1_{model_name}_no_leakage.csv",
                )

            nl_mcs_df = build_chemical_similarity(nl_subset, ensemble, ligands)
            write_results_csv(
                nl_mcs_df,
                mcs_dir / "chemical_similarity_conf1_no_leakage.csv",
            )
            print(f"Wrote no-leakage conf=1 outputs ({len(nl_subset)} LuxRs)")

    mpkt_subset = labels[
        (labels["split"].str.upper() == "REFAHL") & (labels["confidence"] == 1)
    ].copy()
    mpkt_subset = mpkt_subset.drop_duplicates(
        subset=["LuxR_name", "split", "confidence"]
    )
    mpkt_true_lookup = build_true_rank_lookup(mpkt_subset, rankings, ensemble)
    ensemble_true = mpkt_true_lookup[
        mpkt_true_lookup["model_name"] == "Ensemble"
    ][["LuxR_name", "true_rank"]].drop_duplicates()
    mpkt_df = build_mpkt(mpkt_subset, rankings, ensemble_true)
    mpkt_path = mpkt_dir / "mpkt.csv"
    write_results_csv(mpkt_df, mpkt_path)
    mpkt_quartile_df = build_mpkt_quartile_top6_success(mpkt_df)
    mpkt_quartile_path = mpkt_dir / "mpkt_quartile_top6.csv"
    mpkt_quartile_df.to_csv(mpkt_quartile_path, index=False)
    build_mpkt_correlation_test(mpkt_df, mpkt_dir / "mpkt_correlation.csv")

    if no_leakage_names is not None:
        nl_mpkt_subset = mpkt_subset[mpkt_subset["LuxR_name"].isin(no_leakage_names)].copy()
        nl_mpkt_true_lookup = build_true_rank_lookup(nl_mpkt_subset, rankings, ensemble)
        nl_ensemble_true = nl_mpkt_true_lookup[
            nl_mpkt_true_lookup["model_name"] == "Ensemble"
        ][["LuxR_name", "true_rank"]].drop_duplicates()
        nl_mpkt_df = build_mpkt(nl_mpkt_subset, rankings, nl_ensemble_true)
        write_results_csv(nl_mpkt_df, mpkt_dir / "mpkt_no_leakage.csv")
        nl_mpkt_quartile_df = build_mpkt_quartile_top6_success(nl_mpkt_df)
        nl_mpkt_quartile_df.to_csv(
            mpkt_dir / "mpkt_quartile_top6_no_leakage.csv", index=False
        )
        build_mpkt_correlation_test(
            nl_mpkt_df, mpkt_dir / "mpkt_correlation_no_leakage.csv"
        )
        print(f"Wrote no-leakage mpkt outputs ({len(nl_mpkt_subset)} LuxRs)")

    ablation_summary_df = build_ablation_summary(ablation_raw)
    ablation_summary_path = ablation_dir / "ablation_table.csv"
    ablation_summary_df.to_csv(ablation_summary_path, index=False)

    if no_leakage_names is not None:
        nl_labels_with_true = labels[
            labels["true_ahl_id"].notna()
            & (labels["split"].str.upper() == "REFAHL")
            & (labels["confidence"] == 1)
            & labels["LuxR_name"].isin(no_leakage_names)
        ].copy()
        nl_labels_with_true = nl_labels_with_true[
            ~nl_labels_with_true["LuxR_name"].str.endswith(("SHORT", "MED", "LONG", "CORRECT"))
        ]
        nl_abl_df = nl_labels_with_true.copy()
        nl_abl_df["luxr"] = nl_abl_df["LuxR_name"]
        nl_ablation_table = build_ablation_table(
            rankings, nl_abl_df, aggregation=aggregation, rrf_k=rrf_k,
        )
        nl_ablation_table.to_csv(ablation_dir / "ablation_table_raw_no_leakage.csv", index=False)
        nl_ablation_summary = build_ablation_summary(nl_ablation_table)
        nl_ablation_summary.to_csv(ablation_dir / "ablation_table_no_leakage.csv", index=False)
        # Sorted: most models first, within group by AF3 > Boltz2 > B2A > RF3
        sorted_cols = [
            "top 5 success rate", "top 6 success rate", "top 3 success rate",
            "top 1 success rate", "af3 included", "boltz2 included",
            "boltz2 affinity included", "rf3 included",
        ]
        nl_sorted = nl_ablation_table.copy()
        nl_sorted["_n"] = (
            nl_sorted["af3 included"].astype(int)
            + nl_sorted["boltz2 included"].astype(int)
            + nl_sorted["boltz2 affinity included"].astype(int)
            + nl_sorted["rf3 included"].astype(int)
        )
        nl_sorted = nl_sorted.sort_values(
            ["_n", "af3 included", "boltz2 included", "boltz2 affinity included", "rf3 included"],
            ascending=[False, False, False, False, False],
        ).reset_index(drop=True)
        nl_sorted[sorted_cols].to_csv(
            ablation_dir / "sorted_ablation_full_no_leakage.csv", index=False
        )
        print(f"Wrote no-leakage ablation outputs ({len(nl_labels_with_true)} LuxRs)")

    full_rankings_df = build_full_rankings(
        labels, rankings, ensemble, ligands
    )
    full_path = full_dir / "all_splits.csv"
    write_results_csv(full_rankings_df, full_path)

    # Compute mpKT for ALL LuxRs in rankings (not filtered to RefAHL conf1 subset)
    rank_pivot_all = (
        rankings.pivot_table(
            index=["luxr", "ahl"], columns="model", values="rank", aggfunc="first"
        )
        .reset_index()
    )
    model_rank_cols_all = [m for m in MODEL_NAME_MAP if m in rank_pivot_all.columns]
    mpkt_all_recs = []
    for luxr, group in rank_pivot_all.groupby("luxr"):
        mpkt_all_recs.append({
            "LuxR_name": luxr,
            "mpKT": compute_avg_agreement_all(group, model_rank_cols_all),
        })
    mpkt_all_df = pd.DataFrame(mpkt_all_recs)

    excel_models = [
        ("AF3", "af3"),
        ("Boltz2", "boltz2"),
        ("Boltz2Affinity", "boltz2affinity"),
        ("RF3", "rf3"),
        ("ENSEMBLE", None),
    ]
    for display_name, model_key in excel_models:
        excel_path = full_dir / f"{display_name}.xlsx"
        build_ranked_excel(
            model_key, rankings, ensemble, labels, ligands, mpkt_all_df, excel_path
        )



# =============================================================================
# CLI and main
# =============================================================================


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="LuxR benchmark analysis: collect (optional), postprocess, generate CSVs.",
        epilog="If --input-csv is provided, collection is skipped and that file is used for postprocess.",
    )
    # Collection (optional when --input-csv is set)
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=None,
        help="Pre-collected processed CSV; if set, skip collection and use this for postprocess.",
    )
    parser.add_argument(
        "--basepath",
        type=Path,
        default=None,
        help="Base directory for benchmark runs (required only when collection runs).",
    )
    parser.add_argument(
        "--af3-basepath",
        type=Path,
        default=None,
        help="Optional base directory for AF3 results (used only when collection runs).",
    )
    parser.add_argument(
        "--af3-iptm-source",
        choices=("base", "json"),
        default="json",
        help="For collection: use AF3 ipTM from base-folder CSV (base) or compute from JSON (json, default).",
    )
    # Workflow and output paths
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        required=True,
        help="Directory for processed CSV (when collecting) and postprocess outputs.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Root directory for generate outputs (ranking_performance, etc.).",
    )
    parser.add_argument(
        "--processed-name",
        default=PROCESSED_CSV_NAME,
        help="Filename for processed benchmark CSV under workflow-dir (when collecting).",
    )
    # Data paths
    parser.add_argument(
        "--labels-path",
        type=Path,
        required=True,
        help="Path to labels CSV (used by postprocess and generate).",
    )
    parser.add_argument(
        "--ligands-path",
        type=Path,
        required=True,
        help="Path to ligands CSV (used by generate).",
    )
    parser.add_argument(
        "--ablation-path",
        type=Path,
        default=None,
        help="Override ablation table path for generate (default: workflow-dir/ablation_table.csv).",
    )
    parser.add_argument(
        "--no-leakage-csv",
        type=Path,
        default=None,
        help="CSV with no-leakage LuxR subset (refahl_id column). "
             "When provided, conf=1 outputs get additional _no_leakage variants.",
    )
    # Postprocess skip flags
    parser.add_argument("--skip-rankings", action="store_true")
    parser.add_argument("--skip-ensemble", action="store_true")
    parser.add_argument("--skip-ablation", action="store_true")
    parser.add_argument("--skip-done-json", action="store_true")
    parser.add_argument(
        "--skip-generate",
        action="store_true",
        help="Stop after postprocess; don't read rankings/ensemble CSVs or emit results/.",
    )
    parser.add_argument(
        "--expected-ligand-count",
        type=int,
        default=None,
        help="Override the 'all predictions present' guard in compute_all_rankings. "
             "Defaults to the number of rows in the ligands CSV (58 for the paper benchmark).",
    )
    parser.add_argument(
        "--aggregation",
        type=str,
        default=DEFAULT_AGGREGATION,
        choices=AGGREGATIONS,
        help=(
            "How to aggregate per-model ranks into the ensemble score. "
            "'rrf' (default) is Reciprocal Rank Fusion, robust to one bad "
            "model and aligned with top-k success. 'sum'/'mean'/'median'/"
            "'min'/'max' are simple reductions over per-model ranks."
        ),
    )
    parser.add_argument(
        "--rrf-k",
        type=int,
        default=DEFAULT_RRF_K,
        help=f"RRF damping constant (default {DEFAULT_RRF_K}; only used when --aggregation=rrf).",
    )
    return parser.parse_args()


def main() -> None:
    import warnings

    warnings.warn(
        "copal.analysis.pipeline.main() is deprecated; "
        "use 'copal analyze' (Hydra-backed) instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    args = parse_args()
    workflow_dir = args.workflow_dir

    if args.input_csv is not None:
        processed_csv = Path(args.input_csv)
        if not processed_csv.exists():
            raise FileNotFoundError(f"Pre-collected CSV not found: {processed_csv}")
        print(f"Skipping collection; using pre-collected CSV: {processed_csv}")
    else:
        workflow_dir.mkdir(parents=True, exist_ok=True)
        all_df = collect_records(
            args.basepath,
            af3_basepath=args.af3_basepath,
            af3_iptm_source=args.af3_iptm_source,
        )
        processed_csv = workflow_dir / args.processed_name
        all_df.to_csv(processed_csv, index=False)
        print(f"Wrote {processed_csv} ({len(all_df)} rows)")

    if args.expected_ligand_count is not None:
        expected_ligand_count = args.expected_ligand_count
    else:
        expected_ligand_count = len(pd.read_csv(args.ligands_path))
        print(
            f"expected_ligand_count = {expected_ligand_count} (from {args.ligands_path.name})"
        )

    run_postprocess(
        processed_csv,
        workflow_dir,
        args.labels_path,
        skip_rankings=args.skip_rankings,
        skip_ensemble=args.skip_ensemble,
        skip_ablation=args.skip_ablation,
        skip_done_json=args.skip_done_json,
        expected_ligand_count=expected_ligand_count,
        aggregation=args.aggregation,
        rrf_k=args.rrf_k,
    )

    if args.skip_generate:
        print("skip-generate set; stopping after postprocess.")
        return

    rankings_csv = workflow_dir / RANKINGS_CSV_NAME
    ensemble_csv = workflow_dir / ENSEMBLE_CSV_NAME
    if not rankings_csv.exists() or not ensemble_csv.exists():
        missing = [p.name for p in (rankings_csv, ensemble_csv) if not p.exists()]
        print(
            f"Skipping generate: required workflow CSV(s) missing ({', '.join(missing)}). "
            "This typically means --skip-rankings/--skip-ensemble was set, or the "
            "benchmark has fewer than 58 luxr_ahl pairs per receptor."
        )
        return

    run_generate(
        workflow_dir,
        args.output_root,
        args.labels_path,
        args.ligands_path,
        ablation_path=args.ablation_path,
        no_leakage_path=args.no_leakage_csv,
        aggregation=args.aggregation,
        rrf_k=args.rrf_k,
    )


if __name__ == "__main__":
    main()
