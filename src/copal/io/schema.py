"""Receptor and ligand input schemas with friendly validation.

A `receptors.csv` file MUST have the columns:
    receptor_id  — unique string ID for the protein target
    seq          — single-chain amino-acid sequence

Optional columns:
    msa_path        — path to a pre-built MSA (.a3m). If absent, the model
                      runner generates an MSA at inference time.
    chain_count     — integer number of identical chains (default 1, useful
                      for homo-oligomers like LuxR dimers)
    template_paths  — semicolon-separated list of template CIFs

A `ligands.csv` file MUST have the columns:
    ligand_id  — unique string ID
    smiles     — SMILES string

Optional columns:
    inchi, name, formula, mol_weight, etc. (passed through but unused).

Backward-compatible aliases (auto-renamed on load):
    receptor_id ← LuxR_name, refahl_id, target
    seq         ← sequence
    ligand_id   ← AHLid, ahl_id, abbreviation
    smiles      ← Isomeric SMILES, smi
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


_RECEPTOR_ALIASES = {
    "LuxR_name": "receptor_id",
    "refahl_id": "receptor_id",
    "target": "receptor_id",
    "sequence": "seq",
}

_LIGAND_ALIASES = {
    "AHLid": "ligand_id",
    "ahl_id": "ligand_id",
    "abbreviation": "ligand_id",
    "Isomeric SMILES": "smiles",
    "smi": "smiles",
    "SMILES": "smiles",
}


def _apply_aliases(df: pd.DataFrame, aliases: dict[str, str]) -> pd.DataFrame:
    """Rename alias columns to their canonical name, preferring the canonical
    name if it already exists, and using the first present alias otherwise.
    Other already-present aliases for the same canonical name are dropped to
    avoid duplicate columns."""
    cols = list(df.columns)
    keep = list(cols)
    rename: dict[str, str] = {}
    for canonical in set(aliases.values()):
        if canonical in cols:
            for alias, target in aliases.items():
                if target == canonical and alias in cols and alias != canonical:
                    keep.remove(alias)
            continue
        for alias, target in aliases.items():
            if target == canonical and alias in cols:
                rename[alias] = canonical
                for other_alias, other_target in aliases.items():
                    if other_target == canonical and other_alias in cols and other_alias != alias:
                        keep.remove(other_alias)
                break
    df = df[keep]
    return df.rename(columns=rename) if rename else df


def _require_columns(df: pd.DataFrame, required: Iterable[str], label: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label}: missing required column(s) {missing}. "
            f"Got columns: {list(df.columns)}"
        )


def _check_unique(df: pd.DataFrame, col: str, label: str) -> None:
    dupes = df[col][df[col].duplicated()].unique().tolist()
    if dupes:
        raise ValueError(f"{label}: column {col!r} has duplicate values: {dupes}")


def load_receptors_csv(path: str | Path) -> pd.DataFrame:
    """Load and validate a receptors CSV. Returns a normalized DataFrame."""
    df = pd.read_csv(path)
    df = _apply_aliases(df, _RECEPTOR_ALIASES)
    _require_columns(df, ["receptor_id", "seq"], f"receptors CSV {path}")
    df["receptor_id"] = df["receptor_id"].astype(str).str.strip()
    df["seq"] = df["seq"].astype(str).str.strip()
    df = df[df["seq"].str.len() > 0].copy()
    _check_unique(df, "receptor_id", f"receptors CSV {path}")
    return df


def load_ligands_csv(path: str | Path) -> pd.DataFrame:
    """Load and validate a ligands CSV. Returns a normalized DataFrame."""
    df = pd.read_csv(path)
    df = _apply_aliases(df, _LIGAND_ALIASES)
    _require_columns(df, ["ligand_id", "smiles"], f"ligands CSV {path}")
    df["ligand_id"] = df["ligand_id"].astype(str).str.strip()
    df["smiles"] = df["smiles"].astype(str).str.strip()
    df = df[df["smiles"].str.len() > 0].copy()
    _check_unique(df, "ligand_id", f"ligands CSV {path}")
    return df
