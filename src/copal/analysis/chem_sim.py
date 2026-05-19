"""Pairwise ligand chemical-similarity methods + matrix precomputation.

Three MCS-Tanimoto variants over heavy atoms, all using ringMatchesRingOnly=True
and atomCompare=CompareElements:

  - mcs_order:     bondCompare=CompareOrder. Strict on bond orders. Splits the
                   MCS at a single unmatched bond, so saturated/unsaturated chain
                   pairs (e.g. C12 vs C12:1) drop to ~0.5. Historical baseline.
  - mcs_any:       bondCompare=CompareAny.  Ignores bond order entirely.
                   Collapses bond-order-only differences (C12 vs C12:1, 3OC12
                   vs 3OHC12) to 1.0 — over-permissive.
  - mcs_penalized: CompareAny matching, then deduct one atom-equivalent per
                   bond whose order differs between the matched substructures.
                   Smooth degradation on bond-order changes; preserves chain-
                   length and missing-atom sensitivity.

The matrix builder writes a symmetric N×N DataFrame (AHLid-indexed, diagonal=1.0)
suitable for cheap O(1) lookups during benchmark analysis.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from copal.analysis.pipeline import get_rdkit


def _same_canonical(mol_a, mol_b) -> bool:
    Chem, _ = get_rdkit()
    return Chem.MolToSmiles(mol_a, isomericSmiles=True) == Chem.MolToSmiles(
        mol_b, isomericSmiles=True
    )


def _find_mcs(mol_a, mol_b, *, bond_mode: str, timeout: int):
    _, rdFMCS = get_rdkit()
    bond_cmp = {
        "order": rdFMCS.BondCompare.CompareOrder,
        "any": rdFMCS.BondCompare.CompareAny,
    }[bond_mode]
    return rdFMCS.FindMCS(
        [mol_a, mol_b],
        timeout=timeout,
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=bond_cmp,
    )


def mcs_similarity(mol_a, mol_b, *, bond_mode: str = "order", timeout: int = 10) -> float:
    """MCS-Tanimoto over heavy atoms (no bond-mismatch penalty)."""
    if mol_a is None or mol_b is None:
        return float("nan")
    if _same_canonical(mol_a, mol_b):
        return 1.0
    mcs = _find_mcs(mol_a, mol_b, bond_mode=bond_mode, timeout=timeout)
    if mcs is None or mcs.numAtoms == 0:
        return 0.0
    common = float(mcs.numAtoms)
    a = float(mol_a.GetNumHeavyAtoms())
    b = float(mol_b.GetNumHeavyAtoms())
    denom = a + b - common
    if denom <= 0:
        return float("nan")
    return common / denom


def mcs_similarity_penalized(mol_a, mol_b, *, timeout: int = 10) -> float:
    """CompareAny match + 1-atom-equivalent penalty per bond-order mismatch.

    Finds the maximum common atomic skeleton ignoring bond order, then subtracts
    one atom from the effective common count for each MCS bond whose order
    differs between the two molecules. Yields a smooth degradation on
    bond-order changes (one C=C vs C-C ≈ -1 atom equivalent) without
    fragmenting the chain.
    """
    Chem, _ = get_rdkit()
    if mol_a is None or mol_b is None:
        return float("nan")
    if _same_canonical(mol_a, mol_b):
        return 1.0
    mcs = _find_mcs(mol_a, mol_b, bond_mode="any", timeout=timeout)
    if mcs is None or mcs.numAtoms == 0:
        return 0.0

    patt = Chem.MolFromSmarts(mcs.smartsString)
    if patt is None:
        return 0.0
    match_a = mol_a.GetSubstructMatch(patt)
    match_b = mol_b.GetSubstructMatch(patt)
    if not match_a or not match_b:
        return 0.0

    mismatches = 0
    for bond in patt.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ba = mol_a.GetBondBetweenAtoms(match_a[i], match_a[j])
        bb = mol_b.GetBondBetweenAtoms(match_b[i], match_b[j])
        if ba is None or bb is None:
            continue
        if ba.GetBondType() != bb.GetBondType():
            mismatches += 1

    eff = max(0.0, float(mcs.numAtoms) - float(mismatches))
    a = float(mol_a.GetNumHeavyAtoms())
    b = float(mol_b.GetNumHeavyAtoms())
    denom = a + b - eff
    if denom <= 0:
        return float("nan")
    return eff / denom


METHODS: dict[str, Callable] = {
    "mcs_order":     lambda a, b: mcs_similarity(a, b, bond_mode="order"),
    "mcs_any":       lambda a, b: mcs_similarity(a, b, bond_mode="any"),
    "mcs_penalized": mcs_similarity_penalized,
}


def compute_pairwise_matrix(
    ligands_df: pd.DataFrame,
    method: str,
    *,
    smiles_col: str = "smi",
    id_col: str = "AHLid",
) -> pd.DataFrame:
    """Compute the symmetric pairwise similarity matrix for a ligand pool."""
    if method not in METHODS:
        raise ValueError(
            f"Unknown method {method!r}; available: {sorted(METHODS)}"
        )
    fn = METHODS[method]
    Chem, _ = get_rdkit()

    ids: list[str] = []
    mols: list = []
    for aid, smi in zip(ligands_df[id_col], ligands_df[smiles_col]):
        if not isinstance(smi, str):
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        ids.append(str(aid).strip())
        mols.append(mol)

    n = len(ids)
    mat = np.eye(n, dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            s = fn(mols[i], mols[j])
            mat[i, j] = mat[j, i] = s

    df = pd.DataFrame(mat, index=ids, columns=ids)
    df.index.name = id_col
    return df


def write_pairwise_matrix(
    ligands_df: pd.DataFrame,
    method: str,
    out_path: Path,
    *,
    smiles_col: str = "smi",
    id_col: str = "AHLid",
) -> pd.DataFrame:
    """Compute + write the matrix; return it for in-memory use."""
    df = compute_pairwise_matrix(
        ligands_df, method, smiles_col=smiles_col, id_col=id_col
    )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path)
    return df
