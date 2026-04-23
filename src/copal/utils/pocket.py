"""Identify pocket residues by mapping a query sequence onto a reference
protein-ligand structure (default: 4Y17, an *E. coli* SdiA / 3OC8 complex).

The default reference CIF lives in `data/reference/` and can be overridden
with the `COPAL_POCKET_REF_CIF` env var or the `ref_model` argument.
"""
from __future__ import annotations

import os
from pathlib import Path

import biotite.sequence as bseq
import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray


_DEFAULT_REF_CIF = (
    Path(__file__).resolve().parents[3] / "data" / "reference" / "4y17.cif"
)


def _resolve_ref_cif(ref_model: Path | str | None) -> Path:
    if ref_model is not None:
        return Path(ref_model)
    env = os.environ.get("COPAL_POCKET_REF_CIF")
    if env:
        return Path(env)
    return _DEFAULT_REF_CIF


def get_contact_annotations(
    subj: AtomArray,
    obj: AtomArray,
    thres: float,
    annot_subj: np.ndarray,
    annot_obj: np.ndarray,
):
    cell_list = struc.CellList(
        obj[~np.any(np.isnan(obj.coord), axis=-1)], cell_size=thres,
    )
    contacts = cell_list.get_atoms(subj.coord, radius=thres, as_mask=True)
    subj_indices, obj_indices = np.where(contacts)
    pairs = list(zip(annot_subj[subj_indices].tolist(), annot_obj[obj_indices].tolist()))
    seen: set = set()
    out_subj: list = []
    out_obj: list = []
    for pair in pairs:
        if pair in seen:
            continue
        seen.add(pair)
        out_subj.append(pair[0])
        out_obj.append(pair[1])
    return out_subj, out_obj


def get_pocket_res_ids(query_seq: str, ref: AtomArray, is_dimer: bool) -> list[int]:
    matrix = bseq.align.SubstitutionMatrix.std_protein_matrix()
    ref = ref[ref.atom_name == "CA"]
    ref_seq = struc.to_sequence(ref)[0][0]

    alignment = bseq.align.align_optimal(
        bseq.ProteinSequence(query_seq), ref_seq, matrix
    )[0]
    trace = [(x, y) for x, y in alignment.trace if x != -1 and y != -1]

    pocket_ids: list[int] = []
    for x, y in trace:
        if ref[y].is_pocket:
            pocket_ids.append(x)
            if is_dimer:
                pocket_ids.append(x + len(query_seq))
    return list(set(pocket_ids))


def load_reference_pocket(
    ref_model: Path | str | None = None,
    contact_threshold: float = 8.0,
):
    """Load the reference structure and tag pocket residues (residues whose
    CA-atoms come within `contact_threshold` Å of any non-polymer ligand)."""
    from cifutils import parse  # local import: keeps cifutils optional
    from cifutils.enums import ChainType

    ref_path = _resolve_ref_cif(ref_model)
    ref = parse(ref_path)["assemblies"]["1"][0]
    ref = ref[np.isin(ref.chain_id, ["A", "D"])]
    ref_prot = ref[ref.chain_type == ChainType.POLYPEPTIDE_L]
    ref_ligs = ref[ref.chain_type == ChainType.NON_POLYMER]

    pocket_res_ids, _ = get_contact_annotations(
        ref_prot, ref_ligs, contact_threshold, ref_prot.res_id, ref_ligs.res_id,
    )
    ref_prot.set_annotation("is_pocket", np.isin(ref_prot.res_id, pocket_res_ids))
    return ref_prot
