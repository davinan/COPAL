"""Boltz2 manager — writes Boltz2 YAML and emits the per-chunk shell command."""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

from copal.runner.base import Manager, resolve_script_path


class Boltz2Manager(Manager):
    MODEL_ID = "boltz2"
    SCRIPT_ENV_VAR = "COPAL_BOLTZ2_SCRIPT"
    SCRIPT_DEFAULT = "boltz2.sh"

    def __init__(self, run_name: str, out_dir: Path, in_dir: Path, num_models: int = 5):
        super().__init__(run_name, "Boltz2", datetime.now().strftime("%Y-%m-%d"), out_dir, in_dir)
        self.num_models = num_models

    # -- post-prediction accessors -------------------------------------

    def _chain_ids_from_input(self, input_dict: dict) -> np.ndarray:
        """Expand per-sequence chain IDs to a per-token array (needed for PAE masks)."""
        from rdkit import Chem
        chain_ids: list = []
        for sequence in input_dict["sequences"]:
            entity_type = next(iter(sequence.keys()))
            chain_id = sequence[entity_type]["id"]
            if entity_type in {"protein", "rna", "dna"}:
                L = len(sequence[entity_type]["sequence"])
            elif entity_type == "ligand":
                L = Chem.MolFromSmiles(sequence["ligand"]["smiles"]).GetNumAtoms()
            else:
                raise ValueError(f"Unknown entity type: {entity_type}")
            chain_ids.extend([chain_id] * L)
        return np.array(chain_ids)

    def check_done(self, name: str, chunk: str) -> bool:
        return (
            self.out_dir / chunk / f"boltz_results_{chunk}/"
            / "predictions" / name / f"{name}_model_0.cif"
        ).exists()

    def get_structure(self, name: str, chunk: str, return_index: bool = False):
        best_plddt = 0.0
        best_index = 0
        base = self.out_dir / chunk / f"boltz_results_{chunk}/" / "predictions" / name
        best_path = base / f"{name}_model_0.cif"
        for i in range(self.num_models):
            cif_path = base / f"{name}_model_{i}.cif"
            plddt_path = base / f"plddt_{name}_model_{i}.npz"
            if not plddt_path.exists():
                continue
            plddt = np.load(plddt_path)["plddt"].mean()
            if plddt > best_plddt:
                best_plddt = plddt
                best_index = i
                best_path = cif_path
        return (best_path, best_index) if return_index else best_path

    def get_pae_interface(self, name, chain_i: str, chain_j: str):
        chunk, name = name
        recs = []
        for model_index in range(self.num_models):
            pae_path = (
                self.out_dir / chunk / f"boltz_results_{chunk}/" / "predictions"
                / name / f"pae_{name}_model_{model_index}.npz"
            )
            yaml_path = self.in_dir / chunk / f"{name}.yaml"
            if not yaml_path.exists():
                yaml_path = self.old_in_dir / chunk / f"{name}.yaml"
            with open(yaml_path) as fh:
                input_dict = yaml.safe_load(fh)
            chain_ids = self._chain_ids_from_input(input_dict)
            is_interface = chain_ids[:, None] != chain_ids[None, :]
            is_chain_of_interest = (chain_ids == chain_i) | (chain_ids == chain_j)
            is_interface_of_interest = is_chain_of_interest[:, None] | is_chain_of_interest[None, :]
            mask = is_interface & is_interface_of_interest
            pae = np.load(pae_path)["pae"]
            recs.append({
                "name": name,
                "model": self.MODEL_ID,
                "batch_idx": model_index,
                "min_pae_interface": pae[mask].min(),
                "pae_interface": pae[mask].mean(),
                "pred_file": pae_path.parent / f"{name}_model_{model_index}.cif",
            })

            affinity_path = (
                self.out_dir / chunk / f"boltz_results_{chunk}/" / "predictions"
                / name / f"affinity_{name}.json"
            )
            if affinity_path.exists() and model_index == 0:
                with open(affinity_path) as fh:
                    data = json.load(fh)
                recs.append({
                    "name": name,
                    "model": "boltz2affinity",
                    "batch_idx": model_index,
                    "min_pae_interface": data["affinity_pred_value"],
                    "pae_interface": data["affinity_probability_binary"],
                    "pred_file": pae_path.parent / f"{name}_model_{model_index}.cif",
                })
        return recs

    def get_completed_names(self):
        return list({
            (x.parent.parent.parent.parent.name, x.name.split("_")[0])
            for x in self.out_dir.glob("*/*/*/*/*_model_0.cif")
        })

    # -- command generation --------------------------------------------

    def get_model_command(self) -> list[str]:
        script = resolve_script_path(self.SCRIPT_ENV_VAR, self.SCRIPT_DEFAULT)
        cmds = []
        for chunk_dir in self.in_dir.iterdir():
            if not chunk_dir.is_dir() or not any(chunk_dir.iterdir()):
                continue
            (self.out_dir / chunk_dir.name).mkdir(exist_ok=True)
            cmds.append(f"{script} {chunk_dir} {self.out_dir / chunk_dir.name}")
        return cmds

    def organize_data(self):
        pass

    # -- input construction --------------------------------------------

    def make_protein_input(self, ids, prot_info: dict):
        return [
            {"protein": {"id": i, "sequence": prot_info["sequence"], "msa": prot_info["msa_path"]}}
            for i in ids
        ]

    def make_ligand_input(self, ids, lig_info):
        return [{"ligand": {"id": i, "smiles": lig_info}} for i in ids]

    def make_dna_input(self, ids, dna_info: dict):
        return [{"dna": {"id": i, "sequence": dna_info["sequence"]}} for i in ids]

    def make_rna_input(self, ids, rna_info: dict):
        return [{"rna": {"id": i, "sequence": rna_info["sequence"]}} for i in ids]

    def assemble_component_inputs(self, single_input_components, complex_name: str, extra=None):
        payload = {"sequences": single_input_components}
        if extra is not None:
            payload.update(extra)
        return (complex_name, payload)

    # -- IO --------------------------------------------------------------

    def save_inputs(self, inputs):
        for chunk_id, (complex_name, payload) in inputs:
            (self.in_dir / str(chunk_id)).mkdir(exist_ok=True)
            with open(self.in_dir / str(chunk_id) / f"{complex_name}.yaml", "w") as fh:
                yaml.dump(payload, fh)

    def save_single_input(self, pair):
        chunk_id, (complex_name, payload) = pair
        (self.in_dir / str(chunk_id)).mkdir(exist_ok=True)
        with open(self.in_dir / str(chunk_id) / f"{complex_name}.yaml", "w") as fh:
            yaml.dump(payload, fh)

    def get_input_path(self, complex_name: str, chunk: str) -> Path:
        return self.in_dir / chunk / f"{complex_name}.yaml"
