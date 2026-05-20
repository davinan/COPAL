"""AlphaFold 3 manager — writes AF3 JSON and emits the per-chunk shell command."""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np

from copal.runner.base import Manager, resolve_script_path


def _abs_msa_path(msa_path):
    """Absolutize an MSA path for the AF3 JSON.

    AF3 resolves ``unpairedMsaPath`` relative to the input JSON's own
    directory, not the CWD — so a repo-root-relative path (e.g.
    ``data/benchmark_msas/N.a3m``) would not be found. Pass-through when
    falsy (no MSA configured).
    """
    if not msa_path:
        return msa_path
    return str(Path(msa_path).expanduser().resolve())


class AF3Manager(Manager):
    MODEL_ID = "af3"
    SCRIPT_ENV_VAR = "COPAL_AF3_SCRIPT"
    SCRIPT_DEFAULT = "af3.sh"

    def __init__(self, run_name: str, out_dir: Path, in_dir: Path):
        super().__init__(run_name, "AF3", datetime.now().strftime("%Y-%m-%d"), out_dir, in_dir)

    # -- post-prediction accessors -------------------------------------

    def get_structure(self, name: str, chunk: str) -> Path:
        return self.out_dir / chunk / name.lower() / f"{name.lower()}_model.cif"

    def check_done(self, name: str, chunk: str) -> bool:
        return (self.out_dir / chunk / name.lower() / f"{name.lower()}_model.cif").exists()

    def get_pae_interface(self, name, chain_i: str, chain_j: str):
        chunk, name = name
        recs = []
        for sample_dir in (self.out_dir / chunk / name.lower()).iterdir():
            if not sample_dir.is_dir():
                continue
            batch_idx = int(sample_dir.name.split("-")[-1])
            conf_path = sample_dir / f"{name.lower()}_{sample_dir.name}_confidences.json"
            with open(conf_path) as fh:
                af3_conf = json.load(fh)
            pae = np.array(af3_conf["pae"])
            chain_id = np.array(af3_conf["token_chain_ids"])
            is_chain_of_interest = (chain_id == chain_i) | (chain_id == chain_j)
            is_interface_of_interest = is_chain_of_interest[:, None] | is_chain_of_interest[None, :]
            is_interface = chain_id[:, None] != chain_id[None, :]
            mask = is_interface & is_interface_of_interest
            recs.append({
                "name": name,
                "model": self.MODEL_ID,
                "batch_idx": batch_idx,
                "min_pae_interface": pae[mask].min(),
                "pae_interface": pae[mask].mean(),
                "pred_file": sample_dir / f"{name.lower()}_{sample_dir.name}_model.cif",
            })
        return recs

    def get_completed_names(self):
        return list({
            (x.parent.parent.name, x.name.split("_")[0])
            for x in self.out_dir.glob("*/*/*.cif")
            if (x.parent / f"{x.parent.name}_confidences.json").exists()
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
        return [{
            "protein": {
                "id": ids,
                "sequence": prot_info["sequence"],
                "unpairedMsaPath": _abs_msa_path(prot_info["msa_path"]),
                "pairedMsa": "",
                "templates": [],
            }
        }]

    def make_ligand_input(self, ids, lig_info):
        return [{"ligand": {"id": ids, "smiles": lig_info}}]

    def make_dna_input(self, ids, dna_info: dict):
        return [{"dna": {"id": ids, "sequence": dna_info["sequence"]}}]

    def make_rna_input(self, ids, rna_info: dict):
        return [{
            "rna": {
                "id": ids,
                "sequence": rna_info["sequence"],
                "unpairedMsaPath": _abs_msa_path(rna_info["msa_path"]),
                "pairedMsa": "",
                "templates": [],
            }
        }]

    def assemble_component_inputs(self, single_input_components, complex_name: str, extra=None):
        return (
            complex_name,
            {
                "name": complex_name,
                "modelSeeds": [3],
                "sequences": single_input_components,
                "dialect": "alphafold3",
                "version": 3,
            },
        )

    # -- IO --------------------------------------------------------------

    def save_inputs(self, inputs):
        for chunk_id, (complex_name, payload) in inputs:
            (self.in_dir / str(chunk_id)).mkdir(exist_ok=True)
            with open(self.in_dir / str(chunk_id) / f"{complex_name}.json", "w") as fh:
                json.dump(payload, fh)

    def save_single_input(self, pair):
        chunk_id, (complex_name, payload) = pair
        (self.in_dir / str(chunk_id)).mkdir(exist_ok=True)
        with open(self.in_dir / str(chunk_id) / f"{complex_name}.json", "w") as fh:
            json.dump(payload, fh)

    def get_input_path(self, complex_name: str, chunk: str) -> Path:
        return self.in_dir / chunk / f"{complex_name}.json"
