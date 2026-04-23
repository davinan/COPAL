"""RosettaFold3 (Foundry) manager — writes RF3 JSON and emits the per-chunk shell command."""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd

from copal.runner.base import Manager, resolve_script_path


class RF3Manager(Manager):
    MODEL_ID = "rf3"
    SCRIPT_ENV_VAR = "COPAL_RF3_SCRIPT"
    SCRIPT_DEFAULT = "rf3.sh"

    def __init__(self, run_name: str, out_dir: Path, in_dir: Path):
        super().__init__(run_name, "RF3", datetime.now().strftime("%Y-%m-%d"), out_dir, in_dir)

    # -- post-prediction accessors -------------------------------------

    def check_done(self, name: str, chunk: str) -> bool:
        return (self.out_dir / chunk / f"{name}.score").exists()

    def get_structure(self, name: str, chunk: str) -> Path:
        return self.out_dir / chunk / f"{name}_model_0.cif.gz"

    def get_pae_interface(self, name, chain_i: str, chain_j: str):
        chunk, name = name
        df = pd.read_csv(self.out_dir / chunk / f"{name}.score")
        recs = []
        for _, row in df[
            (df.chain_i_interface == f"{chain_i}_1") & (df.chain_j_interface == f"{chain_j}_1")
        ].iterrows():
            recs.append({
                "name": name,
                "model": self.MODEL_ID,
                "batch_idx": row["batch_idx"],
                "min_pae_interface": row["min_pae_interface"],
                "pae_interface": row["pae_interface"],
                "pred_file": self.out_dir / chunk / f"{name}_model_{row['batch_idx']}.cif.gz",
            })
        return recs

    def get_completed_names(self):
        return list({
            (x.parent.name, x.name.split(".score")[0])
            for x in self.out_dir.glob("*/*.score")
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
            {"chain_id": i, "seq": prot_info["sequence"], "msa_path": prot_info["msa_path"]}
            for i in ids
        ]

    def make_ligand_input(self, ids, lig_info):
        return [{"chain_id": i, "smiles": lig_info} for i in ids]

    def make_dna_input(self, ids, dna_info: dict):
        return [{"chain_id": i, "seq": dna_info["sequence"]} for i in ids]

    def make_rna_input(self, ids, rna_info: dict):
        return [
            {"chain_id": i, "seq": rna_info["sequence"], "msa_path": rna_info["msa_path"]}
            for i in ids
        ]

    def assemble_component_inputs(self, single_input_components, complex_name: str, extra=None):
        return (complex_name, {"name": complex_name, "components": single_input_components})

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
