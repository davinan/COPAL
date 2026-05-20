"""Abstract `Manager` for a single model in the ensemble.

Each subclass (AF3, Boltz2, RF3) knows how to:
    - turn protein / ligand / RNA / DNA rows into its native input format
    - assemble them into the model's complex input file (JSON, YAML, …)
    - emit the shell command that runs one chunk of inputs
    - check whether a given prediction exists on disk

Inputs come from up to four CSVs (protein / ligand / RNA / DNA); only
`protein` and `ligand` are used in the default COPAL pipeline.
"""
from __future__ import annotations

import os
import pickle
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from copal.io.schema import load_ligands_csv, load_receptors_csv
from copal.runner.sbatch import make_sbatch_file


class Manager(ABC):
    """Abstract base class for one model's input / command management."""

    def __init__(
        self,
        run_name: str,
        model_name: str,
        date: str,
        out_dir: Path,
        in_dir: Path,
    ):
        self.run_name = run_name
        self.model_name = model_name
        self.date = date
        self.out_dir = Path(out_dir) / model_name
        self.in_dir = Path(in_dir) / model_name
        # Created lazily the first time `skip_existing_structures` archives an
        # already-done prediction's input file — see `make_inputs` below.
        self.old_in_dir = Path(in_dir) / f"{model_name}_done"
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.in_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Abstract surface — subclasses implement these
    # ------------------------------------------------------------------

    @abstractmethod
    def organize_data(self): ...

    @abstractmethod
    def get_model_command(self): ...

    @abstractmethod
    def get_pae_interface(self, name: str, chain_i: str, chain_j: str): ...

    @abstractmethod
    def get_completed_names(self): ...

    @abstractmethod
    def get_structure(self, name: str): ...

    @abstractmethod
    def make_protein_input(self, ids, prot_info: dict): ...

    @abstractmethod
    def make_ligand_input(self, ids, lig_info: dict): ...

    @abstractmethod
    def make_rna_input(self, ids, rna_info: dict): ...

    @abstractmethod
    def make_dna_input(self, ids, dna_info: dict): ...

    @abstractmethod
    def assemble_component_inputs(self, single_input_components, complex_name: str, extra=None): ...

    @abstractmethod
    def save_inputs(self, inputs): ...

    @abstractmethod
    def save_single_input(self, inputs): ...

    @abstractmethod
    def get_input_path(self, complex_name: str, chunk): ...

    @abstractmethod
    def check_done(self, name: str, chunk: str) -> bool: ...

    # ------------------------------------------------------------------
    # Shared helpers
    # ------------------------------------------------------------------

    @staticmethod
    def next_chain_ids(n: int, start_index: int):
        """Yield A, B, C, ... AA, AB, ... chain IDs."""
        ids: list[str] = []
        i = start_index
        while len(ids) < n:
            q, r = divmod(i, 26)
            base = chr(ord("A") + r)
            ids.append(base if q == 0 else chr(ord("A") + (q - 1)) + base)
            i += 1
        return ids, i

    def load_csvs(
        self,
        protein_csv: Path | None = None,
        protein_name_col: str = "protein_name",
        rna_csv: Path | None = None,
        rna_name_col: str = "rna_name",
        dna_csv: Path | None = None,
        dna_name_col: str = "dna_name",
        ligand_csv: Path | None = None,
        ligand_name_col: str = "lig_name",
    ) -> dict:
        prot_df = (
            load_receptors_csv(protein_csv)     # alias-aware: refahl_id/LuxR_name -> receptor_id, sequence -> seq
            if protein_csv
            else pd.DataFrame(columns=[protein_name_col, "sequence", "msa_path"])
        )
        lig_df = (
            load_ligands_csv(ligand_csv)        # alias-aware: AHLid/smi -> ligand_id/smiles
            if ligand_csv
            else pd.DataFrame(columns=[ligand_name_col, "smiles"])
        )
        rna_df = (
            pd.read_csv(rna_csv)
            if rna_csv
            else pd.DataFrame(columns=[rna_name_col, "sequence", "msa_path"])
        )
        dna_df = (
            pd.read_csv(dna_csv)
            if dna_csv
            else pd.DataFrame(columns=[dna_name_col, "sequence"])
        )

        lig_df = lig_df.rename(columns={"smi": "smiles", ligand_name_col: "ligand_name"})
        prot_df = prot_df.rename(columns={"seq": "sequence", protein_name_col: "protein_name"})
        rna_df = rna_df.rename(columns={rna_name_col: "rna_name"})
        dna_df = dna_df.rename(columns={dna_name_col: "dna_name"})

        for col in ("protein_name", "sequence", "msa_path"):
            if col not in prot_df.columns:
                prot_df[col] = None
        for col in ("ligand_name", "smiles"):
            if col not in lig_df.columns:
                lig_df[col] = None
        for col in ("rna_name", "sequence", "msa_path"):
            if col not in rna_df.columns:
                rna_df[col] = None
        for col in ("dna_name", "sequence"):
            if col not in dna_df.columns:
                dna_df[col] = None

        prot_info = {
            str(row["protein_name"]).strip(): {
                "sequence": str(row["sequence"]).strip() if pd.notna(row["sequence"]) else None,
                "msa_path": str(row["msa_path"]).strip() if pd.notna(row["msa_path"]) else None,
            }
            for _, row in prot_df.iterrows()
            if pd.notna(row.get("protein_name"))
        }
        lig_info = {
            str(row["ligand_name"]).strip(): str(row["smiles"]).strip()
            for _, row in lig_df.iterrows()
            if pd.notna(row.get("ligand_name"))
        }
        rna_info = {
            str(row["rna_name"]).strip(): {
                "sequence": str(row["sequence"]).strip() if pd.notna(row["sequence"]) else None,
                "msa_path": str(row["msa_path"]).strip() if pd.notna(row["msa_path"]) else None,
            }
            for _, row in rna_df.iterrows()
            if pd.notna(row.get("rna_name"))
        }
        dna_info = {
            str(row["dna_name"]).strip(): {
                "sequence": str(row["sequence"]).strip() if pd.notna(row["sequence"]) else None,
            }
            for _, row in dna_df.iterrows()
            if pd.notna(row.get("dna_name"))
        }
        return {"prot": prot_info, "lig": lig_info, "rna": rna_info, "dna": dna_info}

    def _input_maker_factory(self, component_type: str, info: dict):
        if component_type in {"prot", "protein"}:
            return self.make_protein_input, info["prot"]
        if component_type in {"lig", "ligand"}:
            return self.make_ligand_input, info["lig"]
        if component_type == "rna":
            return self.make_rna_input, info["rna"]
        if component_type == "dna":
            return self.make_dna_input, info["dna"]
        raise ValueError(f"Unknown component type: {component_type!r}")

    def save_inputs_parallel(self, inputs, num_workers: int = 10):
        print(f"saving inputs parallelized: {len(inputs)}")
        with ThreadPoolExecutor(max_workers=num_workers) as ex:
            list(ex.map(self.save_single_input, inputs))

    def make_sbatch_file(self, **kwargs):
        return make_sbatch_file(
            self.get_model_command(),
            self.in_dir.parent,
            job_name=self.run_name,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # The workhorse: build inputs for a list of chunked component dicts
    # ------------------------------------------------------------------

    def make_inputs(
        self,
        chunked_components: list[dict],
        protein_csv: Path | None = None,
        protein_name_col: str = "protein_name",
        rna_csv: Path | None = None,
        rna_name_col: str = "rna_name",
        dna_csv: Path | None = None,
        dna_name_col: str = "dna_name",
        ligand_csv: Path | None = None,
        ligand_name_col: str = "lig_name",
        skip_existing_inputs: bool = False,
        skip_existing_structures: bool = False,
        use_tqdm: bool = False,
    ) -> None:
        with open(self.out_dir / "chunked_components.pkl", "wb") as fh:
            pickle.dump(chunked_components, fh)

        info = self.load_csvs(
            protein_csv, protein_name_col,
            rna_csv, rna_name_col,
            dna_csv, dna_name_col,
            ligand_csv, ligand_name_col,
        )

        inputs: list = []
        for chunk_idx, components in enumerate(chunked_components):
            chunk = str(chunk_idx)
            iterable = components.items()
            if use_tqdm:
                iterable = tqdm(components.items(), total=len(components), desc="Making inputs")

            for complex_name, (comp_list, extra) in iterable:
                if skip_existing_structures and self.check_done(complex_name, chunk):
                    input_path = self.get_input_path(complex_name, chunk)
                    if input_path.exists():
                        dest_dir = self.old_in_dir / chunk
                        dest_dir.mkdir(parents=True, exist_ok=True)
                        input_path.rename(dest_dir / input_path.name)
                    continue
                if skip_existing_inputs and self.get_input_path(complex_name, chunk).exists():
                    continue

                single_input_components: list = []
                curr_chains = 0
                for component_type, component_name, copies in comp_list:
                    chain_ids, curr_chains = self.next_chain_ids(copies, curr_chains)
                    input_maker, component_info = self._input_maker_factory(component_type, info)
                    single_input_components.extend(
                        input_maker(chain_ids, component_info[component_name])
                    )

                inputs.append(
                    (chunk, self.assemble_component_inputs(
                        single_input_components, complex_name, extra,
                    ))
                )

        print(f"Saving {len(inputs)} for {self.model_name}")
        if len(inputs) > 100:
            self.save_inputs_parallel(inputs)
        else:
            self.save_inputs(inputs)


def resolve_script_path(env_var: str, default_name: str) -> str:
    """Resolve the wrapper shell script path for a model.

    Precedence:
      1. Env var override (`COPAL_AF3_SCRIPT`, `COPAL_BOLTZ2_SCRIPT`, `COPAL_RF3_SCRIPT`).
      2. Packaged template at `copal/runner/scripts/<default_name>`.

    Templates are examples — most users will want to copy one, edit it to
    reference their local AF3/Boltz2/RF3 install, and point the env var
    at the edited copy. See `copal/runner/scripts/README.md`.
    """
    override = os.environ.get(env_var)
    if override:
        return str(Path(override).expanduser().resolve())
    return str(Path(__file__).parent / "scripts" / default_name)
