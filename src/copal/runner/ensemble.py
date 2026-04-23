"""Orchestrate input generation + SLURM submission across one or more models."""
from __future__ import annotations

import pickle
from pathlib import Path

from copal.runner.af3 import AF3Manager
from copal.runner.boltz2 import Boltz2Manager
from copal.runner.local import run_commands_locally, write_cmds_file
from copal.runner.rf3 import RF3Manager
from copal.runner.sbatch import make_sbatch_file

_MODEL_FACTORY = {
    "af3": AF3Manager,
    "boltz2": Boltz2Manager,
    "rf3": RF3Manager,
}


class EnsembleManager:
    """Group one or more model managers under a single working directory."""

    def __init__(self, models: list[str], run_name: str, working_dir: Path | str, chunk_size: int = 1):
        working_dir = Path(working_dir)
        self.out_dir = working_dir / "outs"
        self.in_dir = working_dir / "inputs"
        self.run_name = run_name

        unknown = [m for m in models if m not in _MODEL_FACTORY]
        if unknown:
            raise ValueError(
                f"Unknown model(s) {unknown}; expected subset of {sorted(_MODEL_FACTORY)}"
            )
        unique_models = list(dict.fromkeys(models))  # dedupe, preserve order
        self.model_ids = unique_models
        self.models = [
            _MODEL_FACTORY[m](run_name=run_name, out_dir=self.out_dir, in_dir=self.in_dir)
            for m in unique_models
        ]

    # ------------------------------------------------------------------

    def make_inputs(
        self,
        components,
        components_model_extras: dict | None = None,
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
        model_list_override: dict | None = None,
    ) -> None:
        components_model_extras = components_model_extras or {}
        model_list_override = model_list_override or {}

        with open(self.out_dir / "components.pkl", "wb") as fh:
            pickle.dump(components, fh)

        for model_id, model in zip(self.model_ids, self.models):
            print(f"making inputs for {model_id}")
            list_override = model_list_override.get(model_id)
            extras = components_model_extras.get(model_id)

            if list_override:
                model_components = [
                    {k: v for k, v in c.items() if k in list_override} for c in components
                ]
            else:
                model_components = components

            if extras:
                model_components = [
                    {k: (v[0], extras.get(k)) for k, v in c.items()} for c in model_components
                ]

            print(f"list_override: {len(list_override) if list_override else 0}")
            print(f"model_components: {len(model_components)}")
            model.make_inputs(
                model_components,
                protein_csv=protein_csv,
                protein_name_col=protein_name_col,
                rna_csv=rna_csv,
                rna_name_col=rna_name_col,
                dna_csv=dna_csv,
                dna_name_col=dna_name_col,
                ligand_csv=ligand_csv,
                ligand_name_col=ligand_name_col,
                skip_existing_inputs=skip_existing_inputs,
                skip_existing_structures=skip_existing_structures,
            )

    # ------------------------------------------------------------------

    def check_done(self, name: str, chunk) -> dict[str, bool]:
        return {
            mid: model.check_done(name, str(chunk))
            for mid, model in zip(self.model_ids, self.models)
        }

    def get_all_models_commands(self) -> list[str]:
        cmds: list = []
        for model in self.models:
            cmds.extend(model.get_model_command())
        # Preserve order, drop duplicates
        return list(dict.fromkeys(cmds))

    def make_sbatch_file(
        self,
        partition: str = "gpu",
        time: str = "10:00:00",
        gres: str = "gpu:a4000:1",
        mem: str = "32G",
        max_parallel: int = 100,
        constraint: str | None = None,
        submit_for_me: bool = False,
    ) -> None:
        make_sbatch_file(
            self.get_all_models_commands(),
            self.in_dir.parent,
            partition=partition,
            job_name=self.run_name,
            time=time,
            gres=gres,
            mem=mem,
            max_parallel=max_parallel,
            constraint=constraint,
            submit_for_me=submit_for_me,
        )

    def run_locally(
        self,
        workers: int = 1,
        log_dir: Path | str | None = None,
        dry_run: bool = False,
    ) -> int:
        """Execute model commands on the local machine (no SLURM).

        Writes the ``.cmds`` file alongside SLURM runs would, then runs
        each command via the local shell. Returns the number of failures.
        With ``dry_run=True``, only the ``.cmds`` file is written.
        """
        cmds = self.get_all_models_commands()
        cmds_dir = self.in_dir.parent
        cmds_file = write_cmds_file(cmds, cmds_dir, job_name=self.run_name)
        print("*" * 100)
        print(f"cmds file written to {cmds_file} ({len(cmds)} commands)")
        if dry_run:
            print("dry_run=True; skipping execution.")
            print("*" * 100)
            return 0
        resolved_log_dir = Path(log_dir) if log_dir is not None else (cmds_dir / "logs")
        print(f"Running locally with {workers} worker(s); logs: {resolved_log_dir}")
        print("*" * 100)
        return run_commands_locally(cmds, workers=workers, log_dir=resolved_log_dir)
