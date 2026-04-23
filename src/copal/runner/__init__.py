"""Vendored AF3 / Boltz2 / RF3 orchestration (formerly `model_runner`).

Public API:
    EnsembleManager  — orchestrates a multi-model prediction ensemble
    AF3Manager, Boltz2Manager, RF3Manager — per-model input managers
    make_sbatch_file — build a SLURM array submission
"""
from copal.runner.af3 import AF3Manager
from copal.runner.boltz2 import Boltz2Manager
from copal.runner.ensemble import EnsembleManager
from copal.runner.local import run_commands_locally, write_cmds_file
from copal.runner.rf3 import RF3Manager
from copal.runner.sbatch import make_sbatch_file

__all__ = [
    "EnsembleManager",
    "AF3Manager",
    "Boltz2Manager",
    "RF3Manager",
    "make_sbatch_file",
    "run_commands_locally",
    "write_cmds_file",
]
