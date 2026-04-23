"""COPAL — co-folding ranking of protein-ligand pairs.

Public submodules:
    copal.pipeline    — schedule + submit AF3/Boltz2/RF3 inference
    copal.runner      — in-tree AF3/Boltz2/RF3 managers + SLURM submission
    copal.analysis    — collect predictions, compute rankings, analyse results
    copal.io          — receptor/ligand CSV schema and validation
    copal.cli         — `copal` command-line entry point
"""

__version__ = "0.1.0"
