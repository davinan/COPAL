"""Analysis pipeline for COPAL benchmark.

Modules:
    pipeline       — three-phase collect/postprocess/generate (legacy argparse CLI)
    sets           — BenchmarkSet abstraction; filter helpers
    run            — Hydra-driven single-set entry point
    compute_seqid  — sequence-identity leakage filter (BLAST + global align)
"""

from copal.analysis.pipeline import run_analysis_for_set
from copal.analysis.sets import BenchmarkSet, BenchmarkSetSpec, load_set

__all__ = [
    "BenchmarkSet",
    "BenchmarkSetSpec",
    "load_set",
    "run_analysis_for_set",
]
