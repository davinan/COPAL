"""Hydra-driven entry point for copal analysis on a single benchmark set."""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Sequence

import hydra
from hydra import compose, initialize_config_module
from omegaconf import DictConfig, OmegaConf

from copal.analysis import pipeline
from copal.analysis.sets import from_omega, load_set


log = logging.getLogger(__name__)


# Custom resolver: ${cwd:} returns the current working directory. Works in both
# @hydra.main and compose() contexts, unlike Hydra's built-in ${hydra:runtime.cwd}
# which is only available inside @hydra.main.
if not OmegaConf.has_resolver("cwd"):
    OmegaConf.register_new_resolver("cwd", lambda: os.getcwd())


def _resolve_set(cfg: DictConfig):
    """Compose a BenchmarkSet from the Hydra ``set`` group, allowing null paths."""
    spec = from_omega(cfg.set)
    labels_path = cfg.set.get("labels_path")
    ligands_path = cfg.set.get("ligands_path")
    return load_set(
        spec,
        Path(labels_path) if labels_path else None,
        Path(ligands_path) if ligands_path else None,
    )


def _execute(cfg: DictConfig) -> None:
    """Build the set and run the single-set analysis."""
    log.info("Resolved config:\n%s", OmegaConf.to_yaml(cfg, resolve=True))

    bset = _resolve_set(cfg)
    log.info(
        "Loaded set %r: n_luxrs=%d n_rows=%d",
        bset.name, bset.n_luxrs(), len(bset.labels),
    )

    skip_postprocess = bool(cfg.analysis.skip.get("postprocess", True))
    if not skip_postprocess:
        workflow_dir = Path(cfg.paths.workflow_dir)
        workflow_dir.mkdir(parents=True, exist_ok=True)
        input_csv_value = cfg.paths.get("input_csv")
        if input_csv_value is None:
            raise ValueError(
                "analysis.skip.postprocess=false requires paths.input_csv to "
                "point at a pre-collected predictions CSV."
            )
        pipeline.run_postprocess(
            Path(input_csv_value),
            workflow_dir,
            Path(cfg.set.labels_path),
            skip_ablation=False,
            aggregation=str(cfg.analysis.aggregation),
            rrf_k=int(cfg.analysis.rrf_k),
        )

    chem_sim_path = cfg.analysis.get("chem_sim_matrix")
    pipeline.run_analysis_for_set(
        bset=bset,
        workflow_dir=Path(cfg.paths.workflow_dir),
        output_root=Path(cfg.paths.output_root),
        aggregation=str(cfg.analysis.aggregation),
        rrf_k=int(cfg.analysis.rrf_k),
        per_model_split=tuple(cfg.analysis.per_model_split),
        skip_ablation=bool(cfg.analysis.skip.get("ablation", False)),
        chem_sim_matrix_path=Path(chem_sim_path) if chem_sim_path else None,
    )


@hydra.main(version_base=None, config_path="../config", config_name="config")
def hydra_entry(cfg: DictConfig) -> None:
    """Hydra `@main`-decorated entry. Driven by sys.argv overrides."""
    _execute(cfg)


def run_with_overrides(overrides: Sequence[str] | None = None) -> None:
    """Programmatic entry: compose the config with overrides, then run.

    Used by the `copal analyze` typer subcommand to forward user CLI overrides
    without going through `@hydra.main`'s sys.argv parsing.
    """
    overrides = list(overrides or [])
    with initialize_config_module(version_base=None, config_module="copal.config"):
        cfg = compose(config_name="config", overrides=overrides)
    _execute(cfg)


def _execute_full_ranking(cfg: DictConfig) -> None:
    """Build the set and emit stylized full_ranking .xlsx files."""
    log.info("Resolved config:\n%s", OmegaConf.to_yaml(cfg, resolve=True))

    bset = _resolve_set(cfg)
    log.info(
        "Loaded set %r: n_luxrs=%d n_rows=%d",
        bset.name, bset.n_luxrs(), len(bset.labels),
    )

    pipeline.run_full_ranking_for_set(
        bset=bset,
        workflow_dir=Path(cfg.paths.workflow_dir),
        output_root=Path(cfg.paths.output_root),
        aggregation=str(cfg.analysis.aggregation),
        rrf_k=int(cfg.analysis.rrf_k),
        per_model_split=tuple(cfg.analysis.per_model_split),
    )


def run_full_ranking_with_overrides(overrides: Sequence[str] | None = None) -> None:
    """Programmatic entry: compose the config with overrides, then emit .xlsx files.

    Used by the `copal full-ranking` typer subcommand to forward user CLI
    overrides without going through `@hydra.main`'s sys.argv parsing.
    """
    overrides = list(overrides or [])
    with initialize_config_module(version_base=None, config_module="copal.config"):
        cfg = compose(config_name="config", overrides=overrides)
    _execute_full_ranking(cfg)


def _execute_sim_matrix(cfg: DictConfig) -> None:
    """Precompute a pairwise ligand-similarity matrix for the configured pool."""
    import pandas as pd

    from copal.analysis.chem_sim import write_pairwise_matrix

    log.info("Resolved config:\n%s", OmegaConf.to_yaml(cfg, resolve=True))

    ligands_path = cfg.set.get("ligands_path")
    if not ligands_path:
        raise ValueError(
            "sim-matrix requires a set with ligands_path set "
            "(e.g. set=refAHL_all_conf_no_leakage); 'inference' has no ligand pool."
        )
    ligands_df = pd.read_csv(Path(ligands_path))
    method = str(cfg.sim_matrix.method)
    out_path = Path(cfg.sim_matrix.out)
    log.info("Computing pairwise sim matrix: method=%r ligands=%s out=%s",
             method, ligands_path, out_path)
    df = write_pairwise_matrix(ligands_df, method, out_path)
    log.info("Wrote %s (shape=%s, diagonal=1.0, symmetric)", out_path, df.shape)
    print(f"Wrote {out_path} ({df.shape[0]}x{df.shape[1]} matrix, method={method})")


def run_sim_matrix_with_overrides(overrides: Sequence[str] | None = None) -> None:
    """Programmatic entry: compose the config with overrides, then write matrix."""
    overrides = list(overrides or [])
    with initialize_config_module(version_base=None, config_module="copal.config"):
        cfg = compose(config_name="config", overrides=overrides)
    _execute_sim_matrix(cfg)


if __name__ == "__main__":
    hydra_entry()
