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


def _execute(cfg: DictConfig) -> None:
    """Build the set and run the single-set analysis."""
    log.info("Resolved config:\n%s", OmegaConf.to_yaml(cfg, resolve=True))

    spec = from_omega(cfg.set)
    bset = load_set(
        spec,
        Path(cfg.set.labels_path),
        Path(cfg.set.ligands_path),
    )
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

    pipeline.run_analysis_for_set(
        bset=bset,
        workflow_dir=Path(cfg.paths.workflow_dir),
        output_root=Path(cfg.paths.output_root),
        aggregation=str(cfg.analysis.aggregation),
        rrf_k=int(cfg.analysis.rrf_k),
        per_model_split=tuple(cfg.analysis.per_model_split),
        skip_ablation=bool(cfg.analysis.skip.get("ablation", False)),
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


if __name__ == "__main__":
    hydra_entry()
