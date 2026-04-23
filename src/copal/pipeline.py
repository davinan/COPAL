"""COPAL inference pipeline.

Reads a receptor CSV and a ligand CSV, builds the cross-product of
predictions to run, batches them, and emits per-model inputs (AF3 JSON,
Boltz2 YAML, RF3 JSON) plus a SLURM array submission script using the
in-tree `copal.runner` module.
"""
from __future__ import annotations

import json
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from copal.io.schema import load_ligands_csv, load_receptors_csv
from copal.runner import EnsembleManager

logger = logging.getLogger(__name__)


def build_schedule(receptors_df: pd.DataFrame, ligands_df: pd.DataFrame,
                   batch_size: int, model_names: list[str]) -> pd.DataFrame:
    """Cross-product of receptors × ligands, sliced into batches."""
    records = []
    batch_count = 0
    for _, recv in receptors_df.iterrows():
        for _, lig in ligands_df.iterrows():
            target_name = f"{recv['receptor_id']}_{lig['ligand_id']}"
            records.append({
                "target_name": target_name,
                "protein_name": recv["receptor_id"],
                "ligand_name": lig["ligand_id"],
                "batch_id": batch_count // batch_size,
                **{f"is_done_{m}": False for m in model_names},
            })
            batch_count += 1
    return pd.DataFrame(records)


def enrich_done_status(
    df: pd.DataFrame,
    input_manager,
    model_names: list[str],
    done_lookup: dict | None = None,
) -> pd.DataFrame:
    done_lookup = done_lookup or {}

    def _enrich(row):
        out = row.to_dict()
        target_upper = row["target_name"].upper()
        status = input_manager.check_done(row["target_name"], row["batch_id"])
        for m, v in status.items():
            if target_upper in done_lookup and m in done_lookup.get(target_upper, []):
                v = True
            out[f"is_done_{m}"] = v
        return out

    with ThreadPoolExecutor() as ex:
        rows = list(tqdm(
            ex.map(_enrich, (row for _, row in df.iterrows())),
            total=len(df),
            desc="Checking completed predictions",
        ))
    return pd.DataFrame(rows)


def run_pipeline(
    *,
    receptors_csv: Path,
    ligands_csv: Path,
    out_dir: Path,
    model_names: list[str],
    batch_size: int = 20,
    run_name: str = "copal",
    chunk_size: int = 100,
    done_lookup_path: Path | None = None,
    redo_schedule: bool = False,
    skip_done_check: bool = False,
    executor: str = "slurm",
    workers: int = 1,
    dry_run: bool = False,
    submit: bool = False,
    sbatch_partition: str = "gpu",
    sbatch_gres: str = "gpu",
    sbatch_time: str = "10:00:00",
    sbatch_mem: str = "32G",
    sbatch_constraint: str = "[A100|H200|L40|A4000|A5000|A6000]",
    max_parallel: int = 500,
) -> None:
    """End-to-end: build schedule → check done → write inputs → dispatch runs.

    `executor="slurm"` writes a SLURM array script (default, like before).
    `executor="local"` runs the per-chunk shell commands on the current
    machine, optionally in parallel via ``workers > 1``.
    """
    if executor not in {"slurm", "local"}:
        raise ValueError(f"executor must be 'slurm' or 'local', got {executor!r}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    input_manager = EnsembleManager(
        models=model_names,
        run_name=run_name,
        working_dir=str(out_dir),
        chunk_size=chunk_size,
    )

    scheduled_csv = out_dir / "scheduled_ensemble.csv"
    done_lookup: dict = {}
    if done_lookup_path is not None:
        with open(done_lookup_path) as fh:
            done_lookup = json.load(fh)
        logger.info("Loaded done-lookup with %d entries", len(done_lookup))

    if scheduled_csv.exists() and not redo_schedule:
        schedule_df = pd.read_csv(scheduled_csv)
        logger.info("Reusing %s (%d rows)", scheduled_csv, len(schedule_df))
    else:
        receptors_df = load_receptors_csv(receptors_csv)
        ligands_df = load_ligands_csv(ligands_csv)
        logger.info("Loaded %d receptors and %d ligands", len(receptors_df), len(ligands_df))
        schedule_df = build_schedule(receptors_df, ligands_df, batch_size, model_names)
        schedule_df.to_csv(scheduled_csv, index=False)
        logger.info("Wrote %s (%d rows)", scheduled_csv, len(schedule_df))

    if not skip_done_check:
        schedule_df = enrich_done_status(
            schedule_df, input_manager, model_names,
            done_lookup=done_lookup,
        )
        schedule_df.to_csv(scheduled_csv, index=False)

    for sub in ("inputs", "sarray"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    n_chunks = int(schedule_df["batch_id"].max()) + 1
    chunked_components: list[dict] = [{} for _ in range(n_chunks)]
    for _, row in schedule_df.iterrows():
        chunked_components[int(row["batch_id"])][row["target_name"]] = (
            [
                ("protein", row["protein_name"], 1),
                ("ligand", row["ligand_name"], 1),
            ],
            None,
        )
    chunked_components = [c for c in chunked_components if c]

    model_list_override = {
        m: set(schedule_df.loc[~schedule_df[f"is_done_{m}"], "target_name"])
        for m in model_names
    }
    for m, targets in model_list_override.items():
        logger.info("%d predictions remaining for %s", len(targets), m)

    all_target_names = [n for chunk in chunked_components for n in chunk]

    input_manager.make_inputs(
        components=chunked_components,
        model_list_override=model_list_override,
        components_model_extras={
            "boltz2": {
                name: {"properties": [{"affinity": {"binder": "B"}}]}
                for name in all_target_names
            },
        },
        protein_csv=str(receptors_csv),
        protein_name_col="receptor_id",
        ligand_csv=str(ligands_csv),
        ligand_name_col="ligand_id",
        skip_existing_inputs=False,
        skip_existing_structures=True,
    )

    if executor == "slurm":
        input_manager.make_sbatch_file(
            partition=sbatch_partition,
            gres=sbatch_gres,
            time=sbatch_time,
            mem=sbatch_mem,
            max_parallel=max_parallel,
            constraint=sbatch_constraint,
            submit_for_me=submit,
        )
    else:  # local
        failures = input_manager.run_locally(workers=workers, dry_run=dry_run)
        if failures:
            logger.warning("%d task(s) failed during local execution", failures)


