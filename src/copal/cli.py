"""Top-level COPAL command-line entry point.

Subcommands:
    predict       — schedule + emit AF3/Boltz2/RF3 inputs and SLURM script
    analyze       — collect predictions, postprocess, generate analysis CSVs
    seqid-filter  — sequence-identity / leakage filter (BLAST + global align)
    organize      — move raw model outputs into a hierarchical store
"""
from __future__ import annotations

from pathlib import Path

import typer

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="COPAL — co-folding ranking of protein-ligand pairs.",
)


@app.command()
def predict(
    receptors: Path = typer.Option(..., "--receptors", "-r", exists=True, dir_okay=False,
                                   help="CSV with at least columns: receptor_id, seq."),
    ligands: Path = typer.Option(..., "--ligands", "-l", exists=True, dir_okay=False,
                                 help="CSV with at least columns: ligand_id, smiles."),
    out_dir: Path = typer.Option(..., "--out-dir", "-o", file_okay=False,
                                 help="Working directory for inputs/, sarray/, matches/."),
    models: str = typer.Option("af3,boltz2,rf3", "--models",
                               help="Comma-separated models to run."),
    batch_size: int = typer.Option(20, "--batch-size",
                                   help="Predictions per chunk (= per SLURM task / per local command)."),
    run_name: str = typer.Option("copal", "--run-name"),
    done_lookup: Path | None = typer.Option(None, "--done-lookup", exists=True, dir_okay=False),
    redo_schedule: bool = typer.Option(False, "--redo-schedule"),
    skip_done_check: bool = typer.Option(False, "--skip-done-check"),
    executor: str = typer.Option("slurm", "--executor",
                                 help="How to run the model commands: 'slurm' or 'local'."),
    workers: int = typer.Option(1, "--workers",
                                help="Parallel workers when --executor local (serial if 1)."),
    dry_run: bool = typer.Option(False, "--dry-run",
                                 help="Write inputs + .cmds file but skip execution."),
    submit: bool = typer.Option(False, "--submit",
                                help="Auto-submit the sbatch (slurm executor only)."),
) -> None:
    """Schedule + emit per-model inputs, then dispatch via SLURM or run locally."""
    from copal.pipeline import run_pipeline

    run_pipeline(
        receptors_csv=receptors,
        ligands_csv=ligands,
        out_dir=out_dir,
        model_names=[m.strip() for m in models.split(",") if m.strip()],
        batch_size=batch_size,
        run_name=run_name,
        done_lookup_path=done_lookup,
        redo_schedule=redo_schedule,
        skip_done_check=skip_done_check,
        executor=executor,
        workers=workers,
        dry_run=dry_run,
        submit=submit,
    )


@app.command(
    context_settings={
        "allow_extra_args": True,
        "ignore_unknown_options": True,
    },
)
def analyze(ctx: typer.Context) -> None:
    """Run the Hydra-driven single-set analysis.

    Defaults to ``set=refAHL_all_conf_no_leakage``. Override any config group or
    leaf via standard Hydra syntax, for example:

        copal analyze
        copal analyze set=refAHL_conf1
        copal analyze analysis.aggregation=rrf paths.output_root=/tmp/sweep
        copal analyze 'set.filter.confidences=[1,2]'
    """
    from copal.analysis.run import run_with_overrides

    run_with_overrides(ctx.args)


@app.command("seqid-filter")
def seqid_filter(
    queries: Path = typer.Option(..., "--queries", exists=True, dir_okay=False),
    references: Path = typer.Option(..., "--references", exists=True, dir_okay=False),
    out_table: Path = typer.Option(..., "--out-table"),
    out_leakage: Path | None = typer.Option(None, "--out-leakage"),
    threshold: float = typer.Option(0.40, "--threshold"),
    blast_bin: Path | None = typer.Option(None, "--blast-bin", file_okay=False),
) -> None:
    """BLAST + global-alignment leakage filter between query and reference sequences."""
    import sys

    from copal.analysis.compute_seqid import main as seqid_main

    cli_args = [
        "--queries", str(queries),
        "--references", str(references),
        "--out-table", str(out_table),
        "--threshold", str(threshold),
    ]
    if out_leakage is not None:
        cli_args += ["--out-leakage", str(out_leakage)]
    if blast_bin is not None:
        cli_args += ["--blast-bin", str(blast_bin)]

    saved = sys.argv
    try:
        sys.argv = ["copal-seqid-filter", *cli_args]
        seqid_main()
    finally:
        sys.argv = saved


@app.command()
def organize(
    pred_dir: Path = typer.Option(..., "--pred-dir", exists=True, file_okay=False),
    dst_dir: Path = typer.Option(..., "--dst-dir", file_okay=False),
    log_path: Path | None = typer.Option(None, "--log-path"),
    models: str = typer.Option("boltz2,af3,rf3", "--models"),
    stage_only: bool = typer.Option(False, "--stage-only"),
) -> None:
    """Move raw model outputs into a hierarchical predictions store."""
    import sys

    from copal.organizer import main as organizer_main

    cli_args = [
        "--pred-dir", str(pred_dir),
        "--dst-dir", str(dst_dir),
        "--models", models,
    ]
    if log_path is not None:
        cli_args += ["--log-path", str(log_path)]
    if stage_only:
        cli_args.append("--stage-only")

    saved = sys.argv
    try:
        sys.argv = ["copal-organize", *cli_args]
        organizer_main()
    finally:
        sys.argv = saved


if __name__ == "__main__":
    app()
