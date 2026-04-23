"""Run per-chunk model commands locally instead of via SLURM.

Intended for users who want to try COPAL on a workstation with a single
GPU, or who don't have a SLURM cluster at all.
"""
from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def write_cmds_file(cmds: list[str], cmds_dir: Path, job_name: str = "copal") -> Path:
    """Write ``{cmds_dir}/{job_name}.cmds`` (one shell command per line).

    The `.cmds` file is preserved in both SLURM and local modes so a
    failed run can be resumed or inspected independently.
    """
    cmds_dir = Path(cmds_dir)
    cmds_dir.mkdir(parents=True, exist_ok=True)
    path = cmds_dir / f"{job_name}.cmds"
    payload = "\n".join(cmds) + ("\n" if cmds else "")
    path.write_text(payload)
    return path


def _run_one(idx: int, cmd: str, log_dir: Path | None) -> int:
    tag = f"task-{idx + 1:03d}"
    print(f"[{tag}] $ {cmd}", flush=True)
    stdout = stderr = None
    if log_dir is not None:
        stdout = open(log_dir / f"{tag}.out", "w")
        stderr = open(log_dir / f"{tag}.err", "w")
    try:
        rc = subprocess.call(cmd, shell=True, stdout=stdout, stderr=stderr)
    finally:
        if stdout is not None:
            stdout.close()
        if stderr is not None:
            stderr.close()
    status = "OK" if rc == 0 else f"FAIL (rc={rc})"
    print(f"[{tag}] {status}", flush=True)
    return rc


def run_commands_locally(
    cmds: list[str],
    workers: int = 1,
    log_dir: Path | None = None,
) -> int:
    """Execute each command on the local shell. Returns the number of failures.

    With ``workers > 1``, commands run concurrently via a thread pool.
    Output either streams to the parent terminal (default) or to
    per-task ``.out`` / ``.err`` files under ``log_dir``.
    """
    if not cmds:
        print("No commands to run.")
        return 0
    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)

    failures = 0
    if workers <= 1:
        for i, cmd in enumerate(cmds):
            if _run_one(i, cmd, log_dir) != 0:
                failures += 1
    else:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            rcs = list(ex.map(lambda ic: _run_one(*ic, log_dir), enumerate(cmds)))
        failures = sum(1 for rc in rcs if rc != 0)

    total = len(cmds)
    print(f"\nLocal execution: {total - failures}/{total} tasks succeeded.")
    return failures
