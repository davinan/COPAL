"""Build a SLURM array submission (.cmds + .sbatch + submit helper)."""
from __future__ import annotations

import subprocess
from pathlib import Path


def make_sbatch_file(
    cmds,
    cmds_dir: Path,
    job_name: str = "copal",
    partition: str = "gpu",
    time: str = "10:00:00",
    gres: str = "gpu:a4000:1",
    mem: str = "64G",
    max_parallel: int = 100,
    constraint: str | None = "[A100|L40|H200|A6000|A4000]",
    submit_for_me: bool = False,
):
    """Emit ``{cmds_dir}/{job_name}.cmds`` + ``.sbatch`` + ``submit_sarray.sh``.

    One shell command per line in the .cmds file; the .sbatch uses
    ``SLURM_ARRAY_TASK_ID`` to dispatch a single line per array task.
    """
    cmds_dir = Path(cmds_dir)
    cmds_dir.mkdir(parents=True, exist_ok=True)

    sbatch = [
        "#!/bin/bash",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={cmds_dir/'logs'}/slurm-%X-%a.out",
        f"#SBATCH --error={cmds_dir/'logs'}/slurm-%X-%a.err",
        f"#SBATCH --time={time}",
    ]
    if "bf" not in partition:
        sbatch.append(f"#SBATCH --gres={gres}")
    sbatch.append("#SBATCH --cpus-per-task=2")
    sbatch.append(f"#SBATCH --mem={mem}")
    if constraint is not None:
        sbatch.append(f"#SBATCH --constraint={constraint}")

    cmds_file = cmds_dir / f"{job_name}.cmds"
    sarray_file = cmds_dir / f"{job_name}.sbatch"

    sbatch.append(f'CMD=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {cmds_file})')
    sbatch.append('echo "${CMD}" | bash')

    cmds_file.write_text("\n".join(cmds))
    sarray_file.write_text("\n".join(sbatch))

    submit_helper = cmds_dir / "submit_sarray.sh"
    submit_cmd = f"sbatch -a 1-$(cat {cmds_file}| wc -l)%{max_parallel} {sarray_file}"
    submit_helper.write_text(submit_cmd + "\n")

    print("*" * 100)
    print(f"submittable file written to {sarray_file}")
    if submit_for_me:
        print(f"Submitting with:\n{submit_cmd}")
        subprocess.run(submit_cmd, shell=True, check=False)
    else:
        print(f"Submit with:\n{submit_cmd}")
        print(f"or with the submit_sarray.sh script:\n{submit_helper}")
    print("*" * 100)
