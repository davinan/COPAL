# Quickstart example

Two demo receptors × three demo AHL ligands (the cross-product gives 6
predictions, batched into one SLURM array task).

```
receptors.csv  →  receptor_id, seq
ligands.csv    →  ligand_id, smiles
run.sh         →  invokes `copal predict` with the above
```

Run:

```bash
cd examples/quickstart
./run.sh
```

This writes:

```
out/
├── scheduled_ensemble.csv      # one row per (receptor, ligand)
├── inputs/
│   ├── af3/    *.json
│   ├── boltz2/ *.yaml
│   └── rf3/    *.json
├── sarray/
│   ├── copal.cmds
│   └── copal.sbatch
└── matches/                     # populated after SLURM jobs run
```

After predictions land in `out/matches/`, rank them:

```bash
uv run copal analyze \
    --workflow-dir out/workflow \
    --output-root out/results \
    --labels   ../../data/benchmark_data_all_corrected.csv \
    --ligands  ../../data/Final_confirmed_ligands_58_annotated.csv \
    --basepath out/matches
```

(The two demo receptors are not in the labels CSV, so ranking-performance
metrics will be empty — `analyze` still produces the per-prediction
ranking and ensemble CSVs.)

For a real workflow on your own data: copy `receptors.csv` /
`ligands.csv`, replace the contents (any number of rows), and re-run
`run.sh` pointing at the new files.
