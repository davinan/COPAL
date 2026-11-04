"""Database loading utilities for the RALPH-BLAST inference database.

The shipped database directory (see baselines/ralph_blast/README.md) contains
compressed sources (``*.csv.xz`` / ``*.faa.xz``) plus ``ligands_58.csv``.
``setup_db.sh`` decompresses them to plain ``.csv``/``.faa`` and runs
``makeblastdb``. The loaders below prefer the materialized ``.csv`` and fall
back to reading the shipped ``.csv.xz`` directly (pandas auto-detects ``.xz``),
so they work whether or not setup has run.
"""

import os
import warnings

import pandas as pd

try:
    from rdkit.Chem import PandasTools
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    PandasTools = None


def _resolve(database_path: str, stem: str, ext: str) -> str:
    """Return <database_path>/<stem>.<ext> if present, else the .xz sibling.

    Raises FileNotFoundError if neither exists.
    """
    plain = os.path.join(database_path, f"{stem}.{ext}")
    if os.path.exists(plain):
        return plain
    packed = plain + ".xz"
    if os.path.exists(packed):
        return packed
    raise FileNotFoundError(
        f"Neither {plain} nor {packed} found. Build the database "
        f"(see baselines/ralph_blast/README.md) or run setup_db.sh."
    )


def load_ralph_database(database_path: str = "ralph_database", **_ignored) -> str:
    """Validate and return the local RALPH database directory.

    .. deprecated::
        This used to download the database from Google Drive. The database now
        ships in-repo (compressed) under ``baselines/ralph_blast/ralph_database``;
        run ``baselines/ralph_blast/setup_db.sh`` to materialize it. This shim
        only validates the local path and is kept for backwards compatibility.
    """
    if _ignored:
        warnings.warn(
            "load_ralph_database no longer downloads; the database ships in-repo. "
            "Pass a local database_path and run setup_db.sh. Extra args ignored.",
            DeprecationWarning,
            stacklevel=2,
        )
    if not os.path.isdir(database_path):
        raise FileNotFoundError(
            f"RALPH database directory not found at {database_path}. "
            f"See baselines/ralph_blast/README.md (run setup_db.sh)."
        )
    return database_path


def load_scores(database_path: str = "ralph_database") -> pd.DataFrame:
    """Load the wide RALPH+RefAHL rankings table.

    Columns: ``luxr_name``, ``source`` (RALPH|RefAHL), ``organism``,
    ``true_ahl_id``, ``mpKT``, ``sequence``, and ``{model}_1..58`` (the AHL id at
    each rank) for ``model`` in ensemble/af3/rf3/boltz2/boltz2affinity.
    """
    df = pd.read_csv(_resolve(database_path, "ralph_scores", "csv"), low_memory=False)
    print(f"Loaded scores: {len(df)} receptors")
    return df


def load_model_scores(database_path: str = "ralph_database") -> pd.DataFrame:
    """Load per-pair COPAL scores, indexed by ``(luxr_name, ahl)`` for O(1) lookup.

    The only score column is ``copal`` (ensemble sum-of-ranks).
    Ignore individual-model metrics in older materialized files.
    Return an empty frame gracefully if the score file is absent.
    """
    try:
        path = _resolve(database_path, "ralph_model_scores", "csv")
    except FileNotFoundError:
        return pd.DataFrame()
    df = pd.read_csv(path, usecols=["luxr_name", "ahl", "copal"])
    return df.set_index(["luxr_name", "ahl"]).sort_index()


def load_ligands(database_path: str = "ralph_database") -> pd.DataFrame:
    """Load the 58 AHL ligands with RDKit molecules for visualization.

    Returns a DataFrame with ``AHLid``, ``abbreviation``, ``smi`` and a ``ROMol``
    column. Requires RDKit.
    """
    if not RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required for load_ligands() but is not installed. "
            "Install it with:  pip install rdkit"
        )
    ligands_path = os.path.join(database_path, "ligands_58.csv")
    if not os.path.exists(ligands_path):
        raise FileNotFoundError(f"Ligands file not found at {ligands_path}")

    ligand_df = pd.read_csv(ligands_path)
    PandasTools.AddMoleculeColumnToFrame(ligand_df, smilesCol="smi", molCol="ROMol")
    return ligand_df
