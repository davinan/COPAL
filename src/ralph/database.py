"""Database loading and setup utilities."""

import os
import pandas as pd
from typing import Optional

try:
    from rdkit.Chem import PandasTools
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    PandasTools = None


def load_ralph_database(
    database_path: str = "/content/ralph_database",
    google_drive_file_id: str = "1_1ab-A_0P9D2X0XZS-eeP0uV0HFGeCjk",
    download: bool = True,
) -> str:
    """Download and setup RALPH database.
    
    Args:
        database_path: Path where database should be located
        google_drive_file_id: Google Drive file ID for database download
        download: Whether to download if database doesn't exist
        
    Returns:
        Path to the database directory
    """
    if not os.path.exists(database_path) and download:
        # This will be executed in notebook context
        # The actual download code should be in the notebook
        print(f"Database not found at {database_path}")
        print("Please download the database using the notebook setup cell.")
    
    return database_path


def load_scores(database_path: str = "/content/ralph_database") -> pd.DataFrame:
    """Load RALPH scores DataFrame.
    
    Args:
        database_path: Path to RALPH database directory
        
    Returns:
        DataFrame with RALPH scores
    """
    scores_path = os.path.join(database_path, "ralph_scores.csv")
    if not os.path.exists(scores_path):
        raise FileNotFoundError(f"Scores file not found at {scores_path}")
    
    df = pd.read_csv(scores_path)
    print("Loaded Scores DataFrame")
    return df


def load_ligands(database_path: str = "/content/ralph_database") -> pd.DataFrame:
    """Load ligands DataFrame with molecular structures.
    
    Args:
        database_path: Path to RALPH database directory
        
    Returns:
        DataFrame with ligands, including 'ROMol' column with RDKit molecules
        
    Raises:
        ImportError: If RDKit is not installed. Install with: pip install ralph[rdkit]
                     or: pip install rdkit-pypi
    """
    if not RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required for load_ligands() but is not installed. "
            "Install it with one of the following:\n"
            "  - pip install ralph[rdkit]\n"
            "  - pip install rdkit-pypi\n"
            "  - conda install -c conda-forge rdkit"
        )
    
    ligands_path = os.path.join(database_path, "ligands_58.csv")
    if not os.path.exists(ligands_path):
        raise FileNotFoundError(f"Ligands file not found at {ligands_path}")
    
    ligand_df = pd.read_csv(ligands_path)
    PandasTools.AddMoleculeColumnToFrame(ligand_df, smilesCol='smi', molCol='ROMol')
    print("Added 'ROMol' column to ligand_df with loaded molecules.")
    return ligand_df
