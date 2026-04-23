"""BLASTP query functions for protein sequence search."""

from pathlib import Path
import tempfile
from io import StringIO
import subprocess
import pandas as pd


def write_query_fasta(seq: str, name: str = "query") -> Path:
    """Write a query sequence to a temporary FASTA file.
    
    Args:
        seq: Protein sequence string
        name: Sequence name/identifier for FASTA header
        
    Returns:
        Path to the temporary FASTA file
    """
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False)
    tmp.write(f">{name}\n{seq}\n")
    tmp.close()
    return Path(tmp.name)


def blastp_query(
    query_seq: str,
    db_path: str,
    max_hits: int = 50,
    evalue: float = 1e-5,
) -> pd.DataFrame:
    """Run BLASTP query against a protein database.
    
    Args:
        query_seq: Query protein sequence
        db_path: Path to BLASTP database (without extensions)
        max_hits: Maximum number of hits to return
        evalue: E-value threshold
        
    Returns:
        DataFrame with BLASTP results columns:
        qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore
    """
    query_faa = write_query_fasta(query_seq)

    out_cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]

    cmd = [
        "blastp",
        "-query", str(query_faa),
        "-db", db_path,
        "-outfmt", "6 " + " ".join(out_cols),
        "-max_target_seqs", str(max_hits),
        "-evalue", str(evalue),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )

    if not result.stdout.strip():
        return pd.DataFrame(columns=out_cols)

    df = pd.read_csv(
        StringIO(result.stdout),
        sep="\t",
        names=out_cols,
    )

    return df
