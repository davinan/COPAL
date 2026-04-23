"""Sequence alignment functions using biotite."""

import subprocess
from pathlib import Path
import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.align as align
import biotite.sequence.graphics as graphics


def extract_sequence_from_blastdb(db_path: str, seq_id: str) -> str:
    """Extract a sequence from a BLAST database using blastdbcmd.
    
    Args:
        db_path: Path to BLAST database (without extensions)
        seq_id: Sequence identifier to extract
        
    Returns:
        Protein sequence string
    """
    cmd = [
        "blastdbcmd",
        "-db", db_path,
        "-entry", seq_id,
        "-outfmt", "%s"
    ]
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )
    
    if not result.stdout.strip():
        raise ValueError(f"Sequence {seq_id} not found in database {db_path}")
    
    # Remove any whitespace/newlines
    sequence = result.stdout.strip().replace("\n", "").replace(" ", "")
    return sequence


def plot_sequence_alignment(
    query_seq: str,
    target_seq: str,
    query_name: str = "query",
    target_name: str = "target",
) -> plt.Figure:
    """Create a pairwise sequence alignment plot using biotite.
    
    Args:
        query_seq: Query protein sequence
        target_seq: Target protein sequence
        query_name: Name/label for query sequence
        target_name: Name/label for target sequence
        
    Returns:
        Matplotlib figure with alignment plot
    """
    # Convert sequences to ProteinSequence objects
    sA = seq.ProteinSequence(query_seq)
    sB = seq.ProteinSequence(target_seq)
    
    # Get BLOSUM62 matrix
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    
    # Perform pairwise sequence alignment with affine gap penalty
    # Terminal gaps are not penalized (matching reference notebook)
    alignments = align.align_optimal(
        sA, sB, matrix, gap_penalty=(-10, -1), terminal_penalty=False
    )
    
    if len(alignments) == 0:
        raise ValueError("No alignment found between sequences")
    
    # Draw first alignment
    # The color intensity indicates the similarity
    fig = plt.figure(figsize=(8.0, 2.5))
    ax = fig.add_subplot(111)
    graphics.plot_alignment_similarity_based(
        ax,
        alignments[0],
        matrix=matrix,
        labels=[query_name, target_name],
        show_numbers=True,
        show_line_position=True,
    )
    fig.tight_layout()
    
    return fig
