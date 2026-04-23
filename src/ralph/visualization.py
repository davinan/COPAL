"""Visualization functions for ranking results."""

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from IPython.display import display, HTML
try:
    import ipywidgets as widgets
    from ipywidgets import interact, interactive, fixed, HBox, VBox, Output
    IPYWIDGETS_AVAILABLE = True
except ImportError:
    IPYWIDGETS_AVAILABLE = False
    widgets = None


def clustered_rank_plot(df_plot: pd.DataFrame, gap=2.0, highlight_labels=(), top_k=10):
    """Create a clustered bar plot of rankings.
    
    Args:
        df_plot: DataFrame with columns:
            - label: x tick label
            - cluster: grouping for clusters
            - best_rank: ranking value (smaller is better)
        gap: Gap size between clusters
        highlight_labels: Labels to highlight in red
        top_k: Number of top ranks to display
    """
    df_plot = df_plot.copy()

    max_rank = int(np.nanmax(df_plot["best_rank"]))
    df_plot["height"] = (max_rank + 1) - df_plot["best_rank"]
    df_plot["height"] = df_plot["height"].fillna(0)

    # Order within clusters: by length then unsaturation (so it looks tidy)
    # We extract numeric length + unsat count from the label.
    def parse_len_unsat(lbl):
        m = re.search(r"C(\d+)(?::(\d+))?$", lbl)
        if not m:
            return (999, 999)
        L = int(m.group(1))
        U = int(m.group(2)) if m.group(2) else 0
        return (L, U)

    df_plot["len"], df_plot["unsat"] = zip(*df_plot["label"].map(parse_len_unsat))
    df_plot = df_plot.sort_values(["cluster", "len", "unsat", "label"]).reset_index(drop=True)

    # Build x positions with gaps
    x = 0.0
    xs = []
    labels = []
    colors = []
    for clust, g in df_plot.groupby("cluster", sort=False):
        for _, row in g.iterrows():
            xs.append(x)
            labels.append(row["label"])
            colors.append("red" if row["label"] in set(highlight_labels) else "black")
            x += 1.0
        x += gap

    df_plot["x"] = xs

    # Plot
    fig, ax = plt.subplots(figsize=(14, 3))
    ax.bar(df_plot["x"], df_plot["height"], width=0.85, color=colors)

    # y-axis: show rank numbers even though heights are inverted
    rank_labels = np.arange(1, max_rank + 1)
    tick_positions = (max_rank + 1) - rank_labels
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(rank_labels)
    ax.set_ylabel("Best ranking")
    ax.set_ylim(max_rank - top_k, max_rank)

    # x-axis labels
    ax.set_xticks(df_plot["x"])
    ax.set_xticklabels(labels, rotation=90, fontsize=7)

    # styling
    ax.grid(axis="y", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.show()


def interactive_blastp_results(
    hits_df: pd.DataFrame,
    scores_df: pd.DataFrame = None,
    query_sequence: str = None,
    ligand_df: pd.DataFrame = None,
):
    """Create an interactive widget to explore BLASTP results and view rankings.
    
    Args:
        hits_df: DataFrame with BLASTP results (from blastp_query)
        scores_df: Optional DataFrame with scores (for mapping sseqid to ligand names and sequences)
        query_sequence: Query sequence string (required for alignment plotting)
        ligand_df: Optional DataFrame with ligands including 'ROMol' column for molecule visualization
        
    Returns:
        Interactive widget for exploring hits
    """
    if not IPYWIDGETS_AVAILABLE:
        print("ipywidgets not available. Please install with: pip install ipywidgets")
        return hits_df
    
    from ralph.alignment import plot_sequence_alignment
    
    # Create output area for plots and displays
    output = Output()
    
    # Hit selection dropdown
    hit_dropdown = widgets.Dropdown(
        options=[('Select a hit...', None)] + [(f"{row['sseqid']} (bitscore={row['bitscore']:.1f}, pident={row['pident']:.1f}%)", idx) 
                                                for idx, row in hits_df.iterrows()],
        description='Select Hit:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='600px')
    )
    
    # Top K slider
    # Set max value based on top_rankings length, with a minimum of 100
    max_k = 58
    top_k_slider = widgets.IntSlider(
        value=min(58, max_k) if max_k > 0 else 58,
        min=1,
        max=max(max_k, 100),
        step=1,
        description='Top K:',
        style={'description_width': 'initial'}
    )
    
    def update_display(selected_hit_idx, top_k):
        with output:
            output.clear_output(wait=True)
            
            if selected_hit_idx is None:
                display(HTML("<p>Please select a hit to view alignment, top molecules, and rankings.</p>"))
                return
            
            try:
                selected_hit = hits_df.iloc[selected_hit_idx]
                sseqid = selected_hit['sseqid']
                
                display(HTML(f"<h2>Results for Selected Hit: {sseqid}</h2>"))
                
                # 1. Sequence Alignment Plot
                if query_sequence is None:
                    display(HTML("<p><em>Query sequence not provided. Cannot display alignment.</em></p>"))
                elif scores_df is None:
                    display(HTML("<p><em>Scores DataFrame not provided. Cannot extract target sequence for alignment.</em></p>"))
                else:
                    try:
                        matching_scores = scores_df[scores_df['luxr_name'] == sseqid]
                        if len(matching_scores) > 0 and 'sequence' in matching_scores.columns:
                            target_seq = matching_scores['sequence'].values[0]
                            display(HTML("<h3>1. Sequence Alignment</h3>"))
                            fig = plot_sequence_alignment(
                                query_seq=query_sequence,
                                target_seq=target_seq,
                                query_name="Query",
                                target_name=sseqid
                            )
                            plt.show()
                        else:
                            display(HTML(f"<p><em>Sequence not found in scores_df for {sseqid}.</em></p>"))
                    except Exception as e:
                        display(HTML(f"<p><em>Error displaying alignment: {str(e)}</em></p>"))
                
                # 2. Top K Molecules Display
                # Get rankings from scores_df for the selected hit
                if scores_df is None:
                    display(HTML(f"<h3>2. Top {top_k} Molecules</h3>"))
                    display(HTML("<p><em>Scores DataFrame not provided. Cannot display molecules.</em></p>"))
                else:
                    try:
                        matching_scores = scores_df[scores_df['luxr_name'] == sseqid]
                        if len(matching_scores) == 0:
                            display(HTML(f"<h3>2. Top {top_k} Molecules</h3>"))
                            display(HTML(f"<p><em>No scores found for {sseqid}.</em></p>"))
                        else:
                            display(HTML(f"<h3>2. Top {top_k} Molecules</h3>"))
                            # Get ensemble columns for the selected hit
                            hit_row = matching_scores.iloc[0]
                            
                            # Extract molecule names from ensemble columns (ensemble_1, ensemble_2, ..., ensemble_top_k)
                            top_k_molecule_names = []
                            for i in range(1, top_k + 1):
                                col_name = f'ensemble_{i}'
                                if col_name in hit_row.index:
                                    mol_name = hit_row[col_name]
                                    if pd.notna(mol_name) and str(mol_name).strip():
                                        top_k_molecule_names.append((i, str(mol_name).strip()))
                            
                            if len(top_k_molecule_names) == 0:
                                display(HTML(f"<p><em>No ensemble columns found for {sseqid}.</em></p>"))
                            else:
                                # Extract just the names for display
                                molecule_names = [name for _, name in top_k_molecule_names]
                                
                                if ligand_df is not None:
                                    try:
                                        # Try to find a column that matches molecule names (could be 'AHLid', 'name', etc.)
                                        name_col = None
                                        for col in ['AHLid', 'name', 'id', 'ligand_id']:
                                            if col in ligand_df.columns:
                                                name_col = col
                                                break
                                        
                                        if name_col is not None:
                                            # Filter ligand_df to get molecules matching top K names
                                            topk_molecules = ligand_df[ligand_df[name_col].isin(molecule_names)].copy()
                                            
                                            if 'ROMol' in topk_molecules.columns and len(topk_molecules) > 0:
                                                try:
                                                    from rdkit.Chem.Draw import MolsToGridImage
                                                    
                                                    # Create rank mapping and sort to match ranking order
                                                    rank_to_name = {rank: name for rank, name in top_k_molecule_names}
                                                    name_to_rank = {name: rank for rank, name in top_k_molecule_names}
                                                    topk_molecules['_rank'] = topk_molecules[name_col].map(name_to_rank)
                                                    topk_molecules = topk_molecules.sort_values('_rank').dropna(subset=['ROMol'])
                                                    
                                                    mols = list(topk_molecules['ROMol'])
                                                    legends = []
                                                    for _, row in topk_molecules.iterrows():
                                                        rank = int(row['_rank'])
                                                        name = row[name_col]
                                                        legends.append(f"Rank: {rank}\t[{name}]")
                                                    
                                                    img = MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200), legends=legends)
                                                    display(img)
                                                except ImportError:
                                                    display(HTML("<p><em>RDKit not available. Displaying molecule names instead.</em></p>"))
                                                    display(HTML(f"<p>Top {top_k} molecules: {', '.join(molecule_names)}</p>"))
                                                except Exception as e:
                                                    display(HTML(f"<p><em>Error displaying molecules: {str(e)}</em></p>"))
                                                    display(HTML(f"<p>Top {top_k} molecules: {', '.join(molecule_names)}</p>"))
                                            else:
                                                display(HTML(f"<p>Top {top_k} molecules: {', '.join(molecule_names)}</p>"))
                                        else:
                                            display(HTML(f"<p>Top {top_k} molecules: {', '.join(molecule_names)}</p>"))
                                    except Exception as e:
                                        display(HTML(f"<p><em>Error processing ligand_df: {str(e)}</em></p>"))
                                        display(HTML(f"<p>Top {top_k} molecules: {', '.join(molecule_names)}</p>"))
                                else:
                                    # Display molecule names as text
                                    display(HTML(f"<p>Top {top_k} molecules: {', '.join(molecule_names)}</p>"))
                    except Exception as e:
                        display(HTML(f"<p><em>Error extracting molecules: {str(e)}</em></p>"))
                
                # 3. Clustered Rank Plot
                # Get rankings from scores_df for the selected hit
                if scores_df is None:
                    display(HTML("<h3>3. Rankings Plot</h3>"))
                    display(HTML("<p><em>Scores DataFrame not provided. Cannot display rankings plot.</em></p>"))
                else:
                    try:
                        matching_scores = scores_df[scores_df['luxr_name'] == sseqid]
                        if len(matching_scores) == 0:
                            display(HTML("<h3>3. Rankings Plot</h3>"))
                            display(HTML(f"<p><em>No scores found for {sseqid}.</em></p>"))
                        else:
                            display(HTML("<h3>3. Rankings Plot</h3>"))
                            # Get ensemble columns for the selected hit
                            hit_row = matching_scores.iloc[0]
                            
                            # Extract molecule names and ranks from ensemble columns
                            plot_data = []
                            for i in range(1, top_k + 1):
                                col_name = f'ensemble_{i}'
                                if col_name in hit_row.index:
                                    mol_name = hit_row[col_name]
                                    if pd.notna(mol_name) and str(mol_name).strip():
                                        plot_data.append({
                                            'label': str(mol_name).strip(),
                                            'cluster': 0,
                                            'best_rank': i
                                        })
                            
                            if len(plot_data) == 0:
                                display(HTML(f"<p><em>No ensemble columns found for {sseqid}.</em></p>"))
                            else:
                                # Create DataFrame for clustered_rank_plot
                                df_plot = pd.DataFrame(plot_data)
                                
                                # Highlight top-ranked molecule
                                top_label = df_plot.sort_values('best_rank').iloc[0]['label']
                                
                                clustered_rank_plot(df_plot, top_k=top_k, gap=2.0, highlight_labels=(top_label,))
                    except Exception as e:
                        display(HTML(f"<p><em>Error displaying rankings plot: {str(e)}</em></p>"))
            
            except Exception as e:
                display(HTML(f"<p><em>Error processing selected hit: {str(e)}</em></p>"))
    
    # Create interactive widget
    # Note: interactive() automatically creates widgets for the parameters, so we don't need to add them separately
    interact_widget = interactive(update_display, selected_hit_idx=hit_dropdown, top_k=top_k_slider)
    
    # Create layout with interactive widget (which contains the dropdown and slider) and output
    main_layout = VBox([
        interact_widget,
        output
    ])
    
    return main_layout


def interactive_alignment_viewer(
    hits_df: pd.DataFrame,
    query_sequence: str,
    db_path: str,
    scores_df: pd.DataFrame = None,
    top_rankings: pd.DataFrame = None,
):
    """Create an interactive widget to view sequence alignments and rankings.
    
    Args:
        hits_df: DataFrame with BLASTP results
        query_sequence: Query protein sequence
        db_path: Path to BLAST database (for extracting target sequences)
        scores_df: Optional DataFrame with scores (for mapping sseqid to ligand names)
        top_rankings: Optional DataFrame with aggregated rankings
        
    Returns:
        Interactive widget for viewing alignments
    """
    if not IPYWIDGETS_AVAILABLE:
        print("ipywidgets not available. Please install with: pip install ipywidgets")
        return None
    
    from ralph.alignment import extract_sequence_from_blastdb, plot_sequence_alignment
    
    # Create hit selection dropdown
    hit_options = [('Select a hit...', None)] + [
        (f"{row['sseqid']} (bitscore={row['bitscore']:.1f}, pident={row['pident']:.1f}%)", idx)
        for idx, row in hits_df.iterrows()
    ]
    
    hit_dropdown = widgets.Dropdown(
        options=hit_options,
        description='Select Hit:',
        style={'description_width': 'initial'}
    )
    
    # Output area for alignment plot and ranking info
    alignment_output = Output()
    ranking_output = Output()
    
    def show_alignment(selected_hit_idx):
        alignment_output.clear_output(wait=True)
        ranking_output.clear_output(wait=True)
        
        if selected_hit_idx is None:
            with alignment_output:
                print("Please select a hit to view alignment.")
            return
        
        try:
            selected_hit = hits_df.iloc[selected_hit_idx]
            sseqid = selected_hit['sseqid']
            
            # Extract target sequence from database
            target_seq = scores_df[scores_df['luxr_name'] == sseqid]['sequence'].values[0]
            
            # Create alignment plot
            with alignment_output:
                fig = plot_sequence_alignment(
                    query_seq=query_sequence,
                    target_seq=target_seq,
                    query_name="Query",
                    target_name=sseqid
                )
                plt.show()
            
            # Show ranking information
            if top_rankings is not None and scores_df is not None:
                with ranking_output:
                    matching_scores = scores_df[scores_df['luxr_name'] == sseqid]
                    if len(matching_scores) > 0:
                        # Get unique ligand names from ensemble columns
                        ensemble_cols = [c for c in matching_scores.columns if c.startswith('ensemble_') and not c.endswith('std')]
                        ligand_names = set()
                        for col in ensemble_cols:
                            ligand_names.update(matching_scores[col].dropna().unique())
                        ligand_names = [str(name).strip() for name in ligand_names if name and str(name).strip()]
                        
                        print(f"\n{'='*60}")
                        print(f"Ranking Information for {sseqid}")
                        print(f"{'='*60}")
                        print(f"\nBLASTP Statistics:")
                        print(f"  Bitscore: {selected_hit['bitscore']:.2f}")
                        print(f"  E-value: {selected_hit['evalue']:.2e}")
                        print(f"  % Identity: {selected_hit['pident']:.2f}%")
                        print(f"  Alignment Length: {selected_hit['length']:.0f}")
                        
                        if len(ligand_names) > 0:
                            print(f"\nAssociated Ligands: {', '.join(ligand_names)}")
                            
                            ligand_rankings = top_rankings[top_rankings['name'].isin(ligand_names)]
                            if len(ligand_rankings) > 0:
                                print(f"\nRankings:")
                                print(ligand_rankings[['name', 'ranking', 'score', 'avg_position']].to_string(index=False))
                            else:
                                print("\nNo rankings found for these ligands in top_rankings.")
                        else:
                            print("\nNo ligands found associated with this sequence ID.")
                    else:
                        print(f"No scores found for sequence ID: {sseqid}")
            else:
                with ranking_output:
                    print(f"\n{'='*60}")
                    print(f"BLASTP Statistics for {sseqid}")
                    print(f"{'='*60}")
                    print(f"  Bitscore: {selected_hit['bitscore']:.2f}")
                    print(f"  E-value: {selected_hit['evalue']:.2e}")
                    print(f"  % Identity: {selected_hit['pident']:.2f}%")
                    print(f"  Alignment Length: {selected_hit['length']:.0f}")
                    if top_rankings is None:
                        print("\nRanking information not available (top_rankings not provided)")
        
        except Exception as e:
            with alignment_output:
                print(f"Error displaying alignment: {str(e)}")
            with ranking_output:
                print(f"Error: {str(e)}")
    
    # Create interactive widget
    interact_widget = interactive(show_alignment, selected_hit_idx=hit_dropdown)
    
    # Layout
    layout = VBox([
        interact_widget,
        alignment_output,
        ranking_output
    ])
    
    return layout
