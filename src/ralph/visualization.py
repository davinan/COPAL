"""Visualization functions for ranking results."""

import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import HTML, display

try:
    import ipywidgets as widgets
    from ipywidgets import Output, VBox, interactive
    IPYWIDGETS_AVAILABLE = True
except ImportError:
    IPYWIDGETS_AVAILABLE = False
    widgets = None


# Ranking blocks shown per hit: COPAL ensemble first, then each individual model.
# Column prefixes match the shipped ralph_scores.csv ({prefix}_1..{prefix}_58);
# score columns match ralph_model_scores.csv.
_MODEL_DISPLAY = [
    ("ensemble", "COPAL (ensemble)"),
    ("af3", "AF3"),
    ("rf3", "RF3"),
    ("boltz2", "Boltz2"),
    ("boltz2affinity", "Boltz2 Affinity"),
]
_SCORE_COL = {"ensemble": "copal", "af3": "af3", "rf3": "rf3",
              "boltz2": "boltz2", "boltz2affinity": "boltz2affinity"}
_SCORE_LABEL = {"ensemble": "sum-rank", "af3": "ipTM", "rf3": "ipTM",
                "boltz2": "ipTM", "boltz2affinity": "affinity"}


def _ranked_names(hit_row, prefix, top_k):
    """Return [(rank, ahl_id), ...] from the {prefix}_1..{prefix}_top_k columns."""
    out = []
    for i in range(1, top_k + 1):
        col = f"{prefix}_{i}"
        if col in hit_row.index:
            v = hit_row[col]
            if pd.notna(v) and str(v).strip():
                out.append((i, str(v).strip()))
    return out


def _render_ligand_grid(ranked, ligand_df, per_hit, score_col, score_label):
    """Display a rank-ordered RDKit grid; annotate each ligand with its model score.

    ranked:    [(rank, ahl_id), ...]
    per_hit:   DataFrame of this receptor's per-pair scores, indexed by ahl (or None).
    """
    names = [n for _, n in ranked]
    name_col = next(
        (c for c in ["AHLid", "name", "id", "ligand_id"] if c in ligand_df.columns),
        None,
    )
    if name_col is None:
        display(HTML(f"<p>{', '.join(names)}</p>"))
        return
    sub = ligand_df[ligand_df[name_col].isin(names)].copy()
    if "ROMol" not in sub.columns or len(sub) == 0:
        display(HTML(f"<p>{', '.join(names)}</p>"))
        return
    try:
        from rdkit.Chem.Draw import MolsToGridImage
    except ImportError:
        display(HTML(f"<p>{', '.join(names)}</p>"))
        return
    rankmap = {n: r for r, n in ranked}
    sub["_rank"] = sub[name_col].map(rankmap)
    sub = sub.sort_values("_rank").dropna(subset=["ROMol"])
    legends = []
    for _, row in sub.iterrows():
        nm = row[name_col]
        rk = int(row["_rank"])
        leg = f"#{rk} {nm}"
        if per_hit is not None and score_col in per_hit.columns and nm in per_hit.index:
            sc = per_hit.loc[nm, score_col]
            if pd.notna(sc):
                leg += f"  {score_label}={float(sc):.3g}"
        legends.append(leg)
    img = MolsToGridImage(
        list(sub["ROMol"]), molsPerRow=4, subImgSize=(190, 190), legends=legends
    )
    display(img)


def interactive_blastp_results(
    hits_df: pd.DataFrame,
    scores_df: pd.DataFrame = None,
    query_sequence: str = None,
    ligand_df: pd.DataFrame = None,
    model_scores_df: pd.DataFrame = None,
):
    """Interactive widget to explore BLASTP hits with per-model ligand rankings.

    For the selected hit, shows: the query-vs-hit sequence alignment; and the
    top-K ligands for the model chosen in the dropdown (COPAL ensemble by
    default, or AF3, RF3, Boltz2, Boltz2 Affinity), each ligand annotated with
    that model's score. The header line carries the hit's mpKT and, for RefAHL
    benchmark receptors, the known cognate ligand.

    Args:
        hits_df: BLASTP results (from blastp_query).
        scores_df: wide rankings table (from load_scores).
        query_sequence: query sequence (for the alignment plot).
        ligand_df: ligands with a 'ROMol' column (from load_ligands).
        model_scores_df: per-pair scores indexed by (luxr_name, ahl)
            (from load_model_scores); optional, used to annotate scores.
    """
    if not IPYWIDGETS_AVAILABLE:
        print("ipywidgets not available. Please install with: pip install ipywidgets")
        return hits_df

    from ralph.alignment import plot_sequence_alignment

    output = Output()
    hit_dropdown = widgets.Dropdown(
        options=[('Select a hit...', None)] + [
            (f"{row['sseqid']} (bitscore={row['bitscore']:.1f}, pident={row['pident']:.1f}%)", idx)
            for idx, row in hits_df.iterrows()
        ],
        description='Select Hit:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='600px'),
    )
    top_k_slider = widgets.IntSlider(
        value=10, min=1, max=58, step=1,
        description='Top K:', style={'description_width': 'initial'},
    )
    model_dropdown = widgets.Dropdown(
        options=[(label, prefix) for prefix, label in _MODEL_DISPLAY],
        value='ensemble',
        description='Model:',
        style={'description_width': 'initial'},
    )

    def update_display(selected_hit_idx, model, top_k):
        with output:
            output.clear_output(wait=True)
            if selected_hit_idx is None:
                display(HTML("<p>Select a hit to view its alignment and ligand rankings.</p>"))
                return
            try:
                selected_hit = hits_df.iloc[selected_hit_idx]
                sseqid = selected_hit['sseqid']
                display(HTML(f"<h2>Results for {sseqid}</h2>"))

                if scores_df is None:
                    display(HTML("<p><em>scores_df not provided.</em></p>"))
                    return
                matching = scores_df[scores_df['luxr_name'] == sseqid]
                if len(matching) == 0:
                    display(HTML(f"<p><em>No scores found for {sseqid}.</em></p>"))
                    return
                hit_row = matching.iloc[0]

                # Header line: source, mpKT, known cognate (RefAHL), organism.
                bits = []
                if 'source' in hit_row.index and pd.notna(hit_row['source']):
                    bits.append(f"source: <b>{hit_row['source']}</b>")
                if 'mpKT' in hit_row.index and pd.notna(hit_row['mpKT']):
                    bits.append(f"mpKT: {float(hit_row['mpKT']):.3f}")
                if hit_row.get('source') == 'RefAHL' and pd.notna(hit_row.get('true_ahl_id')):
                    bits.append(f"known cognate: <b>{hit_row['true_ahl_id']}</b>")
                if 'organism' in hit_row.index and pd.notna(hit_row['organism']):
                    bits.append(f"<span style='color:#666'>{str(hit_row['organism'])[:90]}</span>")
                if bits:
                    display(HTML(" &nbsp;|&nbsp; ".join(bits)))

                # 1. Sequence alignment.
                if query_sequence is not None and 'sequence' in hit_row.index:
                    try:
                        display(HTML("<h3>1. Sequence alignment</h3>"))
                        plot_sequence_alignment(
                            query_seq=query_sequence,
                            target_seq=hit_row['sequence'],
                            query_name="Query",
                            target_name=sseqid,
                        )
                        plt.show()
                    except Exception as e:
                        display(HTML(f"<p><em>Alignment error: {e}</em></p>"))

                # Per-hit scores (indexed by ahl) for annotation.
                per_hit = None
                if model_scores_df is not None and len(model_scores_df):
                    try:
                        per_hit = model_scores_df.xs(sseqid, level='luxr_name')
                    except KeyError:
                        per_hit = None

                # 2. Ranked ligands for the model selected in the dropdown.
                label = dict(_MODEL_DISPLAY)[model]
                ranked = _ranked_names(hit_row, model, top_k)
                display(HTML(f"<h3>2. Top {top_k} ligands — {label}</h3>"))
                if not ranked:
                    display(HTML("<p><em>No ranking available for this model.</em></p>"))
                elif ligand_df is not None:
                    _render_ligand_grid(
                        ranked, ligand_df, per_hit,
                        _SCORE_COL[model], _SCORE_LABEL[model],
                    )
                else:
                    display(HTML("<p>" + ", ".join(f"#{r} {n}" for r, n in ranked) + "</p>"))
            except Exception as e:
                display(HTML(f"<p><em>Error processing selected hit: {e}</em></p>"))

    interact_widget = interactive(
        update_display,
        selected_hit_idx=hit_dropdown,
        model=model_dropdown,
        top_k=top_k_slider,
    )
    return VBox([interact_widget, output])


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

    from ralph.alignment import plot_sequence_alignment

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
                        print("\nBLASTP Statistics:")
                        print(f"  Bitscore: {selected_hit['bitscore']:.2f}")
                        print(f"  E-value: {selected_hit['evalue']:.2e}")
                        print(f"  % Identity: {selected_hit['pident']:.2f}%")
                        print(f"  Alignment Length: {selected_hit['length']:.0f}")

                        if len(ligand_names) > 0:
                            print(f"\nAssociated Ligands: {', '.join(ligand_names)}")

                            ligand_rankings = top_rankings[top_rankings['name'].isin(ligand_names)]
                            if len(ligand_rankings) > 0:
                                print("\nRankings:")
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
