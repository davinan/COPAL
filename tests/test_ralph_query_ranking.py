"""Query ranking must use the benchmark's equal-weight reciprocal ranks."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from ralph.ranking import aggregate_blast_hits
from ralph.visualization import interactive_blastp_results


@pytest.fixture
def score_rows():
    return pd.DataFrame({
        'luxr_name': ['HIT_A', 'HIT_B', 'LOW_ID', 'IDENTICAL'],
        'ensemble_1': ['AHL01', 'AHL02', 'AHL03', 'AHL03'],
        'ensemble_2': ['AHL02', 'AHL03', 'AHL02', 'AHL01'],
        'ensemble_3': ['AHL03', 'AHL01', 'AHL01', 'AHL02'],
        'mpKT': [0.01, 0.99, 0.5, 1.0],
    })


@pytest.fixture
def hits():
    return pd.DataFrame({
        'sseqid': ['HIT_A', 'HIT_B', 'LOW_ID', 'IDENTICAL', 'UNAVAILABLE'],
        'pident': [30.0, 99.9, 29.9, 100.0, 50.0],
        'bitscore': [100.0, 300.0, 500.0, 600.0, 150.0],
    })


def test_benchmark_filter_and_unweighted_reciprocal_ranks(hits, score_rows):
    result = aggregate_blast_hits(hits, score_rows)
    assert result.name.tolist() == ['AHL02', 'AHL01', 'AHL03']
    assert result.score.tolist() == pytest.approx([1.5, 1 + 1/3, 1/2 + 1/3])
    assert result['count'].tolist() == [2, 2, 2]
    assert result.avg_position.tolist() == [1.5, 2.0, 2.5]


def test_multiple_blast_alignments_do_not_reweight_a_receptor(hits, score_rows):
    duplicated_hits = pd.concat([hits, hits.iloc[:1]], ignore_index=True)
    pd.testing.assert_frame_equal(
        aggregate_blast_hits(duplicated_hits, score_rows),
        aggregate_blast_hits(hits, score_rows),
    )


@pytest.mark.parametrize('case', ['empty', 'filtered', 'unavailable'])
def test_no_usable_hits(case, hits, score_rows):
    if case == 'empty':
        hits = hits.iloc[:0]
    elif case == 'filtered':
        hits = hits[hits.sseqid.isin(['LOW_ID', 'IDENTICAL'])]
    else:
        hits = hits[hits.sseqid.eq('UNAVAILABLE')]
    assert aggregate_blast_hits(hits, score_rows) is None


def test_published_benchmark_cli_uses_same_calculation(hits, score_rows, monkeypatch):
    path = Path(__file__).resolve().parents[1]/'baselines/ralph_blast/run_ralph_blast.py'
    spec = importlib.util.spec_from_file_location('benchmark_cli', path)
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)
    monkeypatch.setattr(benchmark, 'blastp_query', lambda **_: hits)
    observed = benchmark.run_blast_and_aggregate('SEQUENCE', score_rows, 'unused', 50, 30.0, 100.0)
    pd.testing.assert_frame_equal(observed, aggregate_blast_hits(hits, score_rows))


@pytest.mark.parametrize('has_ranking', [False, True])
def test_query_dropdown_without_selected_hit(hits, score_rows, monkeypatch, has_ranking):
    import ralph.visualization as visual

    rendered = []
    grids = []
    monkeypatch.setattr(visual, 'display', rendered.append)
    monkeypatch.setattr(visual, '_render_ligand_grid', lambda *args: grids.append(args))
    ranking = aggregate_blast_hits(hits, score_rows) if has_ranking else None
    panel = interactive_blastp_results(hits, score_rows, ligand_df=pd.DataFrame(), query_ranking=ranking)
    control = panel.children[0]
    inputs = {child.description: child for child in control.children if hasattr(child, 'description')}
    assert inputs['Select Hit:'].value is None
    inputs['Model:'].value = 'ralph_blast'
    assert inputs['Select Hit:'].disabled
    html = '\n'.join(obj.data for obj in rendered if hasattr(obj, 'data'))
    assert 'RALPH-BLAST ligand ranking for the query' in html
    if has_ranking:
        assert grids[-1][0] == [(1, 'AHL02'), (2, 'AHL01'), (3, 'AHL03')]
        assert grids[-1][3] == 'score'
        inputs['Top K:'].value = 1
        assert grids[-1][0] == [(1, 'AHL02')]
    else:
        assert 'No query ranking is available' in html
        assert not grids
    inputs['Model:'].value = 'ensemble'
    assert not inputs['Select Hit:'].disabled
    inputs['Select Hit:'].value = 0
    assert grids[-1][0][:3] == [(1, 'AHL01'), (2, 'AHL02'), (3, 'AHL03')][:inputs['Top K:'].value]
