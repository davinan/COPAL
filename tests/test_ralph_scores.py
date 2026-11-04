"""Regression coverage for the COPAL-only RALPH score file."""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

from ralph.database import load_model_scores
from ralph.visualization import _smiles_table


@pytest.mark.parametrize('suffix', ['.csv', '.csv.xz'])
@pytest.mark.parametrize('legacy', [False, True])
def test_load_only_copal_scores(tmp_path, suffix, legacy):
    expected = pd.DataFrame({'luxr_name':['RALPH00000'], 'ahl':['AHL01'], 'copal':[12.0]})
    source = expected.assign(af3=0.9, rf3=0.8, boltz2=0.7, boltz2affinity=-1.0) if legacy else expected
    source.to_csv(tmp_path/f'ralph_model_scores{suffix}',index=False)
    observed = load_model_scores(str(tmp_path))
    pd.testing.assert_frame_equal(observed, expected.set_index(['luxr_name','ahl']))


def test_table_omits_removed_score_column():
    row = [(1,'C8','CCC','')]
    html = _smiles_table(row, None)
    assert html.count('<th ') == 3
    assert html.count('<td ') == 3
    assert 'None' not in html
    scored = _smiles_table([(1,'C8','CCC','12')], 'sum-rank')
    assert scored.count('<th ') == 4
    assert 'sum-rank' in scored
    assert '>12</td>' in scored


def test_builder_exports_only_copal_scores(tmp_path, monkeypatch):
    path = Path(__file__).resolve().parents[1]/'baselines/ralph_blast/build_ralph_db.py'
    spec = importlib.util.spec_from_file_location('build_ralph_db',path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    ahls = [f'AHL{i:02d}' for i in range(1,59)]
    wide = pd.DataFrame({'luxr_name':['RALPH00000'],'source':['RALPH'],'organism':['test'],
                         'true_ahl_id':[None],'mpKT':[0.5],'sequence':['AAAA']})
    ranks = pd.DataFrame({f'{model}_{i}':[ahl] for model in builder.MODELS for i,ahl in enumerate(ahls,1)})
    wide = pd.concat([wide,ranks],axis=1)
    scores = pd.DataFrame({'luxr_name':['RALPH00000']*58,'ahl':ahls,'copal':range(4,62),
                           'af3':0.9,'rf3':0.8,'boltz2':0.7,'boltz2affinity':-1.0})
    monkeypatch.setattr(builder,'build_ralph',lambda *_:(wide,scores))
    ligands = tmp_path/'ligands.csv'
    pd.DataFrame({'AHLid':ahls,'abbreviation':ahls,'smi':['CCC']*58}).to_csv(ligands,index=False)
    monkeypatch.setattr(sys,'argv',['build_ralph_db','--xlsx','unused.xlsx','--ralph-v3','unused.csv',
                                   '--bench-rankings','unused.csv','--bench-ensemble','unused.csv',
                                   '--ligands',str(ligands),'--out-dir',str(tmp_path/'out'),'--no-refahl'])
    builder.main()
    observed = pd.read_csv(tmp_path/'out/ralph_model_scores.csv.xz')
    pd.testing.assert_frame_equal(observed,scores[['luxr_name','ahl','copal']])
