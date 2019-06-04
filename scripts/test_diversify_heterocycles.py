"""
Written by Eric Dang.
github: https://github.com/e-dang
email: erickdang@g.ucla.edu
"""

import pytest
from mock import Mock, patch
from rdkit import Chem

from database import MolDatabase
from diversify_heterocycles import DiversifyHeterocycles
from pathlib import Path


@pytest.fixture
def mock_mongodatabase():
    return Mock(spec=MolDatabase)


fp_in = ['chemdraw/pre_monomer/heterocycles_likely1.sdf']
fp_out = 'smiles/pre_monomer/heterocycles.json'
project_dir = Path(__file__).resolve().parents[1]
@pytest.fixture
def div():
    return DiversifyHeterocycles(fp_in, fp_out)


def test_init(monkeypatch, mock_mongodatabase):
    DB = Mock(return_value=mock_mongodatabase)
    monkeypatch.setattr('diversify_heterocycles.MolDatabase', DB)
    div = DiversifyHeterocycles(fp_in, fp_out)
    # assert div.fp_in == [str(project_dir / fp_in[0])]
    # assert div.fp_out == str(project_dir / fp_out)
    assert div.mol_db == mock_mongodatabase
    assert div.groups == ['likely1']


TEST_CASES1 = [
    ('c1ccoc1', f'[CH3:{DiversifyHeterocycles.CONN_MN}]', set(['Cc1ccco1', 'Cc1ccoc1'])),
    ('c1ccoc1', f'[CH3][CH2:{DiversifyHeterocycles.CONN_MN}]', set(['CCc1ccoc1', 'CCc1ccco1'])),
    ('c1ccoc1', f'[CH3][CH2][CH2:{DiversifyHeterocycles.CONN_MN}]', set(['CCCc1ccco1', 'CCCc1ccoc1'])),
    ('c1cc[nH]c1', f'[CH3:{DiversifyHeterocycles.CONN_MN}]', set(['Cc1cc[nH]c1', 'Cn1cccc1', 'Cc1ccc[nH]1'])),
    ('c1cc[nH]c1', f'[CH3][CH2:{DiversifyHeterocycles.CONN_MN}]', set(['CCc1ccc[nH]1', 'CCn1cccc1', 'CCc1cc[nH]c1'])),
    ('c1cc[nH]c1', f'[CH3][CH2][CH2:{DiversifyHeterocycles.CONN_MN}]',
     set(['CCCn1cccc1', 'CCCc1cc[nH]c1', 'CCCc1ccc[nH]1']))
]
@pytest.mark.parametrize('mol, connection, expected', TEST_CASES1)
def test_alternate_connection_point(div, mol, connection, expected):
    mol = Chem.MolFromSmiles(mol)
    assert div.alternate_connection_point(mol, connection) == expected


@pytest.mark.parametrize('mol', ['c1ccoc1'])
@pytest.mark.parametrize('connection', ['[CH3]', '[CH3][CH2]', '[CH3][CH2][CH2]'])
def test_alternate_connection_point_exceptions(div, mol, connection):
    with pytest.raises(SystemExit):
        div.alternate_connection_point(mol, connection)


TEST_CASES2 = [
    (['Cc1ccco1'], 'c1ccoc1', [0, 1], 'likely1', [{'heterocycle': 'Cc1ccco1',
                                                   'parent': 'c1ccoc1',
                                                   'modification': [0, 1],
                                                   'group': 'likely1'}])
]
@pytest.mark.parametrize('mols, parent, modifications, group, expected', TEST_CASES2)
def test_accumulate_mols(div, mols, parent, modifications, group, expected):
    assert div.data == []
    div.accumulate_mols(mols, parent, modifications, group)
    assert div.data == expected


def test_diversify(div):
    with patch('diversify_heterocycles.read_mols') as mock_read_mols:
        mock_read_mols.return_value = [Chem.MolFromSmiles('c1ccoc1')]
        div.connections = [(f'[CH3:{DiversifyHeterocycles.CONN_MN}]', [0, 1])]
        assert div.diversify() == True
        assert sorted(div.data, key=lambda x: x['heterocycle']) == sorted([{'heterocycle': 'Cc1ccoc1',
                                                                            'parent': 'c1ccoc1',
                                                                            'modification': [0, 1],
                                                                            'group': 'likely1'},
                                                                           {'heterocycle': 'Cc1ccco1',
                                                                            'parent': 'c1ccoc1',
                                                                            'modification': [0, 1],
                                                                            'group': 'likely1'}],
                                                                          key=lambda x: x['heterocycle'])


TEST_CASES3 = [
    (True, True),
    (False, False)
]
@pytest.mark.parametrize('json_flag, db_flag', TEST_CASES3, ids=['True', 'False'])
def test_save_data(div, json_flag, db_flag):
    with patch('diversify_heterocycles.Base.write_json') as mock_json, patch('diversify_heterocycles.MolDatabase.insert_heterocycles') as mock_db:
        div.json_flag = json_flag
        div.db_flag = db_flag
        assert div.save_data()
        if json_flag:
            mock_json.assert_called_once_with(div.data, div.fp_out)
        else:
            mock_json.assert_not_called()
        if db_flag:
            mock_db.assert_called_once_with(div.data)
        else:
            mock_db.assert_not_called()
