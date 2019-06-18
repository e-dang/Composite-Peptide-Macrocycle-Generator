"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

import pytest
from mock import Mock, patch
from rdkit import Chem

from macrocycles.scripts.database import MolDatabase
from macrocycles.scripts.side_chain_modifier import SideChainModifier


@pytest.fixture
def mock_mongodatabase():
    return Mock(spec=MolDatabase)


@pytest.fixture
def div():
    return SideChainModifier(['chemdraw/pre_monomer/heterocycles_likely1.sdf'], 'smiles/pre_monomer/heterocycles.json')


# FP_IN = [
#     (1),
#     (1.),
#     ([1, 2]),
#     ([1., 2.]),
#     ((1, 2)),
#     ((1., 2.))
# ]
# @pytest.mark.parametrize('fp_in', FP_IN)
# def test_fp_in_exceptions(div, fp_in):
#     with pytest.raises(TypeError):
#         div.fp_in = fp_in


# FP_OUT = [
#     (1),
#     (1.),
#     (['chemdraw/pre_monomer/side_chains_likely1.sdf']),
#     (('chemdraw/pre_monomer/side_chains_likely1.sdf'))
# ]
# @pytest.mark.parametrize('fp_out', FP_OUT)
# def test_fp_out_exceptions(div, fp_out):
#     with pytest.raises(TypeError):
#         div.fp_out = fp_out


# MOL_DB = [
#     (1),
#     ('test_db'),
#     (()),
#     ([])
# ]
# @pytest.mark.parametrize('mol_db', MOL_DB)
# def test_mol_db_exceptions(div, mol_db):
#     with pytest.raises(TypeError):
#         div.mol_db = mol_db


# GROUPS = [
#     (1),
#     (1.),
#     ([1]),
#     ([1.])
# ]
# @pytest.mark.parametrize('groups', GROUPS)
# def test_groups_exceptions(div, groups):
#     with pytest.raises(TypeError):
#         div.groups = groups


# DATA = [
#     (1),
#     (1.),
#     ('test'),
#     ({}),
#     ([]),
#     ({'parent': 'parent',
#       'modifications': 'modifications',
#       'group': 'group'}),
#     ({'side_chain': 'side_chain',
#       'modifications': 'modifications',
#       'group': 'group'}),
#     ({'side_chain': 'side_chain',
#       'parent': 'parent',
#       'group': 'group'}),
#     ({'side_chain': 'side_chain',
#       'parent': 'parent',
#       'modifications': 'modifications'})
# ]
# @pytest.mark.parametrize('data', DATA)
# def test_data_exceptions(div, data):
#     with pytest.raises(TypeError):
#         div.data = data


TEST_CASES1 = [
    ('c1cc[nH]c1', f'[CH3:{SideChainModifier.CONN_MN}]', set(['Cc1cc[nH]c1', 'Cn1cccc1', 'Cc1ccc[nH]1'])),
    ('c1cc[nH]c1', f'[CH3][CH2:{SideChainModifier.CONN_MN}]', set(['CCc1ccc[nH]1', 'CCn1cccc1', 'CCc1cc[nH]c1'])),
    ('c1cc[nH]c1', f'[CH3][CH2][CH2:{SideChainModifier.CONN_MN}]',
     set(['CCCn1cccc1', 'CCCc1cc[nH]c1', 'CCCc1ccc[nH]1']))
]
@pytest.mark.parametrize('mol, connection, expected', TEST_CASES1)
def test_alternate_connection_point(div, mol, connection, expected):
    """
    Test correct output for modifying side chains
    """
    mol = Chem.MolFromSmiles(mol)
    assert div.alternate_connection_point(mol, connection) == expected


@pytest.mark.parametrize('mol', ['c1ccoc1'])
@pytest.mark.parametrize('connection', ['[CH3]', '[CH3][CH2]', '[CH3][CH2][CH2]'])
def test_alternate_connection_point_exceptions(div, mol, connection):
    """
    Test if SystemExit exception is thrown when connection chains are not atom mapped.
    """
    with pytest.raises(SystemExit):
        div.alternate_connection_point(mol, connection)


TEST_CASES2 = [
    (['Cc1ccco1'], 'c1ccoc1', [0, 1], 'likely1', [{'side_chain': 'Cc1ccco1',
                                                   'parent': 'c1ccoc1',
                                                   'modifications': [0, 1],
                                                   'group': 'likely1'}])
]
@pytest.mark.parametrize('side_chain, parent, modifications, group, expected', TEST_CASES2)
def test_accumulate_mols(div, side_chain, parent, modifications, group, expected):
    """
    Test that data is stored in class instance properly.
    """
    assert div.data == []   # empty before
    div.accumulate_mols(side_chain, parent, modifications, group)
    assert div.data == expected  # filled after


def test_diversify(div):
    """
    Test that function returns true and stores data in class instance.
    """
    with patch('side_chain_modifier.read_mols') as mock_read_mols:
        mock_read_mols.return_value = [Chem.MolFromSmiles('c1ccoc1')]
        div.connections = [(f'[CH3:{SideChainModifier.CONN_MN}]', [0, 1])]
        assert div.diversify() == True
        print(div.data)
        assert sorted(div.data, key=lambda x: x['side_chain']) == sorted([{'side_chain': 'Cc1ccoc1',
                                                                           'parent': 'c1ccoc1',
                                                                           'modifications': [0, 1],
                                                                           'group': 'likely1'},
                                                                          {'side_chain': 'Cc1ccco1',
                                                                           'parent': 'c1ccoc1',
                                                                           'modifications': [0, 1],
                                                                           'group': 'likely1'}],
                                                                         key=lambda x: x['side_chain'])


TEST_CASES3 = [
    (True, True),
    (False, False)
]
@pytest.mark.parametrize('json_flag, db_flag', TEST_CASES3, ids=['True', 'False'])
def test_save_data(div, json_flag, db_flag):
    """
    Test data get saved properly when flags are set.
    """
    with patch('side_chain_modifier.Base.write_json') as mock_json, patch('side_chain_modifier.MolDatabase.insert_side_chains') as mock_db:
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
