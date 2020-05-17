import numpy as np
import pytest
from rdkit import Chem

import cpmg.utils as utils


@pytest.fixture
def string_list_data():
    return ['asdijviu29478', 'dq908asdfghajkw', 'ac90ded']


@pytest.fixture
def string_int_list_data():
    return ['0', '10', '12']


@pytest.fixture
def empty_list_data():
    return []


@pytest.mark.parametrize('data, pred, expected_len1, expected_len2', [
    ([[1, 2, 3], [1, 2, 3, 4], [4, 3, 2, 4], [1, 2, 2], [1, 2, 3, 4]], lambda x: len(x) == 3, 2, 3),
    (map(np.array, [[1, 2, 3], [1, 2, 3, 4], [4, 3, 2, 4], [1, 2, 2], [1, 2, 3, 4]]), lambda x: len(x) == 3, 2, 3),
    (map(tuple, [[1, 2, 3], [1, 2, 3, 4], [4, 3, 2, 4], [1, 2, 2], [1, 2, 3, 4]]), lambda x: len(x) == 3, 2, 3),
])
def test_split(data, pred, expected_len1, expected_len2):

    eq_len, neq_len = utils.split(data, pred=pred)

    assert len(eq_len) == expected_len1
    assert len(neq_len) == expected_len2


def test_get_maximum_len(string_list_data):
    assert(utils.get_maximum(string_list_data, len) == 15)


def test_get_maximum_int(string_int_list_data):
    assert(utils.get_maximum(string_int_list_data, int) == 12)


@pytest.mark.parametrize('data,length', [
    ({'A': 1}, 1),
    ([{'A': 1}, {'B': 2}], 2),
    (({'A': 1}, {'B': 2}), 2),
    (map(lambda x: x, [{'A': 1}, {'B': 2}]), 2),

])
def test_to_list(data, length):
    ls_data = utils.to_list(data)

    assert isinstance(ls_data, (list, tuple))
    assert all(map(lambda x: isinstance(x, dict), ls_data))
    assert len(ls_data) == length


@pytest.mark.parametrize('empty_list_data,func', [('', int), ('', len)], indirect=['empty_list_data'])
def test_get_maximum_empty(empty_list_data, func):
    assert(utils.get_maximum(empty_list_data, func) == None)


@pytest.mark.parametrize('mol,expected_result', [(Chem.MolFromSmiles('C[CH2:1]C'), True), (Chem.MolFromSmiles('CCC'), False)])
def test_has_atom_map_nums(mol, expected_result):
    assert(utils.has_atom_map_nums(mol) == expected_result)


@pytest.mark.parametrize('mol,expected_result', [(Chem.MolFromSmiles('C[CH2:1][CH3:2]'), [(1, 1), (2, 2)]), (Chem.MolFromSmiles('CCC'), [])])
def test_get_atom_map_nums(mol, expected_result):
    assert(utils.get_atom_map_nums(mol) == expected_result)


@pytest.mark.parametrize('mol', [(Chem.MolFromSmiles('C[CH2:1]C')), (Chem.MolFromSmiles('CCC'))])
def test_clear_atom_map_nums(mol):
    utils.clear_atom_map_nums(mol)
    assert(Chem.MolToSmiles(mol) == 'CCC')
    assert(not utils.has_atom_map_nums(mol))


@pytest.mark.parametrize('mol', [(Chem.MolFromSmiles('C[13CH2]C')), (Chem.MolFromSmiles('CCC'))])
def test_clear_isotopes(mol):
    utils.clear_isotopes(mol)
    assert(Chem.MolToSmiles(mol) == 'CCC')
    for atom in mol.GetAtoms():
        assert(atom.GetIsotope() == 0)


@pytest.mark.parametrize('mol,map_num,expected_results', [(Chem.MolFromSmiles('[NH2:1][CH2:2]C[OH:3]'), 1, ('N', 2)), (Chem.MolFromSmiles('[NH2:1][CH2:2]C[OH:3]'), 2, ('C', 2)), (Chem.MolFromSmiles('[NH2:1][CH2:2]C[OH:3]'), 3, ('O', 1))])
def get_atom_with_map_num(mol, map_num, expected_results):
    atom = utils.get_atom_with_map_num(mol, map_num)

    symbol, num_hydrogens = expected_results

    assert(atom.GetSymbol() == symbol)
    assert(atom.GetTotalNumHs() == num_hydrogens)


@pytest.mark.parametrize('mol,map_num', [(Chem.MolFromSmiles('[NH2:1][CH2:2]C[OH:3]'), 4)])
def get_atom_with_map_num_fail(mol, map_num):
    with pytest.raises(RuntimeError):
        atom = utils.get_atom_with_map_num(mol, map_num)


@pytest.mark.parametrize('mol,map_num,symbol', [(Chem.MolFromSmiles('[NH2:1][CH2:2][OH:3]'), 1, 'N'), (Chem.MolFromSmiles('[NH2:1][CH2:2][OH:3]'), 2, 'C'), (Chem.MolFromSmiles('[NH2:1][CH2:2][OH:3]'), 3, 'O')])
def test_remove_atom(mol, map_num, symbol):
    atom = utils.get_atom_with_map_num(mol, map_num)
    prev_num_atoms = len(mol.GetAtoms())

    mol = utils.remove_atom(mol, atom.GetIdx())

    atoms = mol.GetAtoms()
    assert(len(atoms) == prev_num_atoms - 1)
    for atom in atoms:
        assert(atom.GetSymbol() != symbol)
