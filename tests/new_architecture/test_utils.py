import pytest
from rdkit import Chem
import new_architecture.utils as utils


@pytest.fixture
def string_list_data():
    return ['asdijviu29478', 'dq908asdfghajkw', 'ac90ded']


@pytest.fixture
def string_int_list_data():
    return ['0', '10', '12']


@pytest.fixture
def empty_list_data():
    return []


def test_get_maximum_len(string_list_data):
    assert(utils.get_maximum(string_list_data, len) == 15)


def test_get_maximum_int(string_int_list_data):
    assert(utils.get_maximum(string_int_list_data, int) == 12)


@pytest.mark.parametrize('data', [({'A': 1}), ([{'B': 2}])])
def test_to_list(data):
    ls_data = utils.to_list(data)
    assert(isinstance(ls_data, list))


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
