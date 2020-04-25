import pytest
from rdkit import Chem
import new_architecture.utils as utils


@pytest.mark.parametrize('mol,expected_result', [(Chem.MolFromSmiles('C[CH2:1]C'), True), (Chem.MolFromSmiles('CCC'), False)])
def test_has_atom_map_nums(mol, expected_result):
    assert(utils.has_atom_map_nums(mol) == expected_result)


@pytest.mark.parametrize('mol', [(Chem.MolFromSmiles('C[CH2:1]C')), (Chem.MolFromSmiles('CCC'))])
def test_clear_atom_map_nums(mol):
    utils.clear_atom_map_nums(mol)
    assert(Chem.MolToSmiles(mol) == 'CCC')
    assert(not utils.has_atom_map_nums(mol))
