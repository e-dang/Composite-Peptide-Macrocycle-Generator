import pytest
import new_architecture.models as models
from rdkit import Chem

SIDECHAIN_DICT = {'binary': Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2').ToBinary(
), 'kekule': 'CC1=CC(=O)C2=C([NH]1)SC=C2', 'connection': 'methyl', 'shared_id': 'q'}


@pytest.fixture()
def sidechain_from_mol():
    return models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2'), 'methyl', 'q')


@pytest.fixture()
def sidechain_from_dict():
    return models.Sidechain.from_dict(SIDECHAIN_DICT, _id='12498dfhjwu')


def test_sidechain_from_mol(sidechain_from_mol):
    assert(sidechain_from_mol._id == None)
    assert(sidechain_from_mol.binary != None)
    assert(sidechain_from_mol.kekule == 'CC1=CC(=O)C2=C([NH]1)SC=C2')
    assert(sidechain_from_mol.connection == 'methyl')
    assert(sidechain_from_mol.shared_id == 'q')


def test_sidechain_from_dict(sidechain_from_dict):
    assert(sidechain_from_dict._id == '12498dfhjwu')
    assert(sidechain_from_dict.binary != None)
    assert(sidechain_from_dict.kekule == 'CC1=CC(=O)C2=C([NH]1)SC=C2')
    assert(sidechain_from_dict.connection == 'methyl')
    assert(sidechain_from_dict.shared_id == 'q')


def test_sidechain_to_dict(sidechain_from_dict):
    assert(sidechain_from_dict.to_dict() == SIDECHAIN_DICT)
