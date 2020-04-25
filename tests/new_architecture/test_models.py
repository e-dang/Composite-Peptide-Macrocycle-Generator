import pytest
import new_architecture.models as models
from rdkit import Chem
from tests.new_architecture.data.mols import *


@pytest.fixture()
def sidechain_from_mol():
    return models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2'), 'methyl', 'q')


@pytest.fixture()
def sidechain_from_dict():
    return models.Sidechain.from_dict(TEST_SIDECHAIN_3, _id='12498dfhjwu')


@pytest.fixture()
def monomer_from_mol():
    return models.Monomer.from_mol(Chem.MolFromSmiles('NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O'), 'beta3', models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC=CC2=N[NH]C(=O)N12'), 'methyl', 'af'))


@pytest.fixture()
def monomer_from_dict():
    return models.Monomer.from_dict(TEST_MONOMER_1, _id='qr398fhiusd')


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
    assert(sidechain_from_dict.to_dict() == TEST_SIDECHAIN_3)


def test_monomer_from_mol(monomer_from_mol):
    assert(monomer_from_mol._id == None)
    assert(monomer_from_mol.binary != None)
    assert(monomer_from_mol.kekule == 'NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O')
    assert(monomer_from_mol.backbone == 'beta3')
    assert(monomer_from_mol.sidechain == 'af')
    assert(monomer_from_mol.connection == 'methyl')
    assert(monomer_from_mol.is_proline == False)
    assert(monomer_from_mol.imported == False)


def test_monomer_from_dict(monomer_from_dict):
    assert(monomer_from_dict._id == 'qr398fhiusd')
    assert(monomer_from_dict.binary != None)
    assert(monomer_from_dict.kekule == 'O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1')
    assert(monomer_from_dict.backbone == 'alpha')
    assert(monomer_from_dict.sidechain == None)
    assert(monomer_from_dict.connection == None)
    assert(monomer_from_dict.is_proline == True)
    assert(monomer_from_dict.imported == True)


def test_monomer_to_dict(monomer_from_dict):
    assert(monomer_from_dict.to_dict() == TEST_MONOMER_1)
