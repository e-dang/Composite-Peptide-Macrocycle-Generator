import pytest
import new_architecture.models as models
from rdkit import Chem
from tests.new_architecture.data.mols import *


@pytest.fixture()
def template_from_mol():
    return models.Template.from_mol(Chem.MolFromSmiles('CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'))


@pytest.fixture()
def template_from_dict():
    return models.Template.from_dict(TEST_TEMPLATE_1, _id='temp1')


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


@pytest.fixture()
def peptide_from_mol():
    return models.Peptide.from_mol(Chem.MolFromSmiles('NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1'), False, [models.Monomer.from_dict(TEST_MONOMER_1, _id='qr398fhiusd'), models.Monomer.from_dict(TEST_MONOMER_2, _id='acjiafuy892'), models.Monomer.from_dict(TEST_MONOMER_3, _id='afsidvjoasd')])


@pytest.fixture()
def peptide_from_dict():
    return models.Peptide.from_dict(TEST_PEPTIDE_1, _id='aefoi249')

# @pytest.fixture()
# def template_peptide_from_mol():
#     return models.TemplatePeptide.from_mol(Chem.MolFromSmiles('NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1'), False, [models.Monomer.from_dict(TEST_MONOMER_1, _id='qr398fhiusd'), models.Monomer.from_dict(TEST_MONOMER_2, _id='acjiafuy892'), models.Monomer.from_dict(TEST_MONOMER_3, _id='afsidvjoasd')])

# @pytest.fixture()
# def template_peptide_from_dict():
#     return models.TemplatePeptide.from_dict(TEST_PEPTIDE_1, _id='aefoi249')


def test_template_from_mol(template_from_mol):
    assert(template_from_mol._id == None)
    assert(template_from_mol.binary != None)
    assert(template_from_mol.kekule == 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1')


def test_template_from_dict(template_from_dict):
    assert(template_from_dict._id == 'temp1')
    assert(template_from_dict.binary != None)
    assert(template_from_dict.kekule == 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1')


def test_template_to_dict(template_from_dict):
    assert(template_from_dict.to_dict() == TEST_TEMPLATE_1)


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


def test_peptide_from_mol(peptide_from_mol):
    assert(peptide_from_mol._id == None)
    assert(peptide_from_mol.binary != None)
    assert(peptide_from_mol.kekule ==
           'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1')
    assert(peptide_from_mol.has_cap == False)
    assert(len(peptide_from_mol.monomers) == 3)
    for monomer in peptide_from_mol.monomers:
        assert(monomer['_id'] != None)
        assert(monomer['is_proline'] != None)


def test_peptide_from_dict(peptide_from_dict):
    assert(peptide_from_dict._id == 'aefoi249')
    assert(peptide_from_dict.binary != None)
    assert(peptide_from_dict.kekule ==
           'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1')
    assert(peptide_from_dict.has_cap == False)
    assert(len(peptide_from_dict.monomers) == 3)
    for monomer in peptide_from_dict.monomers:
        assert(monomer['_id'] != None)
        assert(monomer['is_proline'] != None)


def test_peptide_to_dict(peptide_from_dict):
    assert(peptide_from_dict.to_dict() == TEST_PEPTIDE_1)


# def test_template_peptide_from_mol(template_peptide_from_mol):
#     assert(template_peptide_from_mol._id == None)
#     assert(template_peptide_from_mol.binary != None)
#     assert(template_peptide_from_mol.kekule ==
#            'C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=CC3=C(O2)SC=C3)CC(=O)NC(CC2=CSC3=CC=CN23)C(=O)NC(CO)C(=O)N2C[C@@H]3C[C@@H]3[C@H]2C(=O)NCC(CC2=C[NH]C3=CC=CN23)C(=O)O)=C1')
#     assert(template_peptide_from_mol.template == 'temp1')
#     assert(template_peptide_from_mol.peptide == '')
