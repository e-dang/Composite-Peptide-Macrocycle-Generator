import pytest
import new_architecture.models as models
import new_architecture.utils as utils
from rdkit import Chem
from tests.new_architecture.data.mols import *
import macrocycles.exceptions as exceptions


@pytest.fixture()
def backbone_from_mol():
    return models.Backbone.from_mol(Chem.MolFromSmiles('N[CH2:1]C(=O)O'))


@pytest.fixture()
def backbone_from_dict():
    return models.Backbone.from_dict(TEST_BACKBONE_1, _id='alpha')


@pytest.fixture()
def connection_from_mol():
    return models.Connection.from_mol(Chem.MolFromSmiles('CC'))


@pytest.fixture()
def connection_from_dict():
    return models.Connection.from_dict(TEST_CONNECTION_1, _id='ethyl')


@pytest.fixture()
def template_from_mol():
    return models.Template.from_mol(Chem.MolFromSmiles('C/C=C/C1=CC(CC[CH:1]=O)=CC=C1'), 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1')


@pytest.fixture()
def template_from_dict():
    return models.Template.from_dict(TEST_TEMPLATE_1, _id='temp1')


@pytest.fixture()
def sidechain_from_mol():
    return models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2'), models.Connection.from_dict(TEST_CONNECTION_1, _id='ethyl'), 'q')


@pytest.fixture()
def sidechain_from_dict():
    return models.Sidechain.from_dict(TEST_SIDECHAIN_3, _id='12498dfhjwu')


@pytest.fixture()
def monomer_from_mol():
    return models.Monomer.from_mol(Chem.MolFromSmiles('NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O'), models.Backbone.from_dict(TEST_BACKBONE_3, _id='beta3'), models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC=CC2=N[NH]C(=O)N12'), models.Connection.from_dict(TEST_CONNECTION_1, _id='ethyl'), 'af'))


@pytest.fixture()
def monomer_from_dict():
    return models.Monomer.from_dict(TEST_MONOMER_1, _id='qr398fhiusd')


@pytest.fixture()
def peptide_from_mol():
    return models.Peptide.from_mol(Chem.MolFromSmiles('NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1'), False, [models.Monomer.from_dict(TEST_MONOMER_1, _id='qr398fhiusd'), models.Monomer.from_dict(TEST_MONOMER_2, _id='acjiafuy892'), models.Monomer.from_dict(TEST_MONOMER_3, _id='afsidvjoasd')])


@pytest.fixture()
def peptide_from_dict():
    return models.Peptide.from_dict(TEST_PEPTIDE_1, _id='aefoi249')


@pytest.fixture()
def template_peptide_from_mol():
    return models.TemplatePeptide.from_mol(Chem.MolFromSmiles('C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1'), models.Template.from_dict(TEST_TEMPLATE_1, _id='temp1'), models.Peptide.from_dict(TEST_PEPTIDE_1, _id='aefoi249'))


@pytest.fixture()
def template_peptide_from_dict():
    return models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_1, _id='afji923')


# @pytest.fixture()
# def macrocycle_from_mol():
    # return models.Macrocycle.from_mol(Chem.MolFromSmiles('C#CC[C@@]12CCCN1C(=O)[C@@H](CC1=CC=CC3=CC=CN13)NC(=O)[C@@H](CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CC[C@@H](CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O'), 'm', models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_1, _id='afji923'), mode)

def test_backbone_from_mol(backbone_from_mol):
    assert(backbone_from_mol._id == None)
    assert(backbone_from_mol.binary != None)
    assert(backbone_from_mol.kekule == 'NCC(=O)O')
    assert(backbone_from_mol.mapped_kekule == 'N[CH2:1]C(=O)O')


def test_backbone_from_dict(backbone_from_dict):
    assert(backbone_from_dict._id == 'alpha')
    assert(backbone_from_dict.binary != None)
    assert(backbone_from_dict.kekule == 'NCC(=O)O')
    assert(backbone_from_dict.mapped_kekule == 'N[CH2:1]C(=O)O')


def test_backbone_to_dict(backbone_from_dict):
    assert(backbone_from_dict.to_dict() == TEST_BACKBONE_1)


@pytest.mark.parametrize('backbone', [(Chem.MolFromSmiles('N[CH2:1]C(=O)O')), (Chem.MolFromSmiles('N[CH2:1]CC(=O)O')), (Chem.MolFromSmiles('NC[CH2:1]C(=O)O'))])
def test_backbone_validate(backbone):
    assert(models.Backbone.validate(backbone))


@pytest.mark.parametrize('backbone', [(Chem.MolFromSmiles('NCC(=O)O')), (Chem.MolFromSmiles('NCCC(=O)O'))])
def test_validate_backbone_fail(backbone):
    with pytest.raises(exceptions.InvalidMolecule):
        models.Backbone.validate(backbone)


def test_backbone_mol(backbone_from_dict):
    backbone = backbone_from_dict.mol
    assert(models.Backbone.validate(backbone))


def test_backbone_eq(backbone_from_dict):
    assert(backbone_from_dict == models.Backbone.from_dict(TEST_BACKBONE_1, _id='shouldnt matter'))
    assert(backbone_from_dict != models.Backbone.from_dict(TEST_BACKBONE_2, _id='shouldnt matter'))


def test_connection_from_mol(connection_from_mol):
    assert(connection_from_mol._id == None)
    assert(connection_from_mol.binary != None)
    assert(connection_from_mol.kekule == 'CC')


def test_connection_from_dict(connection_from_dict):
    assert(connection_from_dict._id == 'ethyl')
    assert(connection_from_dict.binary != None)
    assert(connection_from_dict.kekule == 'C')


def test_connection_to_dict(connection_from_dict):
    assert(connection_from_dict.to_dict() == TEST_CONNECTION_1)


@pytest.mark.parametrize('template', [(Chem.MolFromSmiles('C/C=C/C1=CC(CC[CH:1]=O)=CC=C1')), (Chem.MolFromSmiles('C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)[CH:1]=O)=C1')), (Chem.MolFromSmiles('C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/C)=C1)C[CH:1]=O'))])
def test_template_validate(template):
    assert(models.Template.validate(template))


@pytest.mark.parametrize('template', [(Chem.MolFromSmiles('CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1')), (Chem.MolFromSmiles('CC(C)(C)OC(=O)OC/C=C/C1=CC(C[C@@H](CC=O)C(=O)ON2C(=O)CCC2=O)=C(F)C=C1'))])
def test_validate_template_fail(template):
    with pytest.raises(exceptions.InvalidMolecule):
        models.Template.validate(template)


def test_template_from_mol(template_from_mol):
    assert(template_from_mol._id == None)
    assert(template_from_mol.binary != None)
    assert(template_from_mol.kekule == 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1')
    assert(template_from_mol.mapped_kekule == 'C/C=C/C1=CC(CC[CH:1]=O)=CC=C1')


def test_template_from_dict(template_from_dict):
    assert(template_from_dict._id == 'temp1')
    assert(template_from_dict.binary != None)
    assert(template_from_dict.kekule == 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1')
    assert(template_from_dict.mapped_kekule == 'C/C=C/C1=CC(CC[CH:1]=O)=CC=C1')


def test_template_to_dict(template_from_dict):
    assert(template_from_dict.to_dict() == TEST_TEMPLATE_1)


def test_sidechain_from_mol(sidechain_from_mol):
    assert(sidechain_from_mol._id == None)
    assert(sidechain_from_mol.binary != None)
    assert(sidechain_from_mol.kekule == 'CC1=CC(=O)C2=C([NH]1)SC=C2')
    assert(sidechain_from_mol.connection == 'ethyl')
    assert(sidechain_from_mol.shared_id == 'q')
    assert(isinstance(sidechain_from_mol.attachment_point, int))

    sidechain = sidechain_from_mol.mol
    attachment_point = sidechain.GetAtomWithIdx(sidechain_from_mol.attachment_point)
    assert(attachment_point.GetSymbol() == 'C')
    assert(attachment_point.GetTotalNumHs() == 3)


def test_sidechain_from_dict(sidechain_from_dict):
    assert(sidechain_from_dict._id == '12498dfhjwu')
    assert(sidechain_from_dict.binary != None)
    assert(sidechain_from_dict.kekule == 'CC1=CC(=O)C2=C([NH]1)SC=C2')
    assert(sidechain_from_dict.connection == 'methyl')
    assert(sidechain_from_dict.shared_id == 'q')
    assert(isinstance(sidechain_from_dict.attachment_point, int))

    sidechain = sidechain_from_dict.mol
    attachment_point = sidechain.GetAtomWithIdx(sidechain_from_dict.attachment_point)
    assert(attachment_point.GetSymbol() == 'C')
    assert(attachment_point.GetTotalNumHs() == 3)


def test_sidechain_to_dict(sidechain_from_dict):
    assert(sidechain_from_dict.to_dict() == TEST_SIDECHAIN_3)


def test_sidechain_eq(sidechain_from_dict):
    assert(sidechain_from_dict == models.Sidechain.from_dict(TEST_SIDECHAIN_3, _id='shouldnt matter'))
    assert(sidechain_from_dict != models.Sidechain.from_dict(TEST_SIDECHAIN_2, _id='shouldnt matter'))


def test_sidechain_mapped_mol(sidechain_from_dict):
    sidechain = sidechain_from_dict.mapped_mol
    atom_idx, map_num = zip(*utils.get_atom_map_nums(sidechain))
    assert(len(atom_idx) == len(map_num) == 1)
    assert(atom_idx[0] == sidechain_from_dict.attachment_point)
    assert(map_num[0] == models.Sidechain.MAP_NUM)


@pytest.mark.parametrize('mol,expected_result', [(Chem.MolFromSmiles('N=C(N)NCCCC(N)C(=O)O'), False), (Chem.MolFromSmiles('O=C(O)[C@@H]1C[C@H](O)CN1'), False), (Chem.MolFromSmiles('NC(CC(=O)O)CC1=CC=CO1'), True)])
def test_monomer_is_required(mol, expected_result):
    assert(models.Monomer.is_required(mol) == expected_result)


@pytest.mark.parametrize('mol,expected_result', [(Chem.MolFromSmiles('N=C(N)NCCCC(N)C(=O)O'), False), (Chem.MolFromSmiles('O=C(O)[C@@H]1C[C@H](O)CN1'), True), (Chem.MolFromSmiles('NC(CC(=O)O)CC1=CC=CO1'), False)])
def test_monomer_is_proline(mol, expected_result):
    assert(models.Monomer.is_proline(mol) == expected_result)


def test_monomer_from_mol(monomer_from_mol):
    assert(monomer_from_mol._id == None)
    assert(monomer_from_mol.binary != None)
    assert(monomer_from_mol.kekule == 'NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O')
    assert(monomer_from_mol.backbone['_id'] == 'beta3')
    assert(monomer_from_mol.backbone['kekule'] == 'NCCC(=O)O')
    assert(monomer_from_mol.sidechain == 'af')
    assert(monomer_from_mol.connection == 'ethyl')
    assert(monomer_from_mol.proline == False)
    assert(monomer_from_mol.imported == False)


def test_monomer_from_dict(monomer_from_dict):
    assert(monomer_from_dict._id == 'qr398fhiusd')
    assert(monomer_from_dict.binary != None)
    assert(monomer_from_dict.kekule == 'O=C(O)[C@@H]1C[C@H](OC2=CC=CC=C2)CN1')
    assert(monomer_from_dict.backbone['_id'] == 'alpha')
    assert(monomer_from_dict.backbone['kekule'] == 'NCC(=O)O')
    assert(monomer_from_dict.sidechain == None)
    assert(monomer_from_dict.connection == None)
    assert(monomer_from_dict.proline == True)
    assert(monomer_from_dict.imported == True)


def test_monomer_to_dict(monomer_from_dict):
    assert(monomer_from_dict.to_dict() == TEST_MONOMER_1)


def test_monomer_eq(monomer_from_dict):
    assert(monomer_from_dict == models.Monomer.from_dict(TEST_MONOMER_1, _id='shouldnt matter'))
    assert(monomer_from_dict != models.Monomer.from_dict(TEST_MONOMER_2, _id='shouldnt matter'))


def test_peptide_from_mol(peptide_from_mol):
    assert(peptide_from_mol._id == None)
    assert(peptide_from_mol.binary != None)
    assert(peptide_from_mol.kekule ==
           'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1')
    assert(peptide_from_mol.has_cap == False)
    assert(len(peptide_from_mol.monomers) == 3)
    for monomer in peptide_from_mol.monomers:
        assert(monomer['_id'] != None)
        assert(monomer['proline'] != None)
        with pytest.raises(KeyError):
            monomer['binary']


def test_peptide_from_dict(peptide_from_dict):
    assert(peptide_from_dict._id == 'aefoi249')
    assert(peptide_from_dict.binary != None)
    assert(peptide_from_dict.kekule ==
           'NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1')
    assert(peptide_from_dict.has_cap == False)
    assert(len(peptide_from_dict.monomers) == 3)
    for monomer in peptide_from_dict.monomers:
        assert(monomer['_id'] != None)
        assert(monomer['proline'] != None)


def test_peptide_to_dict(peptide_from_dict):
    assert(peptide_from_dict.to_dict() == TEST_PEPTIDE_1)


def test_peptide_eq(peptide_from_dict):
    assert(peptide_from_dict == models.Peptide.from_dict(TEST_PEPTIDE_1, _id='shouldnt matter'))
    assert(peptide_from_dict != models.Peptide.from_dict(TEST_PEPTIDE_2, _id='shouldnt matter'))


def test_template_peptide_from_mol(template_peptide_from_mol):
    assert(template_peptide_from_mol._id == None)
    assert(template_peptide_from_mol.binary != None)
    assert(template_peptide_from_mol.kekule ==
           'C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1')
    assert(template_peptide_from_mol.template == 'temp1')
    assert(template_peptide_from_mol.peptide['_id'] == 'aefoi249')
    assert(template_peptide_from_mol.peptide['has_cap'] == False)
    assert(len(template_peptide_from_mol.peptide['monomers']) == 3)
    with pytest.raises(KeyError):
        template_peptide_from_mol.peptide['binary']


def test_template_peptide_from_dict(template_peptide_from_dict):
    assert(template_peptide_from_dict._id == 'afji923')
    assert(template_peptide_from_dict.binary != None)
    assert(template_peptide_from_dict.kekule ==
           'C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1')
    assert(template_peptide_from_dict.template == 'temp1')
    assert(template_peptide_from_dict.peptide['_id'] == 'aefoi249')
    assert(template_peptide_from_dict.peptide['has_cap'] == False)
    assert(len(template_peptide_from_dict.peptide['monomers']) == 3)
    with pytest.raises(KeyError):
        template_peptide_from_dict.peptide['binary']


def test_template_peptide_to_dict(template_peptide_from_dict):
    assert(template_peptide_from_dict.to_dict() == TEST_TEMPLATE_PEPTIDE_1)


# def test_macrocycle_from_mol(macrocycle_from_mol):
#     assert(macrocycle_from_mol._id == None)
#     assert(macrocycle_from_mol.binary != None)
#     assert(macrocycle_from_mol.kekule ==
#            'C#CC[C@@]12CCCN1C(=O)[C@@H](CC1=CC=CC3=CC=CN13)NC(=O)[C@@H](CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CC[C@@H](CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O')
#     assert(macrocycle_from_mol.template == 'temp1')
#     assert(macrocycle_from_mol.template_peptide == '')


# def test_macrocycle_from_dict(macrocycle_from_dict):
#     assert(macrocycle_from_dict._id == 'afji923')
#     assert(macrocycle_from_dict.binary != None)
#     assert(macrocycle_from_dict.kekule ==
#            'C#CC[C@@]12CCCN1C(=O)[C@@H](CC1=CC=CC3=CC=CN13)NC(=O)[C@@H](CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CC[C@@H](CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O')
#     assert(macrocycle_from_dict.template == 'temp1')
#     assert(macrocycle_from_dict.peptide['_id'] == 'aefoi249')
#     assert(macrocycle_from_dict.peptide['has_cap'] == False)
#     assert(len(macrocycle_from_dict.peptide['monomers']) == 3)
#     with pytest.raises(KeyError):
#         macrocycle_from_dict.peptide['binary']


# def test_macrocycle_to_dict(macrocycle_from_dict):
#     assert(macrocycle_from_dict.to_dict() == TEST_MACROCYCLE_1)
