import uuid

import numpy as np
import pytest
from rdkit import Chem

import cpmg.exceptions as exceptions
import cpmg.models as models
import cpmg.utils as utils
from data.mols import *


@pytest.fixture(params=['N[CH2:1]C(=O)O', 'N[CH2:1]CC(=O)O', 'NC[CH2:1]C(=O)O'])
def backbone_kekule(request):
    return request.param


@pytest.fixture(params=['C', 'CC'])
def connection_kekule(request):
    return request.param


@pytest.fixture()
def backbone_from_mol(backbone_kekule):
    mol = Chem.MolFromSmiles(backbone_kekule)
    backbone = models.Backbone.from_mol(mol)
    utils.clear_atom_map_nums(mol)
    Chem.Kekulize(mol)
    kekule = Chem.MolToSmiles(mol, kekuleSmiles=True)
    binary = mol.ToBinary()
    return backbone, binary, backbone_kekule, kekule


@pytest.fixture(params=TEST_BACKBONES)
def backbone_from_dict(request):
    _id = str(uuid.uuid4())
    return models.Backbone.from_dict(request.param, _id=_id), request.param, _id


@pytest.fixture()
def connection_from_mol(connection_kekule):
    mol = Chem.MolFromSmiles(connection_kekule)
    connection = models.Connection.from_mol(mol)
    Chem.Kekulize(mol)
    binary = mol.ToBinary()
    return connection, binary, connection_kekule


@pytest.fixture(params=TEST_CONNECTIONS)
def connection_from_dict(request):
    _id = str(uuid.uuid4())
    return models.Connection.from_dict(request.param, _id=_id), request.param, _id


@pytest.fixture()
def template_from_mol():
    kekule = 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'
    oligomerization_kekule = f'C/C=C/C1=CC(CC[CH:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O)=CC=C1'
    friedel_crafts_kekule = f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]'
    tsuji_trost_kekule = f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]'
    pictet_spangler_kekule = f'[O:{models.Template.PS_OXYGEN_MAP_NUM}]=[CH1:{models.Template.PS_CARBON_MAP_NUM}]C[C@H](C[*:{models.Template.WC_MAP_NUM_1}])[CH1:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O'
    template_pictet_spangler_kekule = f'C#CCCC[C@@]([CH1:{models.Template.PS_CARBON_MAP_NUM}]=[O:{models.Template.PS_OXYGEN_MAP_NUM}])(CC(=O)[NH1:{models.Template.TEMPLATE_PS_NITROGEN_MAP_NUM}][*:{models.Template.WC_MAP_NUM_1}])CC1=[CH1:{models.Template.EAS_MAP_NUM}]C=CC([*:{models.Template.WC_MAP_NUM_2}])=C1'
    pyrroloindoline_kekule = f'[*:{models.Template.WC_MAP_NUM_1}]/C=C/[CH3:{models.Template.EAS_MAP_NUM}]'
    aldehyde_cyclization_kekule = f'[O:{models.Template.PS_OXYGEN_MAP_NUM}]=[CH1:{models.Template.PS_CARBON_MAP_NUM}]C[C@H](C[*:{models.Template.WC_MAP_NUM_1}])[CH1:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O'
    template = models.Template.from_mol(Chem.MolFromSmiles(
        kekule), oligomerization_kekule, friedel_crafts_kekule, tsuji_trost_kekule, pictet_spangler_kekule,
        template_pictet_spangler_kekule, pyrroloindoline_kekule, aldehyde_cyclization_kekule)
    return (template, kekule, oligomerization_kekule, friedel_crafts_kekule, tsuji_trost_kekule, pictet_spangler_kekule,
            template_pictet_spangler_kekule, pyrroloindoline_kekule, aldehyde_cyclization_kekule)


@pytest.fixture(params=[TEST_TEMPLATE_1, TEST_TEMPLATE_2, TEST_TEMPLATE_3])
def template_from_dict(request):
    _id = str(uuid.uuid4())
    template = models.Template.from_dict(request.param, _id=_id)
    return template, request.param, _id


@pytest.fixture()
def sidechain_from_mol():
    kekule = 'CC1=CC(=O)C2=C([NH]1)SC=C2'
    mol = Chem.MolFromSmiles(kekule)
    connection_id = str(uuid.uuid4())
    shared_id = str(uuid.uuid4())
    connection = models.Connection.from_dict(TEST_CONNECTION_1, _id=connection_id)
    sidechain = models.Sidechain.from_mol(mol, connection, shared_id)
    Chem.Kekulize(mol)
    return sidechain, mol.ToBinary(), kekule, connection, shared_id


@pytest.fixture(params=TEST_SIDECHAINS)
def sidechain_from_dict(request):
    _id = str(uuid.uuid4())
    return models.Sidechain.from_dict(request.param, _id=_id), request.param, _id


@pytest.fixture()
def monomer_from_mol():
    kekule = 'NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O'
    mol = Chem.MolFromSmiles(kekule)
    backbone_id = str(uuid.uuid4())
    backbone = models.Backbone.from_dict(TEST_BACKBONE_3, _id=backbone_id)
    connection = models.Connection.from_dict(TEST_CONNECTION_1, _id=str(uuid.uuid4()))
    sidechain_shared_id = str(uuid.uuid4())
    sidechain = models.Sidechain.from_mol(Chem.MolFromSmiles(
        'CC1=CC=CC2=N[NH]C(=O)N12'), connection, sidechain_shared_id)
    monomer = models.Monomer.from_mol(mol, backbone, sidechain)
    Chem.Kekulize(mol)
    return monomer, mol.ToBinary(), kekule, backbone, sidechain, connection


@pytest.fixture(params=TEST_MONOMERS_W_IDXS)
def monomer_from_dict(request):
    _id = str(uuid.uuid4())
    return models.Monomer.from_dict(request.param, _id=_id), request.param, _id


@pytest.fixture()
def peptide_from_mol():
    return models.Peptide.from_mol(Chem.MolFromSmiles('NCC(=O)NC(CC(=O)NC(CC(=O)NCC(CC1=CC2=C(OC=C2)S1)C(=O)NC(CC1=CC=C2C=CC=CC=C21)C(=O)O)CC1=CC=CO1)CC1=C2C=CSC2=NS1'), 3, False, [models.Monomer.from_dict(TEST_MONOMER_1, _id='qr398fhiusd'), models.Monomer.from_dict(TEST_MONOMER_2, _id='acjiafuy892'), models.Monomer.from_dict(TEST_MONOMER_3, _id='afsidvjoasd')])


@pytest.fixture()
def peptide_from_dict():
    return models.Peptide.from_dict(TEST_PEPTIDE_1, _id='aefoi249')


@pytest.fixture()
def template_peptide_from_mol():
    return models.TemplatePeptide.from_mol(Chem.MolFromSmiles('C/C=C/C1=CC=CC(CCC(=O)NC(CCC2=C3C=COC3=CS2)C(=O)NC(CCC2=CC=CC3=COC=C23)C(=O)NC(CCC2=CC=C(O)C=C2)C(=O)N2[C@H](C(=O)NC(C(=O)O)C(C)C)C[C@H]3C[C@H]32)=C1'), models.Template.from_dict(TEST_TEMPLATE_1, _id='temp1'), models.Peptide.from_dict(TEST_PEPTIDE_1, _id='aefoi249'))


@pytest.fixture()
def template_peptide_from_dict():
    return models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_1, _id='afji923')


@pytest.fixture()
def macrocycle_from_mol():
    tp_id = str(uuid.uuid4())
    rxn_id = str(uuid.uuid4())
    kekule = 'C#CC[C@@]12CCCN1C(=O)[C@@H](CC1=CC=CC3=CC=CN13)NC(=O)[C@@H](CC(C)C)N1C=C[C@H](CC3=C(F)C=CC(=C3)/C=C/CC3=C4C=CC=CC4=CN3CC[C@@H](CC(=O)O)NC(=O)[C@]34CC[C@H](C3)N4C2=O)C1=O'
    macrocycle = models.Macrocycle.from_mol(Chem.MolFromSmiles(kekule), '', models.TemplatePeptide.from_dict(
        TEST_TEMPLATE_PEPTIDE_1, _id=tp_id), [models.Reaction.from_dict(TEST_REACTION_1, _id=rxn_id)])
    return macrocycle, kekule, TEST_TEMPLATE_PEPTIDE_1['template'], tp_id, rxn_id, TEST_REACTION_1['type'], TEST_TEMPLATE_PEPTIDE_1['peptide']['has_cap']


@pytest.fixture()
def macrocycle_from_dict():
    _id = str(uuid.uuid4())
    macrocycle = models.Macrocycle.from_dict(TEST_MACROCYCLE_1, _id=_id)
    return macrocycle, TEST_MACROCYCLE_1, _id


@pytest.fixture()
def reaction_from_mols():
    rxn_type = 'tsuji_trost'
    smarts = '(c1c([*:3])ccc([OH:4])c1.C(=C/[*:50])\\[CH3:2])>>C(=C/[*:50])\\[CH2:2][O:4]c1ccc([*:3])cc1'
    template = models.Template.from_dict(TEST_TEMPLATE_1, _id=str(uuid.uuid4()))
    reacting_mol = models.Sidechain.from_dict(TEST_SIDECHAIN_1, _id=str(uuid.uuid4()))
    rxn_atom_idx = 5
    reaction = models.Reaction.from_mols(rxn_type, smarts, template, reacting_mol, rxn_atom_idx)
    return reaction, rxn_type, smarts, template, reacting_mol, rxn_atom_idx


@pytest.fixture(params=[TEST_REACTION_1, TEST_REACTION_2])
def reaction_from_dict(request):
    _id = str(uuid.uuid4())
    reaction = models.Reaction.from_dict(request.param, _id=_id)
    return reaction, request.param, _id


@pytest.fixture(params=[TEST_REGIOSQM_PREDICTION_1, TEST_REGIOSQM_PREDICTION_2])
def regiosqm_prediction_from_dict(request):
    _id = str(uuid.uuid4())
    prediction = models.RegioSQMPrediction.from_dict(request.param, _id=_id)
    return prediction, request.param, _id


@pytest.fixture(params=[TEST_PKA_PREDICTION_1, TEST_PKA_PREDICTION_2, TEST_PKA_PREDICTION_3])
def pka_prediction_from_dict(request):
    _id = str(uuid.uuid4())
    prediction = models.pKaPrediction.from_dict(request.param, _id=_id)
    return prediction, request.param, _id


@pytest.fixture()
def peptide_plan_data():
    return [(1, 2, 3), (3, 4, 5, 6)], 3


def test_backbone_from_mol(backbone_from_mol):
    backbone, binary, mapped_kekule, kekule = backbone_from_mol
    assert backbone._id == None
    assert backbone.binary == binary
    assert backbone.kekule == kekule
    assert backbone.mapped_kekule == mapped_kekule


def test_backbone_from_dict(backbone_from_dict):
    backbone, backbone_dict, _id = backbone_from_dict
    assert backbone._id == _id
    assert backbone.binary == backbone_dict['binary']
    assert backbone.kekule == backbone_dict['kekule']
    assert backbone.mapped_kekule == backbone_dict['mapped_kekule']


def test_backbone_to_dict(backbone_from_dict):
    backbone, backbone_dict, _id = backbone_from_dict
    assert backbone.to_dict() == backbone_dict


def test_backbone_validate(backbone_kekule):
    assert models.Backbone.validate(Chem.MolFromSmiles(backbone_kekule))


@pytest.mark.parametrize('backbone', [
    (Chem.MolFromSmiles('NCC(=O)O')),
    (Chem.MolFromSmiles('NCCC(=O)O'))
])
def test_validate_backbone_fail(backbone):
    with pytest.raises(exceptions.InvalidMolecule):
        models.Backbone.validate(backbone)


def test_backbone_mol(backbone_from_dict):
    backbone, _, _ = backbone_from_dict
    backbone = backbone.mol
    assert models.Backbone.validate(backbone)


def test_backbone_eq(backbone_from_dict):
    backbone, backbone_dict, _ = backbone_from_dict

    backbones = deepcopy(TEST_BACKBONES)
    backbones.remove(backbone_dict)

    assert backbone == models.Backbone.from_dict(backbone_dict, _id='shouldnt matter')
    assert backbone != models.Backbone.from_dict(backbones[0], _id='shouldnt matter')


def test_connection_from_mol(connection_from_mol):
    connection, binary, kekule = connection_from_mol
    assert connection._id == None
    assert connection.binary == binary
    assert connection.kekule == kekule


def test_connection_from_dict(connection_from_dict):
    connection, connection_dict, _id = connection_from_dict
    assert connection._id == _id
    assert connection.binary == connection_dict['binary']
    assert connection.kekule == connection_dict['kekule']


def test_connection_to_dict(connection_from_dict):
    connection, connection_dict, _ = connection_from_dict
    assert connection.to_dict() == connection_dict


@pytest.mark.parametrize('func,template', [(models.Template.validate_oligomerization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['oligomerization_kekule'])),
                                           (models.Template.validate_oligomerization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['oligomerization_kekule'])),
                                           (models.Template.validate_oligomerization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['oligomerization_kekule'])),
                                           (models.Template.validate_friedel_crafts_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['friedel_crafts_kekule'])),
                                           (models.Template.validate_friedel_crafts_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['friedel_crafts_kekule'])),
                                           (models.Template.validate_friedel_crafts_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['friedel_crafts_kekule'])),
                                           (models.Template.validate_tsuji_trost_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['tsuji_trost_kekule'])),
                                           (models.Template.validate_tsuji_trost_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['tsuji_trost_kekule'])),
                                           (models.Template.validate_tsuji_trost_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['tsuji_trost_kekule'])),
                                           (models.Template.validate_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['pictet_spangler_kekule'])),
                                           (models.Template.validate_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['pictet_spangler_kekule'])),
                                           (models.Template.validate_template_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['template_pictet_spangler_kekule'])),
                                           (models.Template.validate_pyrroloindoline_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['pyrroloindoline_kekule'])),
                                           (models.Template.validate_pyrroloindoline_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['pyrroloindoline_kekule'])),
                                           (models.Template.validate_pyrroloindoline_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['pyrroloindoline_kekule'])),
                                           (models.Template.validate_aldehyde_cyclization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['aldehyde_cyclization_kekule']))])
def test_template_validate(func, template):
    assert func(template)


@pytest.mark.parametrize('func,template', [(models.Template.validate_oligomerization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['kekule'])),
                                           (models.Template.validate_oligomerization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['kekule'])),
                                           (models.Template.validate_oligomerization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['kekule'])),
                                           (models.Template.validate_friedel_crafts_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['kekule'])),
                                           (models.Template.validate_friedel_crafts_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['kekule'])),
                                           (models.Template.validate_friedel_crafts_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['kekule'])),
                                           (models.Template.validate_tsuji_trost_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['kekule'])),
                                           (models.Template.validate_tsuji_trost_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['kekule'])),
                                           (models.Template.validate_tsuji_trost_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['kekule'])),
                                           (models.Template.validate_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['kekule'])),
                                           (models.Template.validate_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['kekule'])),
                                           (models.Template.validate_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['kekule'])),
                                           (models.Template.validate_template_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['kekule'])),
                                           (models.Template.validate_template_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['kekule'])),
                                           (models.Template.validate_template_pictet_spangler_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['kekule'])),
                                           (models.Template.validate_pyrroloindoline_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['kekule'])),
                                           (models.Template.validate_pyrroloindoline_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['kekule'])),
                                           (models.Template.validate_pyrroloindoline_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['kekule'])),
                                           (models.Template.validate_aldehyde_cyclization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_1['kekule'])),
                                           (models.Template.validate_aldehyde_cyclization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_2['kekule'])),
                                           (models.Template.validate_aldehyde_cyclization_mol,
                                            Chem.MolFromSmiles(TEST_TEMPLATE_3['kekule']))])
def test_template_validate_fail(func, template):
    with pytest.raises(ValueError):
        func(template)


def test_template_from_mol(template_from_mol):
    (template, kekule, oligomerization_kekule, friedel_crafts_kekule, tsuji_trost_kekule, pictet_spangler_kekule,
     template_pictet_spangler_kekule, pyrroloindoline_kekule, aldehyde_cyclization_kekule) = template_from_mol
    assert template._id == None
    assert template.binary != None
    assert isinstance(Chem.Mol(template.binary), Chem.Mol)
    assert template.kekule == kekule
    assert template.oligomerization_kekule == oligomerization_kekule
    assert template.friedel_crafts_kekule == friedel_crafts_kekule
    assert template.tsuji_trost_kekule == tsuji_trost_kekule
    assert template.pictet_spangler_kekule == pictet_spangler_kekule
    assert template.template_pictet_spangler_kekule == template_pictet_spangler_kekule
    assert template.pyrroloindoline_kekule == pyrroloindoline_kekule
    assert template.aldehyde_cyclization_kekule == aldehyde_cyclization_kekule
    assert len(template.__dict__) - 1 == len(template_from_mol)  # -1 for binary field


def test_template_from_dict(template_from_dict):
    template, template_dict, _id = template_from_dict
    assert template._id == _id
    assert template.binary == template_dict['binary']
    assert isinstance(Chem.Mol(template.binary), Chem.Mol)
    assert template.kekule == template_dict['kekule']
    assert template.oligomerization_kekule == template_dict['oligomerization_kekule']
    assert template.friedel_crafts_kekule == template_dict['friedel_crafts_kekule']
    assert template.tsuji_trost_kekule == template_dict['tsuji_trost_kekule']
    assert template.pictet_spangler_kekule == template_dict['pictet_spangler_kekule']
    assert template.template_pictet_spangler_kekule == template_dict['template_pictet_spangler_kekule']
    assert template.pyrroloindoline_kekule == template_dict['pyrroloindoline_kekule']
    assert template.aldehyde_cyclization_kekule == template_dict['aldehyde_cyclization_kekule']
    assert len(template.__dict__) - 1 == len(template_dict)


def test_template_to_dict(template_from_dict):
    template, template_dict, _ = template_from_dict
    assert template.to_dict() == template_dict


def test_sidechain_from_mol(sidechain_from_mol):
    sidechain, binary, kekule, connection, shared_id = sidechain_from_mol
    sidechain_mol = sidechain.mol

    assert sidechain._id == None
    assert sidechain.binary == binary
    assert sidechain.kekule == kekule
    assert sidechain.connection == connection.kekule
    assert sidechain.shared_id == shared_id
    assert sidechain.attachment_point in list(range(len(sidechain_mol.GetAtoms())))

    attachment_point = sidechain_mol.GetAtomWithIdx(sidechain.attachment_point)
    assert attachment_point.GetSymbol() == 'C'
    assert attachment_point.GetTotalNumHs() == 3


def test_sidechain_from_dict(sidechain_from_dict):
    sidechain, sidechain_dict, _id = sidechain_from_dict
    sidechain_mol = sidechain.mol

    assert sidechain._id == _id
    assert sidechain.binary == sidechain_dict['binary']
    assert sidechain.kekule == sidechain_dict['kekule']
    assert sidechain.connection == sidechain_dict['connection']
    assert sidechain.shared_id == sidechain_dict['shared_id']
    assert sidechain.attachment_point in list(range(len(sidechain_mol.GetAtoms())))

    attachment_point = sidechain_mol.GetAtomWithIdx(sidechain.attachment_point)
    assert attachment_point.GetSymbol() == 'C'
    assert attachment_point.GetTotalNumHs() == 3


@pytest.mark.parametrize('sidechain_kekule', [(doc['kekule']) for doc in TEST_SIDECHAINS])
def test_sidechain_validate(sidechain_kekule):
    sidechain = Chem.MolFromSmiles(sidechain_kekule)
    attachment_point = models.Sidechain.validate(sidechain)
    attachment_atom = sidechain.GetAtomWithIdx(attachment_point)

    assert attachment_atom.GetSymbol() == 'C'
    assert attachment_atom.GetTotalNumHs() == 3


def test_sidechain_to_dict(sidechain_from_dict):
    sidechain, sidechain_dict, _ = sidechain_from_dict
    assert sidechain.to_dict() == sidechain_dict


def test_sidechain_eq(sidechain_from_dict):
    sidechain, sidechain_dict, _ = sidechain_from_dict

    sidechains = deepcopy(TEST_SIDECHAINS)
    sidechains.remove(sidechain_dict)

    assert sidechain == models.Sidechain.from_dict(sidechain_dict, _id='shouldnt matter')
    assert sidechain != models.Sidechain.from_dict(sidechains[0], _id='shouldnt matter')


def test_sidechain_mapped_mol(sidechain_from_dict):
    sidechain, _, _ = sidechain_from_dict

    mapped_mol = sidechain.mapped_mol
    atom_idx, map_num = zip(*utils.get_atom_map_nums(mapped_mol))

    assert len(atom_idx) == len(map_num) == 1
    assert atom_idx[0] == sidechain.attachment_point
    assert map_num[0] == models.Sidechain.MAP_NUM


@pytest.mark.parametrize('mol,expected_result', [(Chem.MolFromSmiles('N=C(N)NCCCC(N)C(=O)O'), False), (Chem.MolFromSmiles('O=C(O)[C@@H]1C[C@H](O)CN1'), False), (Chem.MolFromSmiles('NC(CC(=O)O)CC1=CC=CO1'), True)])
def test_monomer_is_required(mol, expected_result):
    assert(models.Monomer.is_required(mol) == expected_result)


@pytest.mark.parametrize('mol,expected_result', [(Chem.MolFromSmiles('N=C(N)NCCCC(N)C(=O)O'), False), (Chem.MolFromSmiles('O=C(O)[C@@H]1C[C@H](O)CN1'), True), (Chem.MolFromSmiles('NC(CC(=O)O)CC1=CC=CO1'), False)])
def test_monomer_is_proline(mol, expected_result):
    assert(models.Monomer.is_proline(mol) == expected_result)


def test_monomer_from_mol(monomer_from_mol):
    monomer, binary, kekule, backbone, sidechain, connection = monomer_from_mol
    assert monomer._id == None
    assert monomer.binary == binary
    assert monomer.kekule == kekule
    assert monomer.backbone['_id'] == backbone._id
    assert monomer.backbone['kekule'] == backbone.kekule
    assert monomer.sidechain == sidechain.shared_id
    assert monomer.connection == connection.kekule
    assert monomer.proline == False
    assert monomer.imported == False
    assert monomer.index == None


def test_monomer_from_dict(monomer_from_dict):
    monomer, monomer_dict, _id = monomer_from_dict
    assert monomer._id == _id
    assert monomer.binary == monomer_dict['binary']
    assert monomer.kekule == monomer_dict['kekule']
    assert monomer.backbone['_id'] == monomer_dict['backbone']['_id']
    assert monomer.backbone['kekule'] == monomer_dict['backbone']['kekule']
    assert monomer.sidechain == monomer_dict['sidechain']
    assert monomer.connection == monomer_dict['connection']
    assert monomer.index == monomer_dict['index']


def test_monomer_to_dict(monomer_from_dict):
    monomer, monomer_dict, _ = monomer_from_dict
    assert monomer.to_dict() == monomer_dict


def test_monomer_eq(monomer_from_dict):
    monomer, monomer_dict, _id = monomer_from_dict

    monomers = deepcopy(TEST_MONOMERS_W_IDXS)
    monomers.remove(monomer_dict)

    assert monomer == models.Monomer.from_dict(monomer_dict, _id='shouldnt matter')
    assert monomer != models.Monomer.from_dict(monomers[0], _id='shouldnt matter')


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


def test_template_peptide_eq(template_peptide_from_dict):
    assert(template_peptide_from_dict == models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_1, _id='doesnt matter'))
    assert(template_peptide_from_dict != models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_2, _id='doesnt matter'))


@pytest.mark.parametrize('smiles', [('C1CCCCCCCCC1'), ('C1CCCCCCCCCC1')])
def test_macrocycle_validate(smiles):
    assert(models.Macrocycle.validate(Chem.MolFromSmiles(smiles)))


@pytest.mark.parametrize('smiles', [('C1CCCCCCCC1'), ('CCCCCCCCCCC')])
def test_macrocycle_validate_fail(smiles):
    with pytest.raises(exceptions.InvalidMolecule):
        models.Macrocycle.validate(Chem.MolFromSmiles(smiles))


def test_macrocycle_from_mol(macrocycle_from_mol):
    macrocycle, kekule, template_id, tp_id, rxn_id, rxn_type, has_cap = macrocycle_from_mol
    assert(macrocycle._id == None)
    assert(macrocycle.binary != None)
    assert(isinstance(Chem.Mol(macrocycle.binary), Chem.Mol))
    assert(macrocycle.kekule == kekule)
    assert(macrocycle.has_cap == has_cap)
    assert(macrocycle.template == template_id)
    assert(macrocycle.template_peptide == tp_id)
    assert(macrocycle.reactions == [{'_id': rxn_id, 'type': rxn_type}])


def test_macrocycle_from_dict(macrocycle_from_dict):
    macrocycle, macrocycle_dict, _id = macrocycle_from_dict
    assert(macrocycle._id == _id)
    assert(macrocycle.binary != None)
    assert(macrocycle.kekule == macrocycle_dict['kekule'])
    assert(macrocycle.has_cap == macrocycle_dict['has_cap'])
    assert(macrocycle.template == macrocycle_dict['template'])
    assert(macrocycle.template_peptide == macrocycle_dict['template_peptide'])
    assert(macrocycle.reactions == macrocycle_dict['reactions'])


def test_macrocycle_to_dict(macrocycle_from_dict):
    macrocycle, macrocycle_dict, _ = macrocycle_from_dict
    assert(macrocycle.to_dict() == macrocycle_dict)


def test_reaction_from_mols(reaction_from_mols):
    reaction, rxn_type, smarts, template, reacting_mol, rxn_atom_idx = reaction_from_mols
    assert(reaction._id == None)
    assert(reaction.type == rxn_type)
    assert(reaction.smarts == smarts)
    assert(reaction.binary != None)
    assert(isinstance(AllChem.ChemicalReaction(reaction.binary), AllChem.ChemicalReaction))
    assert(reaction.template == template._id)
    assert(isinstance(reaction.reacting_mol, dict) or reaction.reacting_mol is None)
    assert(reaction.reacting_mol == {'_id': reacting_mol.shared_id, 'kekule': reacting_mol.kekule})
    assert(reaction.rxn_atom_idx == rxn_atom_idx)


def test_reaction_from_dict(reaction_from_dict):
    reaction, reaction_dict, _id = reaction_from_dict
    assert(reaction._id == _id)
    assert(reaction.type == reaction_dict['type'])
    assert(reaction.smarts == reaction_dict['smarts'])
    assert(reaction.binary != None)
    assert(isinstance(AllChem.ChemicalReaction(reaction.binary), AllChem.ChemicalReaction))
    assert(reaction.template == reaction_dict['template'])
    assert(isinstance(reaction.reacting_mol, dict) or reaction.reacting_mol is None)
    assert(reaction.reacting_mol == reaction_dict['reacting_mol'])
    assert(reaction.rxn_atom_idx == reaction_dict['rxn_atom_idx'])


def test_reaction_to_dict(reaction_from_dict):
    reaction, reaction_dict, _ = reaction_from_dict
    assert(reaction.to_dict() == reaction_dict)


def test_regiosqm_prediction_from_dict(regiosqm_prediction_from_dict):
    prediction, prediction_dict, _id = regiosqm_prediction_from_dict
    assert(prediction._id == _id)
    assert(prediction.predictions == prediction_dict['predictions'])
    assert(prediction.reacting_mol == prediction_dict['reacting_mol'])
    assert(prediction.solvent == prediction_dict['solvent'])


def test_regiosqm_prediction_to_dict(regiosqm_prediction_from_dict):
    prediction, prediction_dict, _ = regiosqm_prediction_from_dict
    assert(prediction.to_dict() == prediction_dict)


def test_pka_prediction_from_dict(pka_prediction_from_dict):
    prediction, prediction_dict, _id = pka_prediction_from_dict
    assert(prediction._id == _id)
    assert(prediction.predictions == prediction_dict['predictions'])
    assert(prediction.reacting_mol == prediction_dict['reacting_mol'])
    assert(prediction.solvent == prediction_dict['solvent'])


def test_pka_prediction_to_dict(pka_prediction_from_dict):
    prediction, prediction_dict, _ = pka_prediction_from_dict
    assert(prediction.to_dict() == prediction_dict)


@pytest.mark.parametrize('length', [(3), (4), (5)])
def test_peptide_plan_constructor(length):
    plan = models.PeptidePlan(length)

    assert plan.reg_length == length
    assert plan.cap_length == length + 1


@pytest.mark.parametrize('data, peptide_length, expected_len, expected_data', [
    ([(1, 2, 3), (1, 2, 3), (3, 4, 5)], 3, 2, [(1, 2, 3), (3, 4, 5)])
])
def test_peptide_plan_add(data, peptide_length, expected_len, expected_data):
    plan = models.PeptidePlan(peptide_length)

    for tup in data:
        plan.add(tup)

    assert(len(plan) == expected_len)
    assert(plan.combinations == set(expected_data))


@pytest.mark.parametrize('data, peptide_length', [
    ([(1, 2, 3)], 4),
    ([(1, 2, 3, 4, 5)], 3)
])
def test_peptide_plan_add_fail(data, peptide_length):
    plan = models.PeptidePlan(peptide_length)

    with pytest.raises(RuntimeError):
        for tup in data:
            plan.add(tup)


def test_peptide_plan_iter(peptide_plan_data):
    data, peptide_length = peptide_plan_data
    plan = models.PeptidePlan(peptide_length)

    for tup in data:
        plan.add(tup)

    for combination in plan:
        assert(combination in data)
        data.remove(combination)


def test_peptide_plan_data(peptide_plan_data):
    data, peptide_length = peptide_plan_data
    plan = models.PeptidePlan(peptide_length)

    for tup in data:
        plan.add(tup)

    plan_data = plan.data()

    assert isinstance(plan_data, plan.PeptidePlanData)
    assert isinstance(plan_data.reg_combos, np.ndarray)
    assert isinstance(plan_data.cap_combos, np.ndarray)
    assert np.array_equal(plan_data.reg_combos, np.array([data[0]]))
    assert np.array_equal(plan_data.cap_combos, np.array([data[1]]))
    assert plan_data.length == peptide_length
