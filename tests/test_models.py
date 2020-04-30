import uuid

import pytest
from rdkit import Chem

import cpmg.exceptions as exceptions
import cpmg.models as models
import cpmg.utils as utils
from data.mols import *


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
    kekule = 'CC(C)(C)OC(=O)OC/C=C/C1=CC(CCC(=O)ON2C(=O)CCC2=O)=CC=C1'
    oligomerization_kekule = f'C/C=C/C1=CC(CC[CH:{models.Template.OLIGOMERIZATION_MAP_NUM}]=O)=CC=C1'
    friedel_crafts_kekule = f'[*:{models.Template.FRIEDEL_CRAFTS_WC_MAP_NUM}]/C=C/[CH3:{models.Template.FRIEDEL_CRAFTS_EAS_MAP_NUM}]'
    template = models.Template.from_mol(Chem.MolFromSmiles(kekule), oligomerization_kekule, friedel_crafts_kekule)
    return template, kekule, oligomerization_kekule, friedel_crafts_kekule


@pytest.fixture(params=[TEST_TEMPLATE_1, TEST_TEMPLATE_2, TEST_TEMPLATE_3])
def template_from_dict(request):
    _id = str(uuid.uuid4())
    template = models.Template.from_dict(request.param, _id=_id)
    return template, request.param, _id


@pytest.fixture()
def sidechain_from_mol():
    return models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2'), models.Connection.from_dict(TEST_CONNECTION_1, _id='ethyl'), 'q')


@pytest.fixture()
def sidechain_from_dict():
    return models.Sidechain.from_dict(TEST_SIDECHAIN_3, _id='12498dfhjwu')


@pytest.fixture()
def monomer_from_mol():
    return models.Monomer.from_mol(Chem.MolFromSmiles('NCC(CC1=CC=CC2=N[NH]C(=O)N12)C(=O)O'), models.Backbone.from_dict(TEST_BACKBONE_3, _id='beta3'), models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC=CC2=N[NH]C(=O)N12'), models.Connection.from_dict(TEST_CONNECTION_1, _id='methyl'), 'af'))


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

@pytest.fixture()
def reaction_from_mols():
    rxn_type = 'tsuji_trost'
    smarts = '(c1c([*:3])ccc([OH:4])c1.C(=C/[*:50])\\[CH3:2])>>C(=C/[*:50])\\[CH2:2][O:4]c1ccc([*:3])cc1'
    template = models.Template.from_dict(TEST_TEMPLATE_1, _id=str(uuid.uuid4()))
    reacting_mol = models.Sidechain.from_dict(TEST_SIDECHAIN_1, _id=str(uuid.uuid4()))
    rxn_atom_idx = 5
    reaction = models.Reaction.from_mols(rxn_type, smarts, template, reacting_mol, rxn_atom_idx)
    return reaction, rxn_type, smarts, template, reacting_mol, rxn_atom_idx


@pytest.fixture()
def reaction_from_dict():
    _id = str(uuid.uuid4())
    reaction = models.Reaction.from_dict(TEST_REACTION_1, _id=_id)
    return reaction, TEST_REACTION_1, _id


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


@pytest.mark.parametrize('template', [({'oligomerization_kekule': 'C/C=C/C1=CC(CC[CH:1]=O)=CC=C1',
                                        'friedel_crafts_kekule': '[*:201]/C=C/[CH3:200]'}),
                                      ({'oligomerization_kekule': 'C/C=C/C1=CC=C(F)C(C[C@@H](CC=O)[CH:1]=O)=C1',
                                        'friedel_crafts_kekule': '[*:201]/C=C/[CH3:200]'}),
                                      ({'oligomerization_kekule': 'C#CCCC[C@](C=O)(CC1=CC=CC(/C=C/C)=C1)C[CH:1]=O',
                                        'friedel_crafts_kekule': '[*:201]/C=C/[CH3:200]'})])
def test_template_validate(template):
    assert(models.Template.validate(template))


@pytest.mark.parametrize('template', [({'oligomerization_kekule': 'C/C=C/C1=CC(CCC=O)=CC=C1',
                                        'friedel_crafts_kekule': '[*:201]/C=C/[CH3:200]'}),
                                      ({'oligomerization_kekule': 'C/C=C/C1=CC(CC[CH:1]=O)=CC=C1',
                                          'friedel_crafts_kekule': '[*:201]/C=C/C'}),
                                      ({'oligomerization_kekule': 'C/C=C/C1=CC(CC[CH:1]=O)=CC=C1',
                                          'friedel_crafts_kekule': '*/C=C/[CH3:200]'})])
def test_validate_template_fail(template):
    with pytest.raises(exceptions.InvalidMolecule):
        models.Template.validate(template)


def test_template_from_mol(template_from_mol):
    template, kekule, oligomerization_kekule, friedel_crafts_kekule = template_from_mol
    assert(template._id == None)
    assert(template.binary != None)
    assert(isinstance(Chem.Mol(template.binary), Chem.Mol))
    assert(template.kekule == kekule)
    assert(template.oligomerization_kekule == oligomerization_kekule)
    assert(template.friedel_crafts_kekule == friedel_crafts_kekule)


def test_template_from_dict(template_from_dict):
    template, template_dict, _id = template_from_dict
    assert(template._id == _id)
    assert(template.binary == template_dict['binary'])
    assert(isinstance(Chem.Mol(template.binary), Chem.Mol))
    assert(template.kekule == template_dict['kekule'])
    assert(template.oligomerization_kekule == template_dict['oligomerization_kekule'])


def test_template_to_dict(template_from_dict):
    template, template_dict, _ = template_from_dict
    assert(template.to_dict() == template_dict)


def test_sidechain_from_mol(sidechain_from_mol):
    assert(sidechain_from_mol._id == None)
    assert(sidechain_from_mol.binary != None)
    assert(sidechain_from_mol.kekule == 'CC1=CC(=O)C2=C([NH]1)SC=C2')
    assert(sidechain_from_mol.connection == 'C')
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
    assert(sidechain_from_dict.connection == 'C')
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
    assert(monomer_from_mol.connection == 'C')
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


def test_template_peptide_eq(template_peptide_from_dict):
    assert(template_peptide_from_dict == models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_1, _id='doesnt matter'))
    assert(template_peptide_from_dict != models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_2, _id='doesnt matter'))

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

def test_reaction_from_mols(reaction_from_mols):
    reaction, rxn_type, smarts, template, reacting_mol, rxn_atom_idx = reaction_from_mols
    assert(reaction._id == None)
    assert(reaction.type == rxn_type)
    assert(reaction.smarts == smarts)
    assert(reaction.binary != None)
    assert(isinstance(AllChem.ChemicalReaction(reaction.binary), AllChem.ChemicalReaction))
    assert(reaction.template == template._id)
    assert(reaction.reacting_mol == reacting_mol.shared_id)
    assert(reaction.rxn_atom_idx == rxn_atom_idx)


def test_reaction_from_dict(reaction_from_dict):
    reaction, reaction_dict, _id = reaction_from_dict
    assert(reaction._id == _id)
    assert(reaction.type == reaction_dict['type'])
    assert(reaction.smarts == reaction_dict['smarts'])
    assert(reaction.binary != None)
    assert(isinstance(AllChem.ChemicalReaction(reaction.binary), AllChem.ChemicalReaction))
    assert(reaction.template == reaction_dict['template'])
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
