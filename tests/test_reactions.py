import pytest
import uuid
import cpmg.reactions as rxns
import cpmg.models as models
from cpmg.exceptions import InvalidMolecule
from data.mols import *

SIDECHAIN = models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC=C(O)C=C1'),
                                      models.Connection.from_dict(TEST_CONNECTION_1), str(uuid.uuid4()))
MONOMER = models.Monomer.from_dict(TEST_MONOMER_1, _id=str(uuid.uuid4()))
TEMPLATE_1 = models.Template.from_dict(TEST_TEMPLATE_1, _id=str(uuid.uuid4()))
TEMPLATE_2 = models.Template.from_dict(TEST_TEMPLATE_2, _id=str(uuid.uuid4()))
TEMPLATE_3 = models.Template.from_dict(TEST_TEMPLATE_3, _id=str(uuid.uuid4()))

FC_RESULT_SMARTS_1 = f'(Oc1[cH:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]cc([*:{rxns.FriedelCrafts.NUCLEOPHILE_WC_MAP_NUM}])cc1.C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>Oc1ccc([*:{rxns.FriedelCrafts.NUCLEOPHILE_WC_MAP_NUM}])c[c:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]1[CH2:{models.Template.EAS_MAP_NUM}]/C=C/[*:{models.Template.WC_MAP_NUM_1}]'
FC_RESULT_SMARTS_2 = f'(O=[C:{rxns.FriedelCrafts.BACKBONE_CARBOXYL_MAP_NUM}]([C@@H]1C[C@H](Oc2cc[cH:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]cc2)CN1[*:{rxns.FriedelCrafts.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.FriedelCrafts.C_TERM_WILDCARD_MAP_NUM}].C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>O=[C:{rxns.FriedelCrafts.BACKBONE_CARBOXYL_MAP_NUM}]([C@@H]1C[C@H](Oc2cc[c:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]([CH2:{models.Template.EAS_MAP_NUM}]/C=C/[*:{models.Template.WC_MAP_NUM_1}])cc2)CN1[*:{rxns.FriedelCrafts.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.FriedelCrafts.C_TERM_WILDCARD_MAP_NUM}]'
TT_RESULT_SMARTS_1 = f'(c1c([OH:{rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM}])ccc([*:{rxns.TsujiTrost.NUCLEOPHILE_WC_MAP_NUM}])c1.C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH2:{models.Template.EAS_MAP_NUM}][O:{rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM}]c1ccc([*:{rxns.TsujiTrost.NUCLEOPHILE_WC_MAP_NUM}])cc1'


@pytest.mark.parametrize('model_mol,template,atom_idx,expected', [(SIDECHAIN, TEMPLATE_1, 3, FC_RESULT_SMARTS_1),
                                                                  (SIDECHAIN, TEMPLATE_2, 3, FC_RESULT_SMARTS_1),
                                                                  (SIDECHAIN, TEMPLATE_3, 3, FC_RESULT_SMARTS_1),
                                                                  (MONOMER, TEMPLATE_1, 10, FC_RESULT_SMARTS_2),
                                                                  (MONOMER, TEMPLATE_2, 10, FC_RESULT_SMARTS_2),
                                                                  (MONOMER, TEMPLATE_3, 10, FC_RESULT_SMARTS_2)])
def test_friedel_crafts(model_mol, template, atom_idx, expected):
    reaction = rxns.FriedelCrafts()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM)

    rxn_smarts = reaction.generate(reacting_mol, template, reacting_atom, model_mol)

    assert(rxn_smarts == expected)


@pytest.mark.parametrize('model_mol,template,atom_idx', [(SIDECHAIN, TEMPLATE_1, 5),
                                                         (SIDECHAIN, TEMPLATE_2, 5),
                                                         (SIDECHAIN, TEMPLATE_3, 5),
                                                         (MONOMER, TEMPLATE_1, 13),
                                                         (MONOMER, TEMPLATE_2, 13),
                                                         (MONOMER, TEMPLATE_3, 13)])
def test_friedel_crafts_fail(model_mol, template, atom_idx):
    reaction = rxns.FriedelCrafts()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(reacting_mol, template, reacting_atom, model_mol)


@pytest.mark.parametrize('model_mol,template,atom_idx,expected', [(SIDECHAIN, TEMPLATE_1, 5, TT_RESULT_SMARTS_1),
                                                                  (SIDECHAIN, TEMPLATE_2, 5, TT_RESULT_SMARTS_1),
                                                                  (SIDECHAIN, TEMPLATE_3, 5, TT_RESULT_SMARTS_1)])
def test_tsuji_trost(model_mol, template, atom_idx, expected):
    reaction = rxns.TsujiTrost()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM)

    rxn_smarts = reaction.generate(reacting_mol, template, reacting_atom, model_mol)
    print(rxn_smarts)
    assert(rxn_smarts == expected)


@pytest.mark.parametrize('model_mol,template,atom_idx', [(SIDECHAIN, TEMPLATE_1, 4),
                                                         (SIDECHAIN, TEMPLATE_2, 4),
                                                         (SIDECHAIN, TEMPLATE_3, 4)])
def test_tsuji_trost_fail(model_mol, template, atom_idx):
    reaction = rxns.TsujiTrost()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM)

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(reacting_mol, template, reacting_atom, model_mol)
