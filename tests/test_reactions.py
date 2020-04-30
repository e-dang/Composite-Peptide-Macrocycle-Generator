import pytest
import uuid
import cpmg.reactions as rxns
import cpmg.models as models
from cpmg.exceptions import InvalidMolecule
from data.mols import *

SIDECHAIN_1 = models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC=C(O)C=C1'),
                                        models.Connection.from_dict(TEST_CONNECTION_1), str(uuid.uuid4()))
SIDECHAIN_2 = models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=CC=C[NH]1'),
                                        models.Connection.from_dict(TEST_CONNECTION_1), str(uuid.uuid4()))
SIDECHAIN_3 = models.Sidechain.from_mol(Chem.MolFromSmiles('CC1=C[NH]C2=CC=CC=C12'),
                                        models.Connection.from_dict(TEST_CONNECTION_1), str(uuid.uuid4()))
MONOMER = models.Monomer.from_dict(TEST_MONOMER_1, _id=str(uuid.uuid4()))
TEMPLATE_1 = models.Template.from_dict(TEST_TEMPLATE_1, _id=str(uuid.uuid4()))
TEMPLATE_2 = models.Template.from_dict(TEST_TEMPLATE_2, _id=str(uuid.uuid4()))
TEMPLATE_3 = models.Template.from_dict(TEST_TEMPLATE_3, _id=str(uuid.uuid4()))

FC_RESULT_SMARTS_1 = [f'(Oc1[cH:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]cc([*:{rxns.FriedelCrafts.NUCLEOPHILE_WC_MAP_NUM}])cc1.C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>Oc1ccc([*:{rxns.FriedelCrafts.NUCLEOPHILE_WC_MAP_NUM}])c[c:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]1[CH2:{models.Template.EAS_MAP_NUM}]/C=C/[*:{models.Template.WC_MAP_NUM_1}]']
FC_RESULT_SMARTS_2 = [f'(O=[C:{rxns.FriedelCrafts.BACKBONE_CARBOXYL_MAP_NUM}]([C@@H]1C[C@H](Oc2cc[cH:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]cc2)CN1[*:{rxns.FriedelCrafts.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.FriedelCrafts.C_TERM_WILDCARD_MAP_NUM}].C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>O=[C:{rxns.FriedelCrafts.BACKBONE_CARBOXYL_MAP_NUM}]([C@@H]1C[C@H](Oc2cc[c:{rxns.FriedelCrafts.NUCLEOPHILE_EAS_MAP_NUM}]([CH2:{models.Template.EAS_MAP_NUM}]/C=C/[*:{models.Template.WC_MAP_NUM_1}])cc2)CN1[*:{rxns.FriedelCrafts.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.FriedelCrafts.C_TERM_WILDCARD_MAP_NUM}]']
TT_RESULT_SMARTS_1 = [f'(c1c([OH:{rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM}])ccc([*:{rxns.TsujiTrost.NUCLEOPHILE_WC_MAP_NUM}])c1.C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{models.Template.EAS_MAP_NUM}])>>C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH2:{models.Template.EAS_MAP_NUM}][O:{rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM}]c1ccc([*:{rxns.TsujiTrost.NUCLEOPHILE_WC_MAP_NUM}])cc1']
PS_RESULT_SMARTS_1 = [f'(O=[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}]([C@@H](C[*:{models.Template.WC_MAP_NUM_1}])C[CH:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]=O)[NH:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]C(Cc1[nH]cc[cH:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]1)[*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])>>O=[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}]1[C@@H](C[*:{models.Template.WC_MAP_NUM_1}])C[C:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]2[c:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]3cc[nH]c3CC([*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])[N:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]12']
PS_RESULT_SMARTS_2 = [f'(C#CCCC[C@@](C[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}](=O)[NH:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]C(Cc1[nH]cc[cH:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]1)[*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])(C[*:{models.Template.WC_MAP_NUM_1}])[CH:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]=O)>>C#CCCC[C@]1(C[*:{models.Template.WC_MAP_NUM_1}])C[C:{rxns.PictetSpangler.TEMPLATE_OLIGOMERIZATION_MAP_NUM}](=O)[N:{rxns.PictetSpangler.BACKBONE_NITROGEN_MAP_NUM}]2C([*:{rxns.PictetSpangler.C_TERM_WILDCARD_MAP_NUM}])Cc3[nH]cc[c:{rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM}]3[C:{rxns.PictetSpangler.TEMPLATE_ALDEHYDE_C_MAP_NUM}]12']
PYR_RESULT_SMARTS_1 = [f'(c1ccc2c(c1)[nH][cH:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}][c:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]2CC([NH:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}][*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}].C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}])>>C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH2:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}][C@:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]12CC([*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}])[N:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}]([*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[C@@H:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}]1Nc1ccccc12',
                       f'(c1ccc2c(c1)[nH][cH:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}][c:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]2CC([NH:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}][*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}].C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH3:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}])>>C(=C/[*:{models.Template.WC_MAP_NUM_1}])\\[CH2:{rxns.Pyrroloindoline.TEMPLATE_EAS_MAP_NUM}][C@@:{rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM}]12CC([*:{rxns.Pyrroloindoline.C_TERM_WILDCARD_MAP_NUM}])[N:{rxns.Pyrroloindoline.BACKBONE_NITROGEN_MAP_NUM}]([*:{rxns.Pyrroloindoline.N_TERM_WILDCARD_MAP_NUM}])[C@H:{rxns.Pyrroloindoline.ADJ_CARBON_MAP_NUM}]1Nc1ccccc12']
TPS_RESULTS_SMARTS_1 = [f'(C#CCCC[C@@](CC(=O)[NH:{rxns.TemplatePictetSpangler.NITROGEN_MAP_NUM}][*:{models.Template.WC_MAP_NUM_1}])(Cc1cc([*:{models.Template.WC_MAP_NUM_2}])cc[cH:{rxns.TemplatePictetSpangler.EAS_MAP_NUM}]1)[CH:{rxns.TemplatePictetSpangler.ALDEHYDE_C_MAP_NUM}]=O)>>C#CCCC[C@]12CC(=O)[N:{rxns.TemplatePictetSpangler.NITROGEN_MAP_NUM}]([*:{models.Template.WC_MAP_NUM_1}])[C@H:{rxns.TemplatePictetSpangler.ALDEHYDE_C_MAP_NUM}]1[c:{rxns.TemplatePictetSpangler.EAS_MAP_NUM}]1ccc([*:{models.Template.WC_MAP_NUM_2}])cc1C2']


@pytest.mark.parametrize('model_mol,template,atom_idx,expected', [(SIDECHAIN_1, TEMPLATE_1, 3, FC_RESULT_SMARTS_1),
                                                                  (SIDECHAIN_1, TEMPLATE_2, 3, FC_RESULT_SMARTS_1),
                                                                  (SIDECHAIN_1, TEMPLATE_3, 3, FC_RESULT_SMARTS_1),
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


@pytest.mark.parametrize('model_mol,template,atom_idx', [(SIDECHAIN_1, TEMPLATE_1, 5),
                                                         (SIDECHAIN_1, TEMPLATE_2, 5),
                                                         (SIDECHAIN_1, TEMPLATE_3, 5),
                                                         (MONOMER, TEMPLATE_1, 13),
                                                         (MONOMER, TEMPLATE_2, 13),
                                                         (MONOMER, TEMPLATE_3, 13)])
def test_friedel_crafts_fail(model_mol, template, atom_idx):
    reaction = rxns.FriedelCrafts()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(reacting_mol, template, reacting_atom, model_mol)


@pytest.mark.parametrize('model_mol,template,atom_idx,expected', [(SIDECHAIN_1, TEMPLATE_1, 5, TT_RESULT_SMARTS_1),
                                                                  (SIDECHAIN_1, TEMPLATE_2, 5, TT_RESULT_SMARTS_1),
                                                                  (SIDECHAIN_1, TEMPLATE_3, 5, TT_RESULT_SMARTS_1)])
def test_tsuji_trost(model_mol, template, atom_idx, expected):
    reaction = rxns.TsujiTrost()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM)

    rxn_smarts = reaction.generate(reacting_mol, template, reacting_atom, model_mol)

    assert(rxn_smarts == expected)


@pytest.mark.parametrize('model_mol,template,atom_idx', [(SIDECHAIN_1, TEMPLATE_1, 4),
                                                         (SIDECHAIN_1, TEMPLATE_2, 4),
                                                         (SIDECHAIN_1, TEMPLATE_3, 4)])
def test_tsuji_trost_fail(model_mol, template, atom_idx):
    reaction = rxns.TsujiTrost()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.TsujiTrost.NUCLEOPHILE_EAS_MAP_NUM)

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(reacting_mol, template, reacting_atom, model_mol)


@pytest.mark.parametrize('model_mol,template,atom_idx,expected', [(SIDECHAIN_2, TEMPLATE_2, 2, PS_RESULT_SMARTS_1),
                                                                  (SIDECHAIN_2, TEMPLATE_3, 2, PS_RESULT_SMARTS_2)])
def test_pictet_spangler(model_mol, template, atom_idx, expected):
    reaction = rxns.PictetSpangler()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM)

    rxn_smarts = reaction.generate(reacting_mol, template, reacting_atom, model_mol)

    assert(rxn_smarts == expected)


@pytest.mark.parametrize('model_mol,template,atom_idx', [(SIDECHAIN_2, TEMPLATE_1, 2), (SIDECHAIN_2, TEMPLATE_2, 3)])
def test_pictet_spangler_fail(model_mol, template, atom_idx):
    reaction = rxns.PictetSpangler()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.PictetSpangler.NUCLEOPHILE_EAS_MAP_NUM)

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(reacting_mol, template, reacting_atom, model_mol)


@pytest.mark.parametrize('model_mol,expected', [(TEMPLATE_3, TPS_RESULTS_SMARTS_1)])
def test_template_pictet_spangler(model_mol, expected):
    reaction = rxns.TemplatePictetSpangler()

    rxn_smarts = reaction.generate(model_mol)

    assert(rxn_smarts == expected)


@pytest.mark.parametrize('model_mol,template,atom_idx,expected', [(SIDECHAIN_3, TEMPLATE_1, 1, PYR_RESULT_SMARTS_1),
                                                                  (SIDECHAIN_3, TEMPLATE_2, 1, PYR_RESULT_SMARTS_1),
                                                                  (SIDECHAIN_3, TEMPLATE_3, 1, PYR_RESULT_SMARTS_1)])
def test_pyrroloindoline(model_mol, template, atom_idx, expected):
    reaction = rxns.Pyrroloindoline()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM)

    rxn_smarts = reaction.generate(reacting_mol, template, reacting_atom, model_mol)

    assert(rxn_smarts == expected)


@pytest.mark.parametrize('model_mol,template,atom_idx', [(SIDECHAIN_3, TEMPLATE_1, 2),
                                                         (SIDECHAIN_3, TEMPLATE_2, 2),
                                                         (SIDECHAIN_3, TEMPLATE_3, 2)])
def test_pyrroloindoline_fail(model_mol, template, atom_idx):
    reaction = rxns.Pyrroloindoline()
    reacting_mol = model_mol.mol
    reacting_atom = reacting_mol.GetAtomWithIdx(atom_idx)
    reacting_atom.SetAtomMapNum(rxns.Pyrroloindoline.NUCLEOPHILE_EAS_MAP_NUM)

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(reacting_mol, template, reacting_atom, model_mol)
