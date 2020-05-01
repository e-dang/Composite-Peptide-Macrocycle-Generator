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


@pytest.mark.parametrize('model_mol', [(TEMPLATE_1), (TEMPLATE_2)])
def test_template_pictet_spangler_fail(model_mol):
    reaction = rxns.TemplatePictetSpangler()

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(model_mol)


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


@pytest.mark.parametrize('model_mol,expected', [(TEMPLATE_2, ALD_RESULT_SMARTS_1)])
def test_aldehyde_cyclization(model_mol, expected):
    reaction = rxns.AldehydeCyclization()

    rxn_smarts = reaction.generate(model_mol)

    assert(rxn_smarts == expected)


@pytest.mark.parametrize('model_mol', [(TEMPLATE_1), (TEMPLATE_3)])
def test_aldehyde_cyclization_fail(model_mol):
    reaction = rxns.AldehydeCyclization()

    with pytest.raises(InvalidMolecule):
        _ = reaction.generate(model_mol)
