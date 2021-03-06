import importlib
import uuid
from copy import deepcopy

import mock
import pytest

import cpmg.decorators
import cpmg.filters
import cpmg.models as models
import cpmg.reactions as rxns
from data.mols import *

mock.patch('cpmg.filters.tpsa_filter', lambda x: x).start()
mock.patch('cpmg.filters.rotatable_bond_filter', lambda x: x).start()
mock.patch('cpmg.filters.molecular_weight_filter', lambda x: x).start()
mock.patch('cpmg.filters.aldehyde_filter', lambda x: x).start()
mock.patch('cpmg.decorators.apply_stereochemistry', lambda x: x).start()
mock.patch('cpmg.decorators.methylate', lambda x: x).start()
mock.patch('cpmg.decorators.carboxyl_to_amide', lambda x: x).start()

import cpmg.generators as generators  # nopep8
importlib.reload(generators)

SIDECHAIN_1 = models.Sidechain.from_dict(TEST_SIDECHAIN_1, _id=str(uuid.uuid4()))
TEMPLATE_1 = models.Template.from_dict(TEST_TEMPLATE_1, _id='temp1')
TEMPLATE_2 = models.Template.from_dict(TEST_TEMPLATE_2, _id='temp2')
TEMPLATE_3 = models.Template.from_dict(TEST_TEMPLATE_3, _id='temp3')
CONNECTION = models.Connection.from_dict(TEST_CONNECTION_2, _id='ethyl')
BACKBONES = [models.Backbone.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_BACKBONE_1, 'alpha'), (TEST_BACKBONE_2, 'beta2'), (TEST_BACKBONE_3, 'beta3')]]
MONOMERS_3 = [models.Monomer.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_MONOMER_3, '98asfh'), (TEST_MONOMER_4, 'adjha82'), (TEST_MONOMER_5, 'admaiof7')]]
MONOMERS_4 = [models.Monomer.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_MONOMER_1, 'ad98fh'), (TEST_MONOMER_2, 'sdwd89cvh'), (TEST_MONOMER_3, '98asfh'), (TEST_MONOMER_2, 'sdwd89cvh')]]
MONOMERS_5 = deepcopy(MONOMERS_4)
MONOMERS_6 = deepcopy(MONOMERS_4)
MONOMERS_7 = deepcopy(MONOMERS_4)
MONOMERS_5.append(models.Monomer.from_dict(TEST_MONOMER_3, _id='98asfh'))
MONOMERS_6.append(models.Monomer.from_dict(TEST_MONOMER_4, _id='af9f7w'))
MONOMERS_7.append(models.Monomer.from_dict(TEST_MONOMER_6, _id='12fakwd'))
PEPTIDE_3 = models.Peptide.from_dict(TEST_PEPTIDE_2, _id='qadfioj2')
PEPTIDE_4_CAP = models.Peptide.from_dict(TEST_PEPTIDE_4_CAP, _id='acd98efh')
PEPTIDE_4 = models.Peptide.from_dict(TEST_PEPTIDE_3, _id='124fiawd')
PEPTIDE_5 = models.Peptide.from_dict(TEST_PEPTIDE_4, _id='af882hd')
PEPTIDE_6 = models.Peptide.from_dict(TEST_PEPTIDE_5, _id='fasef00')
TEMPLATE_PEPTIDE_1 = models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_2, _id='afsfg-30i')
TEMPLATE_PEPTIDE_2 = models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_3, _id='fspva8d')
TEMPLATE_PEPTIDE_3 = models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_4, _id='afvm,svu9')
TEMPLATE_PEPTIDE_4 = models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_5, _id='faopsvp98')
TEMPLATE_PEPTIDE_5 = models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_6, _id='ca0-fiuqe2')
TEMPLATE_PEPTIDE_6 = models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_7, _id='cas.f-3aw')
TEMPLATE_PEPTIDE_7 = models.TemplatePeptide.from_dict(TEST_TEMPLATE_PEPTIDE_8, _id=str(uuid.uuid4()))
REACTION_1 = models.Reaction.from_mols(rxns.FriedelCrafts.TYPE, FC_RESULT_SMARTS_1[0], TEMPLATE_1, SIDECHAIN_1, 3)
REACTION_2 = models.Reaction.from_mols(rxns.FriedelCrafts.TYPE, FC_RESULT_SMARTS_1[0], TEMPLATE_2, SIDECHAIN_1, 3)
REACTION_3 = models.Reaction.from_mols(rxns.FriedelCrafts.TYPE, FC_RESULT_SMARTS_1[0], TEMPLATE_3, SIDECHAIN_1, 3)
REACTION_4 = models.Reaction.from_mols(rxns.TsujiTrost.TYPE, TT_RESULT_SMARTS_1[0], TEMPLATE_1, SIDECHAIN_1, 5)
REACTION_5 = models.Reaction.from_mols(rxns.TsujiTrost.TYPE, TT_RESULT_SMARTS_1[0], TEMPLATE_2, SIDECHAIN_1, 5)
REACTION_6 = models.Reaction.from_mols(rxns.TsujiTrost.TYPE, TT_RESULT_SMARTS_1[0], TEMPLATE_3, SIDECHAIN_1, 5)
REACTION_7 = models.Reaction.from_mols(rxns.AldehydeCyclization.TYPE, ALD_RESULT_SMARTS_1[0], TEMPLATE_2, None, None)
REACTION_8 = models.Reaction.from_mols(rxns.TemplatePictetSpangler.TYPE,
                                       TPS_RESULTS_SMARTS_1[0], TEMPLATE_3, None, None)
REACTION_9 = models.Reaction.from_dict(TEST_REACTION_2, _id=str(uuid.uuid4()))
REACTION_10 = models.Reaction.from_dict(TEST_REACTION_3, _id=str(uuid.uuid4()))
MACROCYCLE_1 = models.Macrocycle.from_dict(TEST_MACROCYCLE_1, _id=str(uuid.uuid4()))
REGIOSQM_PREDICTION = models.RegioSQMPrediction.from_dict(TEST_REGIOSQM_PREDICTION_1, _id=str(uuid.uuid4()))
PKA_PREDICTION = models.pKaPrediction.from_dict(TEST_PKA_PREDICTION_1, _id=str(uuid.uuid4()))


@pytest.mark.parametrize('sidechain, connections, expected_result', [
    (models.Sidechain.from_dict(TEST_SIDECHAIN_1, _id='1afdw'), [CONNECTION], models.Sidechain.from_mol(
        Chem.MolFromSmiles('CCc1ccc(O)cc1'), CONNECTION, TEST_SIDECHAIN_1['shared_id']))
])
def test_sidechain_modifier(sidechain, connections, expected_result):
    generator = generators.SidechainModifier()

    sidechains = generator.generate(sidechain, connections)

    assert len(sidechains) == 1
    assert sidechains[0] == expected_result


@pytest.mark.parametrize('sidechain, backbones, num_results, expected_results', [
    (models.Sidechain.from_dict(TEST_SIDECHAIN_4, _id='af'), BACKBONES, 3, MONOMERS_3)
])
def test_monomer_generator(sidechain, backbones, num_results, expected_results):
    generator = generators.MonomerGenerator()

    monomers = generator.generate(sidechain, backbones)

    assert len(monomers) == num_results
    for monomer, expected_result in zip(monomers, expected_results):
        assert monomer == expected_result


@pytest.mark.parametrize('monomers,peptide_length,expected_result', [
    (MONOMERS_3, 3, PEPTIDE_3),
    (MONOMERS_4, 4, PEPTIDE_4),
    (MONOMERS_5, 5, PEPTIDE_5),
    (MONOMERS_5, 4, PEPTIDE_4_CAP)
])
def test_peptide_generator(monomers, peptide_length, expected_result):
    generator = generators.PeptideGenerator()

    peptides = generator.generate(monomers, peptide_length)

    assert len(peptides) == 1
    assert peptides[0] == expected_result


@pytest.mark.parametrize('monomers, peptide_length', [
    (MONOMERS_6, 4)
])
def test_peptide_generator_fail(monomers, peptide_length):
    generator = generators.PeptideGenerator()

    with pytest.raises(RuntimeError):
        generator.generate(monomers, peptide_length)


@pytest.mark.parametrize('peptide, templates, num_results, expected_results', [
    (PEPTIDE_6, [TEMPLATE_1, TEMPLATE_2, TEMPLATE_3], 6, [TEMPLATE_PEPTIDE_1, TEMPLATE_PEPTIDE_2, TEMPLATE_PEPTIDE_3,
                                                          TEMPLATE_PEPTIDE_4, TEMPLATE_PEPTIDE_5, TEMPLATE_PEPTIDE_6])
])
def test_template_peptide_generator(peptide, templates, num_results, expected_results):
    generator = generators.TemplatePeptideGenerator()

    template_peptides = generator.generate(peptide, templates)

    assert len(template_peptides) == num_results
    template_peptides.sort(key=lambda x: x.kekule)
    expected_results.sort(key=lambda x: x.kekule)
    for template_peptide, expected_result in zip(template_peptides, expected_results):
        assert template_peptide == expected_result


@pytest.mark.parametrize('template_peptide, reaction_combos, num_results, expected_results', [
    (TEMPLATE_PEPTIDE_7, [[REACTION_9, REACTION_10]], 1, [TEST_MACROCYCLE_1['kekule']])
])
def test_macrocycle_generator(template_peptide, reaction_combos, num_results, expected_results):
    generator = generators.MacrocycleGenerator()

    macrocycles = generator.generate(template_peptide, reaction_combos)

    assert len(macrocycles) == num_results
    for macrocycle, expected_result in zip(macrocycles, expected_results):
        assert macrocycle.kekule == expected_result


@pytest.mark.parametrize('nucleophile, templates, impl, num_results, expected_results', [
    (SIDECHAIN_1, [TEMPLATE_1, TEMPLATE_2, TEMPLATE_3],
     rxns.FriedelCrafts(), 3, [REACTION_1, REACTION_2, REACTION_3]),
    (SIDECHAIN_1, [TEMPLATE_1, TEMPLATE_2, TEMPLATE_3], rxns.TsujiTrost(), 3, [REACTION_4, REACTION_5, REACTION_6])
])
def test_intermolecular_reaction_generator(nucleophile, templates, impl, num_results, expected_results):
    with mock.patch('cpmg.generators.repo.BackboneRepository.load', return_value=[models.Backbone.from_dict(TEST_BACKBONE_1, _id=str(uuid.uuid4()))]), \
            mock.patch('cpmg.generators.filters.proxies.repo.RegioSQMRepository.load', return_value=[REGIOSQM_PREDICTION]), \
            mock.patch('cpmg.generators.filters.proxies.repo.pKaRepository.load', return_value=[PKA_PREDICTION]):
        generator = generators.InterMolecularReactionGenerator(impl)
        reactions = generator.generate(nucleophile, templates)

    assert len(reactions) == num_results
    for reaction, expected_result in zip(reactions, expected_results):
        assert reaction == expected_result


@pytest.mark.parametrize('reacting_mol, impl, num_results, expected_results', [
    (TEMPLATE_2, rxns.AldehydeCyclization(), 1, [REACTION_7]),
    (TEMPLATE_3, rxns.TemplatePictetSpangler(), 1, [REACTION_8])
])
def test_intramolecular_reaction_generator(reacting_mol, impl, num_results, expected_results):
    generator = generators.IntraMolecularReactionGenerator(impl)

    reactions = generator.generate(reacting_mol)

    assert len(reactions) == num_results
    for reaction, expected_result in zip(reactions, expected_results):
        assert reaction == expected_result


def test_get_all_generator_strings():
    generator_strings = set(generators.get_all_generator_strings())

    assert generator_strings == {generators.SidechainModifier.STRING,
                                 generators.MonomerGenerator.STRING,
                                 generators.PeptidePlanGenerator.STRING,
                                 generators.PeptideGenerator.STRING,
                                 generators.TemplatePeptideGenerator.STRING,
                                 generators.MacrocycleGenerator.STRING,
                                 generators.ConformerGenerator.STRING,
                                 generators.InterMolecularReactionGenerator.STRING,
                                 generators.IntraMolecularReactionGenerator.STRING}


@pytest.mark.parametrize('generator', [
    (generators.SidechainModifier),
    (generators.MonomerGenerator),
    (generators.PeptidePlanGenerator),
    (generators.PeptideGenerator),
    (generators.ConformerGenerator),
    (generators.TemplatePeptideGenerator),
    (generators.MacrocycleGenerator)
])
def test_create_generator_from_string(generator):
    produced_generator = generators.create_generator_from_string(generator.STRING)

    assert isinstance(produced_generator, generator)


@pytest.mark.parametrize('generator, args', [
    (generators.InterMolecularReactionGenerator, [rxns.FriedelCrafts()]),
    (generators.InterMolecularReactionGenerator, None),
    (generators.IntraMolecularReactionGenerator, [rxns.TemplatePictetSpangler()]),
    (generators.IntraMolecularReactionGenerator, None)
])
def test_create_generator_from_string_reactions(generator, args):
    if args is None:
        produced_generator = generators.create_generator_from_string(generator.STRING)
    else:
        produced_generator = generators.create_generator_from_string(generator.STRING, *args)

    assert isinstance(produced_generator, generator)
    if args is not None:
        assert produced_generator.impl == args
    else:
        assert len(produced_generator.impl) != 0
        assert all([isinstance(impl, (rxns.InterMolecularReaction, rxns.IntraMolecularReaction))
                    for impl in produced_generator.impl])


def test_create_generator_from_string_fail():
    with pytest.raises(ValueError) as err:
        generators.create_generator_from_string('dne')
        exception = str(err.value)

        assert 'dne' in exception
        assert 'generator' in exception
