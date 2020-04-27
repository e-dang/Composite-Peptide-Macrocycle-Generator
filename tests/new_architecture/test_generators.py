import pytest
import new_architecture.generators as generators
import new_architecture.models as models
from tests.new_architecture.data.mols import *
from copy import deepcopy

CONNECTION = models.Connection.from_dict(TEST_CONNECTION_2, _id='ethyl')
BACKBONES = [models.Backbone.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_BACKBONE_1, 'alpha'), (TEST_BACKBONE_2, 'beta2'), (TEST_BACKBONE_3, 'beta3')]]
MONOMERS_3 = [models.Monomer.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_MONOMER_3, '98asfh'), (TEST_MONOMER_4, 'adjha82'), (TEST_MONOMER_5, 'admaiof7')]]
MONOMERS_4 = [models.Monomer.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_MONOMER_1, 'ad98fh'), (TEST_MONOMER_2, 'sdwd89cvh'), (TEST_MONOMER_3, '98asfh'), (TEST_MONOMER_2, 'sdwd89cvh')]]
MONOMERS_5 = deepcopy(MONOMERS_4)
MONOMERS_6 = deepcopy(MONOMERS_4)
MONOMERS_5.append(models.Monomer.from_dict(TEST_MONOMER_3, _id='98asfh'))
MONOMERS_6.append(models.Monomer.from_dict(TEST_MONOMER_4, _id='af9f7w'))
PEPTIDE_3 = models.Peptide.from_dict(TEST_PEPTIDE_2, _id='qadfioj2')
PEPTIDE_4_CAP = models.Peptide.from_dict(TEST_PEPTIDE_4_CAP, _id='acd98efh')
PEPTIDE_4 = models.Peptide.from_dict(TEST_PEPTIDE_3, _id='124fiawd')
PEPTIDE_5 = models.Peptide.from_dict(TEST_PEPTIDE_4, _id='af882hd')


@pytest.mark.parametrize('data,expected_result', [((models.Sidechain.from_dict(TEST_SIDECHAIN_1, _id='1afdw'), [CONNECTION]), models.Sidechain.from_mol(Chem.MolFromSmiles('CCc1ccc(O)cc1'), CONNECTION, TEST_SIDECHAIN_1['shared_id']))])
def test_sidechain_modifier(data, expected_result):
    generator = generators.SidechainModifier()

    sidechains = generator.generate(data)

    assert(len(sidechains) == 1)
    assert(sidechains[0] == expected_result)


@pytest.mark.parametrize('data,num_results,expected_results', [((models.Sidechain.from_dict(TEST_SIDECHAIN_4, _id='af'), BACKBONES), 3, MONOMERS_3)])
def test_monomer_generator(data, num_results, expected_results):
    generator = generators.MonomerGenerator()

    monomers = generator.generate(data)

    assert(len(monomers) == num_results)
    for monomer, expected_result in zip(monomers, expected_results):
        assert(monomer == expected_result)


@pytest.mark.parametrize('monomers,peptide_length,expected_result', [(MONOMERS_3, 3, PEPTIDE_3), (MONOMERS_4, 4, PEPTIDE_4), (MONOMERS_5, 5, PEPTIDE_5), (MONOMERS_5, 4, PEPTIDE_4_CAP)])
def test_peptide_generator(monomers, peptide_length, expected_result):
    generator = generators.PeptideGenerator(peptide_length)

    peptides = generator.generate(monomers)

    assert(len(peptides) == 1)
    assert(peptides[0] == expected_result)


@pytest.mark.parametrize('monomers,peptide_length', [(MONOMERS_6, 4)])
def test_peptide_generator_fail(monomers, peptide_length):
    generator = generators.PeptideGenerator(peptide_length)

    with pytest.raises(RuntimeError):
        peptides = generator.generate(monomers)
