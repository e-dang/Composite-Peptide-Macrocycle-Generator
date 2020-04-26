import pytest
import new_architecture.generators as generators
import new_architecture.models as models
from tests.new_architecture.data.mols import *

CONNECTION = models.Connection.from_dict(TEST_CONNECTION_2, _id='ethyl')
@pytest.mark.parametrize('data,expected_result', [((models.Sidechain.from_dict(TEST_SIDECHAIN_1, _id='1afdw'), [CONNECTION]), models.Sidechain.from_mol(Chem.MolFromSmiles('CCc1ccc(O)cc1'), CONNECTION, TEST_SIDECHAIN_1['shared_id']))])
def test_sidechain_modifier(data, expected_result):
    generator = generators.SidechainModifier()

    result = generator.generate(data)

    assert(len(result) == 1)
    assert(result[0] == expected_result)


BACKBONES = [models.Backbone.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_BACKBONE_1, 'alpha'), (TEST_BACKBONE_2, 'beta2'), (TEST_BACKBONE_3, 'beta3')]]
MONOMERS = [models.Monomer.from_dict(doc, _id=_id) for doc, _id in [(
    TEST_MONOMER_3, '98asfh'), (TEST_MONOMER_4, 'adjha82'), (TEST_MONOMER_5, 'admaiof7')]]


@pytest.mark.parametrize('data,num_results,expected_results', [((models.Sidechain.from_dict(TEST_SIDECHAIN_4, _id='af'), BACKBONES), 3, MONOMERS)])
def test_monomer_generator(data, num_results, expected_results):
    generator = generators.MonomerGenerator()

    results = generator.generate(data)

    assert(len(results) == num_results)
    for result, expected_result in zip(results, expected_results):
        assert(result == expected_result)
