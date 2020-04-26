import pytest
import new_architecture.generators as generators
import new_architecture.models as models
from tests.new_architecture.data.mols import *

CONNECTION = models.Connection.from_dict(TEST_CONNECTION_2, _id='ethyl')
@pytest.mark.parametrize('data,expected_result', [((models.Sidechain.from_dict(TEST_SIDECHAIN_1, _id='1afdw'), [CONNECTION]), models.Sidechain.from_mol(Chem.MolFromSmiles('CCc1ccc(O)cc1'), CONNECTION, TEST_SIDECHAIN_1['shared_id']))])
def test_sidechain_modifier(data, expected_result):
    generator = generators.SideChainConnectionModifier()

    result = generator.generate(data)

    assert(len(result) == 1)
    assert(result[0] == expected_result)
