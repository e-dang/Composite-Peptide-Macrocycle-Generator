import pytest
import new_architecture.validators as validators
import macrocycles.exceptions as exceptions
from rdkit import Chem


@pytest.mark.parametrize('backbone', [(Chem.MolFromSmiles('N[CH2:1]C(=O)O')), (Chem.MolFromSmiles('N[CH2:1]CC(=O)O')), (Chem.MolFromSmiles('NC[CH2:1]C(=O)O'))])
def test_validate_backbone(backbone):
    assert(validators.validate_backbone(backbone))


@pytest.mark.parametrize('backbone', [(Chem.MolFromSmiles('NCC(=O)O')), (Chem.MolFromSmiles('NCCC(=O)O'))])
def test_validate_backbone_fail(backbone):
    with pytest.raises(exceptions.InvalidMolecule):
        validators.validate_backbone(backbone)
