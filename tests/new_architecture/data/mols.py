from rdkit import Chem

TEST_SIDECHAIN_1 = {'type': 'sidechain', 'binary': Chem.MolFromSmiles('CC1=CC=C(O)C=C1').ToBinary(
), 'kekule': 'CC1=CC=C(O)C=C1', 'connection': 'methyl', 'shared_id': 'a'}
TEST_SIDECHAIN_2 = {'type': 'sidechain', 'binary': Chem.MolFromSmiles('CC1=COC2=NC(=O)[NH]C=C12').ToBinary(
), 'kekule': 'CC1=COC2=NC(=O)[NH]C=C12', 'connection': 'methyl', 'shared_id': 's'}
TEST_SIDECHAIN_3 = {'binary': Chem.MolFromSmiles('CC1=CC(=O)C2=C([NH]1)SC=C2').ToBinary(
), 'kekule': 'CC1=CC(=O)C2=C([NH]1)SC=C2', 'connection': 'methyl', 'shared_id': 'q'}
