# import os

# import pytest
# from rdkit import Chem

# import macrocycles.exceptions as exceptions
# import macrocycles.utils as utils
# import macrocycles.config as config


# @pytest.mark.parametrize('mols,map_nums,stereo,clear_map_nums,expected_smiles', [
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'),
#      (1, 2), None, True, 'O=C(O)C(Cc1ccc(O)cc1)[N]C(=O)[C@@H]1C[C@H](Oc2ccccc2)CN1'),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'), (1, 2),
#      None, False, 'O=C(O)C(Cc1ccc(O)cc1)[N:1][C:2](=O)[C@@H]1C[C@H](Oc2ccccc2)CN1'),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'),
#      (1, 2), 'CW', True, 'O=C(O)C(Cc1ccc(O)cc1)[N]C(=O)[C@@H]1C[C@H](Oc2ccccc2)CN1'),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'), (1, 2),
#      'CW', False, 'O=C(O)C(Cc1ccc(O)cc1)[N:1][C:2](=O)[C@@H]1C[C@H](Oc2ccccc2)CN1'),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'),
#      (1, 2), 'CCW', True, 'O=C(O)C(Cc1ccc(O)cc1)[N]C(=O)[C@@H]1C[C@H](Oc2ccccc2)CN1'),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'), (1, 2),
#      'CCW', False, 'O=C(O)C(Cc1ccc(O)cc1)[N:1][C:2](=O)[C@@H]1C[C@H](Oc2ccccc2)CN1'),
#     (('[CH3:2]C1=CC=C(O)C=C1', 'N[CH2:1]C(=O)[OH]'), (1, 2), None, True, 'NC(Cc1ccc(O)cc1)C(=O)O'),
#     (('[CH3:2]C1=CC=C(O)C=C1', 'N[CH2:1]C(=O)[OH]'), (1, 2), None, False, 'N[CH:1](C(=O)O)[CH2:2]c1ccc(O)cc1'),
#     (('[CH3:2]C1=CC=C(O)C=C1', 'N[CH2:1]C(=O)[OH]'), (1, 2), 'CW', True, 'N[C@H](Cc1ccc(O)cc1)C(=O)O'),
#     (('[CH3:2]C1=CC=C(O)C=C1', 'N[CH2:1]C(=O)[OH]'), (1, 2), 'CW', False, 'N[C@@H:1](C(=O)O)[CH2:2]c1ccc(O)cc1'),
#     (('[CH3:2]C1=CC=C(O)C=C1', 'N[CH2:1]C(=O)[OH]'), (1, 2), 'CCW', True, 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'),
#     (('[CH3:2]C1=CC=C(O)C=C1', 'N[CH2:1]C(=O)[OH]'), (1, 2), 'CCW', False, 'N[C@H:1](C(=O)O)[CH2:2]c1ccc(O)cc1'),
#     (('[CH3:1]/C=C/C1=CC([CH2:2]CC=O)=CC=C1', ), (1, 2), None, True, 'O=CCC1CC=Cc2cccc1c2'),
#     (('[CH3:1]/C=C/C1=CC([CH2:2]CC=O)=CC=C1', ), (1, 2), None, False, 'O=CC[CH:2]1c2cccc(c2)C=C[CH2:1]1'),
#     (('[CH3:1]/C=C/C1=CC([CH2:2]CC=O)=CC=C1', ), (1, 2), 'CW', True, 'O=CC[C@H]1CC=Cc2cccc1c2'),
#     (('[CH3:1]/C=C/C1=CC([CH2:2]CC=O)=CC=C1', ), (1, 2), 'CW', False, 'O=CC[C@@H:2]1c2cccc(c2)C=C[CH2:1]1'),
#     (('[CH3:1]/C=C/C1=CC([CH2:2]CC=O)=CC=C1', ), (1, 2), 'CCW', True, 'O=CC[C@@H]1CC=Cc2cccc1c2'),
#     (('[CH3:1]/C=C/C1=CC([CH2:2]CC=O)=CC=C1', ), (1, 2), 'CCW', False, 'O=CC[C@H:2]1c2cccc(c2)C=C[CH2:1]1')
# ], ids=['monomer proline clear',
#         'monomer proline',
#         'monomer proline CW clear',
#         'monomer proline CW',
#         'monomer proline CCW clear',
#         'monomer proline CCW',
#         'sidechain connection clear',
#         'sidechain connection',
#         'sidechain connection CW clear',
#         'sidechain connection CW',
#         'sidechain connection CCW clear',
#         'sidechain connection CCW',
#         'intramolecular clear',
#         'intramolecular',
#         'intramolecular CW clear',
#         'intramolecular CW',
#         'intramolecular CCW clear',
#         'intramolecular CCW'])
# def test_connect_mols(mols, map_nums, stereo, clear_map_nums, expected_smiles):
#     """
#     Test for utils.connect_mols().

#     Args:
#         mols (iterable[str]): The SMILES strings of the molecules to connect
#         map_nums (iterable[int]): The atom map numbers of the atoms to connect.
#         stereo (str): The stereochemistry to apply. Can be either 'CW' or 'CCW'.
#         clear_map_nums (bool): Whether to clear the map numbers or not.
#         expected_smiles (str): The expected SMILES string from merging the molecules.
#     """

#     mols = list(map(Chem.MolFromSmiles, mols))

#     new_mol = utils.connect_mols(*mols, map_nums=map_nums, stereo=stereo, clear_map_nums=clear_map_nums)

#     print(Chem.MolToSmiles(new_mol))
#     # assert Chem.MolToSmiles(new_mol) == expected_smiles


# @pytest.mark.parametrize('mols,map_nums', [
#     ([], (1, 2)),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1',
#       'O=[CH:3][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'), (1, 2))
# ], ids=['zero mols', 'three mols'])
# def test_connect_mols_fail_wrong_number_mols(mols, map_nums):
#     """
#     Test of utils.connect_mols() that fails due to incorrect number of molecules given.

#     Args:
#         mols (iterable[str]): The SMILES strings of the molecules.
#         map_nums (iterable[int]): The map numbers on the atoms to connect.
#     """

#     mols = list(map(Chem.MolFromSmiles, mols))

#     with pytest.raises(exceptions.MergeError) as error:
#         _ = utils.connect_mols(*mols, map_nums=map_nums)
#         assert str(error.value) == 'Can only merge 1 or 2 molecules at a time.'


# @pytest.mark.parametrize('mols,map_nums', [
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)C[NH:2]1'), (1, 2)),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:3][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'), (1, 2))
# ], ids=['two identical map nums', 'missing map number'])
# def test_connect_mols_fail_wrong_map_numbers(mols, map_nums):
#     """
#     Test of utils.connect_mols() that fails due to duplicate or missing map numbers.

#     Args:
#         mols (iterable[str]): The SMILES strings of the molecules.
#         map_nums (iterable[int]): The map numbers on the atoms to connect.
#     """

#     mols = list(map(Chem.MolFromSmiles, mols))

#     with pytest.raises(exceptions.MergeError) as error:
#         _ = utils.connect_mols(*mols, map_nums=map_nums)
#         assert str(error.value) == 'Could not find 2 atoms with the given map numbers. Check for duplicate map numbers \
#                                or that the map numbers are present on the molecules.'


# @pytest.mark.parametrize('mols,map_nums', [
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)CN1'), (1, 2, 3)),
#     (('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 'O=[CH:2][C@@H]1C[C@H](OC2=CC=CC=C2)C[NH:2]1'), (1,)),
# ], ids=['3 map nums', 'one map num'])
# def test_connect_mols_fail_wrong_number_of_map_numbers(mols, map_nums):
#     """
#     Test of utils.connect_mols() that fails due to incorrect number of map numbers given.

#     Args:
#         mols (iterable[str]): The SMILES strings of the molecules.
#         map_nums (iterable[int]): The map numbers on the atoms to connect.
#     """

#     mols = list(map(Chem.MolFromSmiles, mols))

#     with pytest.raises(exceptions.MergeError) as error:
#         _ = utils.connect_mols(*mols, map_nums=map_nums)
#         assert str(error.value) == 'Can only specify 2 distinct map numbers at a time.'


# @pytest.mark.parametrize('smiles,atom_map_num,expected_element,expected_hydrogens', [
#     ('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 1, 'N', 2),
#     ('[NH2:1]C(CC1=CC=C(O)C=C1)[C:2](=O)O', 2, 'C', 0)
# ], ids=['single map num', 'multiple map nums'])
# def test_find_atom(smiles, atom_map_num, expected_element, expected_hydrogens):
#     """
#     Test for utils.find_atom() when given SMILES strings with the requested atom map number.

#     Args:
#         smiles (str): The smiles string of a molecule.
#         atom_map_num (int): The atom map number to find on the molecule.
#         expected_element (str): The symbol of the atom tagged with the specified atom map number.
#         expected_hydrogens (int): The number of hydrogens on the atom tagged with the specified atom map number.
#     """

#     mol = Chem.MolFromSmiles(smiles)

#     atom = utils.find_atom(mol, atom_map_num)

#     assert atom.GetAtomMapNum() == atom_map_num
#     assert atom.GetSymbol() == expected_element
#     assert atom.GetTotalNumHs() == expected_hydrogens


# @pytest.mark.parametrize('smiles,atom_map_num', [
#     ('[NH2:1]C(CC1=CC=C(O)C=C1)C(=O)O', 2),
#     ('NC(CC1=CC=C(O)C=C1)C(=O)O', 1)
# ], ids=['wrong map num', 'zero map nums'])
# def test_find_atom_fail(smiles, atom_map_num):
#     """
#     Test for utils.find_atoms() when given SMILES strings either missing the requested atom map number or void of any
#     atom map numbers.

#     Args:
#         smiles (str): The SMILES string.
#         atom_map_num (int): The atom map number of the atom to find.
#     """

#     mol = Chem.MolFromSmiles(smiles)

#     atom = utils.find_atom(mol, atom_map_num)

#     assert atom == None


# @pytest.mark.parametrize('smiles', [('NC(CC1=CC=C(O)C=C1)C(=O)O')], ids=['tyrosine'])
# def test_atom_to_wildcard(smiles):
#     """
#     Test for utils.atom_to_wildcard().

#     Args:
#         smiles (str): The SMILES string of the molecule to test on.
#     """

#     mol = Chem.MolFromSmiles(smiles)
#     for atom in mol.GetAtoms():
#         utils.atom_to_wildcard(atom)

#         assert atom.GetAtomicNum() == 0
#         assert atom.GetIsotope() == 0
#         assert atom.GetFormalCharge() == 0
#         assert not atom.GetIsAromatic()
#         assert atom.GetNumExplicitHs() == 0


# @pytest.mark.parametrize('factors,expected_values', [
#     ([(1, 2), (3, 4)], [[1, 3], [1, 4], [2, 3], [2, 4]])
# ], ids=['two iterables'])
# def test_random_order_cartesian_product_correct_product(factors, expected_values):
#     """
#     Test for utils.random_order_cartesian_product(). Tests that the given cartesian product contains all the expected
#     elements.

#     Args:
#         factors (iterable[iterables]): The iterables used to form the cartesian product.
#         expected_values (iterable[iterable]): The expected cartesian product.
#     """

#     product = utils.random_order_cartesian_product(*factors)
#     for element in product:
#         assert element in expected_values


# @pytest.mark.parametrize('factors,expected_number_values', [
#     ([(1, 2, 3, 4, 5), (1, 2, 3, 4, 5), (1, 2, 3, 4, 5)], 125),
#     ([range(10), range(10)], 100)
# ], ids=['length 125', 'length 100'])
# def test_random_order_cartesian_product_length_only(factors, expected_number_values):
#     """
#     Test for utils.random_order_cartesian_product(), only verifying that the size of the returned product is the
#     expected size of the cartesian product.

#     Args:
#         factors (iterable[iterable]): The iterables used to form the cartesian product.
#         expected_number_values (int): The expected size of the cartesian product.
#     """

#     product = list(utils.random_order_cartesian_product(*factors))

#     assert len(product) == expected_number_values


# @pytest.mark.parametrize('filepath,file_nums,expected_value', [
#     ('data.txt', (1,), 'data_1.txt'),
#     ('data.json', (1,), 'data_1.json'),
#     ('test/data.json', (1, 2), 'test/data_1_2.json')
# ], ids=['txt extension', 'json extension', 'multiple file nums'])
# def test_attach_file_num(filepath, file_nums, expected_value):
#     """
#     Test for utils.attach_file_num() that is given properly formatted data.

#     Args:
#         filepath (str): The filepath.
#         file_nums (iterable[ints]): The file numbers to attach to the filepath.
#         expected_value (str): The expected filepath with attached file numbers.
#     """

#     new_filepath = utils.attach_file_num(filepath, *file_nums)

#     assert new_filepath == expected_value


# @pytest.mark.parametrize('filepath,file_nums', [
#     ('data.txt.so', (1,)),
#     ('data', (1,)),
# ], ids=['one extension', 'zero extensions'])
# def test_attach_file_num_fail(filepath, file_nums):
#     """
#     Test for utils.attach_file_num() that is expected to fail and raise an exception due to improper formatted data.

#     Args:
#         filepath (str): The filepath.
#         file_nums (iterable[int]): The file numbers to attach to the filepath.
#     """

#     with pytest.raises(exceptions.InvalidFilePath):
#         _ = utils.attach_file_num(filepath, *file_nums)


# @pytest.fixture()
# def example_filepath():
#     """
#     Creates three files in the tests/data/ directory with attached file numbers and returns the base file path, then
#     cleans up by removing the three files.
#     """

#     paths = [f'{config.PROJECT_DIR}/tests/data/example_file_0.txt',
#              f'{config.PROJECT_DIR}/tests/data/example_file_1.txt',
#              f'{config.PROJECT_DIR}/tests/data/example_file_2.txt']

#     # create files
#     for filepath in paths:
#         with open(filepath, 'w+') as _:
#             pass

#     # return filebase name
#     yield f'{config.PROJECT_DIR}/tests/data/example_file.txt'

#     # delete files
#     for filepath in paths:
#         os.remove(filepath)


# def test_file_rotator(example_filepath):
#     """
#     Test of utils.file_rotator().

#     Args:
#         example_filepath (pytest.Fixure): The example_filepath pytest fixure.
#     """

#     new_filepath = utils.file_rotator(example_filepath)

#     assert new_filepath == f'{config.PROJECT_DIR}/tests/data/example_file_3.txt'


# def test_get_file_num_range(example_filepath):
#     """
#     Test of utils.get_file_num_range().

#     Args:
#         example_filepath (pytest.Fixure): The example_filepath pytest fixure.
#     """

#     low, high = utils.get_file_num_range(example_filepath)

#     assert low == 0
#     assert high == 3
