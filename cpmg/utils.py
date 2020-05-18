import numpy as np
from rdkit import Chem
from types import GeneratorType

from cpmg.exceptions import MergeError


def split(data, pred):
    """
    https://stackoverflow.com/questions/8793772/how-to-split-a-sequence-according-to-a-predicate
    """

    yes, no = [], []
    for d in data:
        if pred(d):
            yes.append(d)
        else:
            no.append(d)

    return [yes, no]


def get_maximum(data, func):
    try:
        return np.max(list(map(func, data)))
    except ValueError:
        return None


def to_list(data):
    if isinstance(data, (list, tuple)):
        return data

    if isinstance(data, (GeneratorType, map, filter)):
        return list(data)

    return [data]


def connect_mols(*mols, map_nums, stereo=None, clear_map_nums=True):
    """
    Function for combining either one or two molecules at the positions marked by atom map numbers. This function also
    applies the specified stereochemistry at the connected position if applicable, and can clear the map numbers from
    the molecule after making the connection if desired.

    Args:
        map_nums (iterable): An iterable containing two integer elements that specify the map numbers on the atoms
            to connect.
        stereo (str, optional): The stereochemistry to place on the new connection. Can either be 'CW' or 'CCW'.
            Defaults to None.
        clear_map_nums (bool, optional): Whether to clear atom numbers after making the connection or not.
            Defaults to True.

    Raises:
        MergeError: Raised if either no molecules or more than two are provided.

    Returns:
        RDKit Mol: The result of connecting the molecule(s) at the specified positions.
    """

    def update_hydrogen_counts(atom, clear_map_nums):
        """
        Inner method for clearing the atom map number and updating hydrogen counts.
        """

        if clear_map_nums:
            atom.SetAtomMapNum(0)

        if atom.GetNumExplicitHs() != 0:
            atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)

        return atom

    if len(mols) < 1 or len(mols) > 2:
        raise MergeError('Can only merge 1 or 2 molecules at a time.')

    if len(map_nums) != 2:
        raise MergeError('Can only specify 2 distinct map numbers at a time.')

    # find atoms that will form a bond together and update hydrogen counts
    combo, *mols = mols
    for mol in mols:
        combo = Chem.CombineMols(combo, mol)

    # find atoms that will form a bond together and update hydrogen counts
    combo = Chem.RWMol(combo)
    Chem.SanitizeMol(combo)
    try:
        atom1, atom2 = [update_hydrogen_counts(atom, clear_map_nums)
                        for atom in combo.GetAtoms() if atom.GetAtomMapNum() in map_nums]
    except ValueError:
        raise exceptions.MergeError('Could not find 2 atoms with the given map numbers. Check for duplicate map numbers'
                                    ' or that the map numbers are present on the molecules.')

    # create bond, remove hydrogens, and sanitize
    combo.AddBond(atom1.GetIdx(), atom2.GetIdx(), order=Chem.BondType.SINGLE)
    Chem.RemoveHs(combo)
    Chem.SanitizeMol(combo)

    # add stereochemistry as specified
    stereo_center = atom1 if atom1.GetHybridization() == Chem.HybridizationType.SP3 and atom1.GetTotalNumHs() != 2 else atom2
    if stereo == 'CCW':
        stereo_center.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    elif stereo == 'CW':
        stereo_center.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

    return Chem.MolFromSmiles(Chem.MolToSmiles(combo))


def has_atom_map_nums(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            return True

    return False


def get_atom_map_nums(mol):
    atom_map_nums = []
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num != 0:
            atom_map_nums.append((atom.GetIdx(), map_num))

    return atom_map_nums


def clear_atom_map_nums(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def clear_isotopes(mol):
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)


def get_atom_with_map_num(mol, map_num):
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == map_num:
            return atom

    raise RuntimeError(f'Atom map number {map_num} not present on molecule {Chem.MolToSmiles(mol)}')


def remove_atom(mol, atom_idx):
    mol = Chem.RWMol(mol)
    mol.RemoveAtom(atom_idx)
    return mol


def get_classmembers(module):
    import inspect
    import sys
    return inspect.getmembers(sys.modules[module], inspect.isclass)


def get_module_strings(module):
    def get_module_strings_closure():
        return [member.STRING for _, member in get_classmembers(module) if hasattr(member, 'STRING')]

    return get_module_strings_closure


def get_filtered_classmembers(module, pred):
    return [member for _, member in get_classmembers(module) if pred(member)]


def create_factory_function_closure(module, obj_type):
    def factory_function_closure(string, *args):
        for _, member in get_classmembers(module):
            try:
                if string == member.STRING:
                    return member(*args)
            except AttributeError:
                pass

        raise ValueError(f'Unrecognized {obj_type} string: {string}')

    return factory_function_closure
