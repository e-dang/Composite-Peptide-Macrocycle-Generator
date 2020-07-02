import csv
import functools
import json
import os
import random
from types import GeneratorType

import numpy as np
from rdkit import Chem

from cpmg.exceptions import MergeError
from cpmg.ranges import Key


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

    if isinstance(data, (GeneratorType, map, filter, Key)):
        return list(data)

    return [data]


def random_order_cartesian_product(*factors):
    """
    Randomly samples the cartesian product of the factors without repeats.
    Code from https://stackoverflow.com/questions/48686767/how-to-sample-from-cartesian-product-without-repetition

    Yields:
        list: A list of containing a random sample from the cartesian produce of the factors.
    """

    amount = functools.reduce(lambda prod, factor: prod * len(factor), factors, 1)
    index_linked_list = [None, None]
    for max_index in reversed(range(amount)):
        index = random.randint(0, max_index)
        index_link = index_linked_list
        while index_link[1] is not None and index_link[1][0] <= index:  # pylint: disable=E1136
            index += 1
            index_link = index_link[1]
        index_link[1] = [index, index_link[1]]
        items = []
        for factor in factors:
            items.append(factor[index % len(factor)])
            index //= len(factor)
        yield items


def random_sample_cartesian_product(*factors, sample_size=0):
    """
    Randomly samples the cartesian product of the factors without repeats. The explicit sample size makes the execution
    faster than the random_order_cartesian_product variant.

    Code from https://stackoverflow.com/questions/48686767/how-to-sample-from-cartesian-product-without-repetition

    Yields:
        list: A list of containing a random sample from the cartesian produce of the factors.
    """

    amount = functools.reduce(lambda prod, factor: prod * len(factor), factors, 1)
    for max_index in random.sample(list(range(amount)), sample_size):
        items = []
        for factor in factors:
            items.append(factor[max_index % len(factor)])
            max_index //= len(factor)
        yield items


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
        raise MergeError('Could not find 2 atoms with the given map numbers. Check for duplicate map numbers'
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


def atom_to_wildcard(atom):
    """
    Changes the atom to a wildcard atom.

    Args:
        atom (RDKit Atom): The atom to change to a wildcard.
    """

    atom.SetAtomicNum(0)
    atom.SetIsotope(0)
    atom.SetFormalCharge(0)
    atom.SetIsAromatic(False)
    atom.SetNumExplicitHs(0)


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


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in Python, i.e. will suppress all print, even
    if the print originates in a compiled C/Fortran sub-function. This will not suppress raised exceptions, since
    exceptions are printed to stderr just before a script exits, and after the context manager has exited (at least,
    I think that is why it lets exceptions through).

    Code taken from: https://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
    '''

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


def load_text(filepath):
    with open(filepath, 'r') as file:
        return file.readlines()


def save_text(data, filepath):
    with open(filepath, 'w') as file:
        for line in data:
            file.write(line)


def load_json(filepath):
    with open(filepath, 'r') as file:
        return json.load(file)


def save_json(data, filepath):
    with open(filepath, 'w') as file:
        json.dump(data, file)


def load_csv(filepath):
    with open(filepath, 'r') as file:
        return list(csv.reader(file, delimiter=','))


def attach_file_num(filepath, *file_nums, delimiter='_'):
    """
    Function that inserts an underscore and the specified file number to the file name given in the filepath.

    Args:
        filepath (str): The filepath containing the file name to augment.
        file_nums (iterable): The file numbers to attach to the end of the file name in the filepath.

    Returns:
        str: The augmented filepath.
    """

    if filepath.count('.') != 1:
        raise OSError('The given filepath needs to have exactly one \'.\' character.')

    new_fp, ext = filepath.split('.')

    for file_num in file_nums:
        if file_num is not None:
            new_fp += delimiter + str(file_num)

    new_fp += '.' + ext

    return new_fp


def rotate_file(filepath):
    """
    Function that takes a filepath and attaches a underscore and a number before its extension to ensure uniqueness of
    the file name given by the filepath.

    Args:
        filepath (str): The path to the file.

    Returns:
        str: The augmented filepath.
    """

    idx = 0
    while True:
        new_fp = attach_file_num(filepath, idx)
        idx += 1
        if not (os.path.exists(new_fp) and os.path.isfile(new_fp)):
            return new_fp
