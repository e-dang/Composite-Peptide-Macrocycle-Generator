import numpy as np
from rdkit import Chem


def get_maximum(data, func):
    try:
        return np.max(list(map(func, data)))
    except ValueError:
        return None


def to_list(data):
    if isinstance(data, dict):
        return [data]

    return data


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
