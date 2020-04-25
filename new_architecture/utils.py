import numpy as np


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


def clear_atom_map_nums(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
