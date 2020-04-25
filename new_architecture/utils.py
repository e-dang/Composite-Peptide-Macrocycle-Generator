from rdkit import Chem


def has_atom_map_nums(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            return True

    return False


def clear_atom_map_nums(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
