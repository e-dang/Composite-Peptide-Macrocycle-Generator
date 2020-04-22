from rdkit import Chem
from rdkit.Chem import AllChem
import os

CYCLOPEP_DIR = os.path.join('.', 'data', 'cyclopeptides')


def load_pdb(filepath):
    return Chem.MolFromPDBFile(os.path.join(CYCLOPEP_DIR, filepath), removeHs=False)


def remove_confs(mol, max_id):
    for conf in mol.GetConformers():
        if conf.GetId() >= max_id:
            mol.RemoveConformer(conf.GetId())


def align_confs(mol):
    for conf in mol.GetConformers():
        AllChem.AlignMol(mol, mol, prbCid=conf.GetId(), refCid=0)


def overlay(mol1, mol2):
    print(len(mol1.GetAtoms()), mol1.GetNumConformers())
    print(len(mol2.GetAtoms()), mol2.GetNumConformers())
    Chem.SanitizeMol(mol1)
    Chem.SanitizeMol(mol2)
    assert (Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2))
    return AllChem.AlignMol(mol1, mol2, prbCid=0, refCid=0, maxIters=1000)


def get_ring_atoms(mol1):
    ring_atoms = [ring for ring in mol1.GetRingInfo().AtomRings() if
                  len(ring) >= 10]

    return list(set().union(*ring_atoms))


def overlay_ring_atoms(mol1, mol2, mol2_conf_id):
    mol1_atoms = get_ring_atoms(mol1)
    mol2_atoms = get_ring_atoms(mol2)
    return AllChem.AlignMol(mol1, mol2, prbCid=0, refCid=mol2_conf_id, atomMap=list(zip(mol1_atoms, mol2_atoms)), maxIters=1000)


def renumber_atoms(mol1, mol2):
    match = mol1.GetSubstructMatch(mol2)
    mol1 = AllChem.RenumberAtoms(mol1, match)
    for atom1, atom2 in zip(mol1.GetAtoms(), mol2.GetAtoms()):
        assert(atom1.GetSymbol() == atom2.GetSymbol())
    Chem.SanitizeMol(mol1)
    Chem.SanitizeMol(mol2)
    return mol1, mol2


if __name__ == "__main__":
    semi_emp_fp = os.path.join('semi_empirical', 'semi_emp_1_10_lowest.pdb')
    conf_bpp_fp = os.path.join('confbuster++', 'cyclopeptide_1.pdb')
    semi_emp = load_pdb(semi_emp_fp)
    conf_bpp = load_pdb(conf_bpp_fp)
    remove_confs(semi_emp, 10)
    remove_confs(conf_bpp, 1)
    semi_emp, conf_bpp = renumber_atoms(semi_emp, conf_bpp)
    align_confs(semi_emp)
    align_confs(conf_bpp)
    # print(overlay(conf_bpp, semi_emp))
    print(overlay_ring_atoms(conf_bpp, semi_emp, 9))
    # Chem.MolToPDBFile(semi_emp, 'semi_emp_1.pdb')
    # Chem.MolToPDBFile(conf_bpp, 'confb_1.pdb')
