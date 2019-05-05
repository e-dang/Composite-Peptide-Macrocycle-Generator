from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import argparse
from pathlib import Path
from utils import read_multiple_sdf


def alternate_connection_point(mol):

    print(Chem.MolToSmiles(mol))
    attached = set()
    methyl = Chem.MolFromSmiles('C')
    methyl.GetAtoms()[0].SetAtomMapNum(2)
    mols = []
    mod_mol = Chem.RWMol(mol)
    patt = Chem.MolFromSmarts('[CH3]*')  # methyl carbon is attachment point

    for i in range(len(mol.GetAtoms())):
        matches = mod_mol.GetSubstructMatches(patt, useChirality=False)
        for pairs in matches:
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(mod_mol, atom_idx)
                if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 3:  # this atom is methyl carbon

                    # add neighbor idx to list of atoms that have had methyl carbon attached
                    try:
                        attached.add(Chem.Atom.GetNeighbors(atom)[0].GetIdx())
                        mod_mol.RemoveAtom(atom_idx)
                    except:
                        print('Error, that atom has already been added to the set!')
                        print('Molecule:', Chem.MolToSmiles(mol))

        # find next eligible atom to attach methyl to
        found = False
        for atom in mod_mol.GetAtoms():
            if atom.GetIdx() in attached:
                continue
            elif atom.GetSymbol() == 'C' and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(1)
                found = True
                break
            elif atom.GetSymbol() == 'N' and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(1)
                found = True
                break
            elif atom.GetSymbol() == 'O' and atom.GetTotalNumHs() != 0:
                atom.SetAtomMapNum(1)
                found = True
                break

        # no more eligible atoms that haven't already had methyl attached
        if not found:
            break

        # prepare for attachment
        combo = Chem.RWMol(Chem.CombineMols(mod_mol, methyl))

        mol_atom = None
        methyl_atom = None
        for atom in combo.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                mol_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)
            elif atom.GetAtomMapNum() == 2:
                methyl_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)

        # create bond
        combo.AddBond(mol_atom, methyl_atom, order=Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(combo)
        mols.append(Chem.MolToSmiles(combo))
        mod_mol = combo

    return mols


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--in', dest='input', nargs='+',
                        default=['sidechains_cust.sdf'], help='The sdf file containing monomer side chains.')
    parser.add_argument('-o', '--out', dest='out', default='diverse_side_chains.txt', help='The output text file.')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='chemdraw/', help='The input filepath relative to script')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/side_chains/',
                        help='The ouput filepath relative to script')

    args = parser.parse_args()

    # get absolute filepath to input and output files
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.input]
    fp_out = str(base_path / args.fp_out / args.out)

    mols = read_multiple_sdf(fp_in)
    for mol in mols:
        mod_mol = alternate_connection_point(mol)
        print(mod_mol)
        for smiles in mod_mol:
            Draw.MolToImage(Chem.MolFromSmiles(smiles)).show()
        exit()


if __name__ == '__main__':
    main()
