import argparse
import random as rand
from itertools import product
from pathlib import Path
from time import time

from rdkit import Chem
from rdkit.Chem import Draw

# num_mols = 10

# get custom amino acids
# cust_aa = []
# with open('/Users/ericdang/Documents/UCLA_Research/smiles/custom_aa/cust_S_aa.txt', 'r') as f:
#     lines = f.readlines()
#     rand_ind = sorted([rand.randint(0, len(lines)) for x in range(20)])
#     for ind , line in enumerate(lines):
#         if ind == rand_ind[0]:
#             cust_aa.append(Chem.MolFromSmiles(line))
#             rand_ind.pop(0)
#
# # get natural amino acids
# codes = 'A R N D C E Q G H I L K M F P S T W Y V'.split(' ')
# nat_aa = [Chem.rdmolfiles.MolFromFASTA(code) for code in codes]
#
# # create connections
# with open('/Users/ericdang/Documents/UCLA_Research/smiles/peptides/demo_10.txt', 'w') as f:
#     for x in range(num_mols):
#         rand.seed(time())
#         peptide = None
#         for i in range(1, 6):
#             if i % 2 == 0:
#                 rand_ind = rand.randint(0, len(cust_aa) - 1)
#                 aa = cust_aa[rand_ind]
#             else:
#                 rand_ind = rand.randint(0, len(nat_aa) - 1)
#                 aa = nat_aa[rand_ind]
#
#             if i == 1:  # start peptide
#                 peptide = aa
#                 continue
#
#             patt = Chem.MolFromSmarts('NCC(=O)O')
#             aa_match = aa.GetSubstructMatches(patt)
#             peptide_match = peptide.GetSubstructMatches(patt)
#
#             # set atom map number for nitrogen connection point of amino acid
#             aa_old_attach_idx = None
#             for pairs in aa_match:
#                 for atom_idx in pairs:
#                     atom = Chem.Mol.GetAtomWithIdx(aa, atom_idx)
#                     if atom.GetSymbol() == 'N':
#                         atom.SetAtomMapNum(1)
#                         aa_old_attach_idx = atom_idx
#
#             # set atom map number for carboxyl carbon connection point of peptide and remove carboxyl oxygen
#             pep_old_attach_idx = None
#             peptide = Chem.RWMol(peptide)
#             for pairs in peptide_match:
#                 for atom_idx in pairs:
#                     atom = Chem.Mol.GetAtomWithIdx(peptide, atom_idx)
#                     if atom.GetSymbol() == 'O' and Chem.Atom.GetTotalNumHs(atom) == 1:
#                         peptide.RemoveAtom(atom_idx)
#                     elif atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) != 2:
#                         atom.SetAtomMapNum(2)
#                         pep_old_attach_idx = atom_idx
#             # peptide = Chem.rdchem.Mol.GetMol(peptide)
#
#             combo = Chem.RWMol(Chem.CombineMols(aa, peptide))
#
#             pep_atom = None
#             aa_atom = None
#             for atom in combo.GetAtoms():
#                 if atom.GetAtomMapNum() == 1:
#                     aa_atom = atom.GetIdx()
#                     atom.SetAtomMapNum(0)
#                     Chem.Mol.GetAtomWithIdx(aa, aa_old_attach_idx).SetAtomMapNum(0)
#                 elif atom.GetAtomMapNum() == 2:
#                     pep_atom = atom.GetIdx()
#                     atom.SetAtomMapNum(0)
#                     Chem.Mol.GetAtomWithIdx(peptide, pep_old_attach_idx).SetAtomMapNum(0)
#
#             combo.AddBond(aa_atom, pep_atom, order=Chem.rdchem.BondType.SINGLE)
#             Chem.SanitizeMol(combo)
#             # peptide = Chem.rdchem.Mol.GetMol(combo)
#             peptide = combo
#         f.write(Chem.MolToSmiles(peptide))
#         f.write('\n')


def get_monomers(fp_in):
    """
    Read monomer SMILES strings from all files defined by fp_in into list and shuffle list.
    """

    monomers = []
    for file in fp_in:
        with open(file, 'r') as f:
            monomers.extend(f.readlines())

    rand.shuffle(monomers)

    return monomers


def get_required(fp_req):
    """
    Read in required monomer SMILES strings from file defined by fp_req into a set.
    """

    required = set()
    with open(fp_req, 'r') as f:
        for smiles in f.readlines():
            required.add(smiles)

    return required


def connect_monomers(monomers):
    """
    Takes a list of monomers and connects them to form a peptide and returns the peptide SMILES string.
    """

    # begin conneting each monomer in n_monomers
    peptide = None
    for i, monomer in enumerate(monomers):

        monomer = Chem.MolFromSmiles(monomer)

        # start peptide with first monomer
        if i == 0:
            peptide = monomer
            continue

        # TODO: make option for passing in patterns for other backbone structures
        patt = Chem.MolFromSmarts('NCC(=O)O')   # alpha amino acid backbone
        monomer_match = monomer.GetSubstructMatches(patt)
        peptide_match = peptide.GetSubstructMatches(patt)

        # set atom map number for nitrogen connection point of monomer
        monomer_old_attach_idx = None
        for pairs in monomer_match:
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(monomer, atom_idx)
                if atom.GetSymbol() == 'N':
                    atom.SetAtomMapNum(1)
                    monomer_old_attach_idx = atom_idx

        # set atom map number for carboxyl carbon connection point of peptide and remove carboxyl oxygen
        pep_old_attach_idx = None
        peptide = Chem.RWMol(peptide)
        for pairs in peptide_match:
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(peptide, atom_idx)
                if atom.GetSymbol() == 'O' and Chem.Atom.GetTotalNumHs(atom) == 1:
                    peptide.RemoveAtom(atom_idx)
                elif atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) != 2:
                    atom.SetAtomMapNum(2)
                    pep_old_attach_idx = atom_idx

        # combine mols to enable modification
        combo = Chem.RWMol(Chem.CombineMols(monomer, peptide))

        # get atom indicies on new combo mol and remove old atom map numbers
        peptide_atom = None
        monomer_atom = None
        for atom in combo.GetAtoms():
            if atom.GetAtomMapNum() == 1:
                monomer_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)
                Chem.Mol.GetAtomWithIdx(monomer, monomer_old_attach_idx).SetAtomMapNum(0)
            elif atom.GetAtomMapNum() == 2:
                peptide_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)
                Chem.Mol.GetAtomWithIdx(peptide, pep_old_attach_idx).SetAtomMapNum(0)

        combo.AddBond(monomer_atom, peptide_atom, order=Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(combo)
        peptide = combo

    return Chem.MolToSmiles(peptide)


def main():
    parser = argparse.ArgumentParser(
        description='Connects [n] number of monomers to form a peptide. Takes input file(s) containing '
        'monomer SMILES strings and outputs to a text file the peptides as SMILES strings.')
    parser.add_argument('num_monomers', type=int, help='The number of monomers per peptide.')
    parser.add_argument('-n', '--num_pep', dest='num_peptides', type=int,
                        default=None, help='The number of peptides to make. Defaults to all permutations.')
    parser.add_argument('-i', '--in', dest='in_file', nargs='+',
                        default=['custom_S.txt'], help='The input text file(s) containing monomer SMILES strings')
    parser.add_argument('-o', '--out', dest='out_file', default='demo_10.txt',
                        help='The output text file to write the peptide SMILES strings')
    parser.add_argument('-r', '--req', dest='required', nargs='+', default=['required.txt'],
                        help='The text file defining monomers that a peptide is required to have at least one of.')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='smiles/monomers', help='The filepath to the input files')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/peptides',
                        help='The filepath to the output file')

    args = parser.parse_args()

    # get filepath to infile, outfile, and requirement file
    # filepath to base project directory
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.in_file]
    fp_out = str(base_path / args.fp_out / args.out_file)
    fp_req = [str(base_path / args.fp_in / file) for file in args.required]

    print(args.num_monomers)
    print(args.num_peptides)
    print(args.in_file)
    print(args.out_file)
    print(fp_in)
    print(fp_out)

    exit()
    # get all monomers and create all possible "combinations" (actually cartesian product)
    monomers = get_monomers(fp_in)
    mono_combs = product(monomers, args.num_monomers)    # cartesian product

    # create set of required monomers
    required = get_required(fp_req)

    # create peptides
    with open(fp_out, 'w') as f:
        for count, monomers in enumerate(mono_combs):

            # only produce num_peptides peptides
            if count >= args.num_peptides:
                break

            # check to see if list of monomers has at least one required monomer/side chain
            mono_set = set(monomers)
            if not mono_set.intersection(required):
                continue

            # create peptide and write to SMILES to file
            f.write(connect_monomers(monomers))
            f.write('\n')


if __name__ == '__main__':
    rand.seed(time())
    main()
