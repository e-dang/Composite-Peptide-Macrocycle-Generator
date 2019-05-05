import argparse
import random as rand
from itertools import product, chain
from pathlib import Path
from time import time
from tqdm import tqdm

from rdkit import Chem


def get_monomers(fp_in):
    """
    Read monomer SMILES strings from all files defined by fp_in into list, shuffle list, and compress into generator

    Args:
        fp_in (list): List of full filepaths to all input files

    Returns:
        tuple - (generator, int): Generator contains all monomer SMILES strings, and length is the number of monomers
    """

    monomers = []
    length = None
    for file in fp_in:
        with open(file, 'r') as f:
            monomers.extend(list(f.readlines()))

    rand.shuffle(monomers)
    length = len(monomers)
    monomers = chain(monomers)

    return monomers, length


def get_required(fp_req):
    """
    Read in required monomer SMILES strings from file defined by fp_req into a set.

    Args:
        fp_req (list): Peptides are required to contain at least one monomer from the list of required monomers. This
            argument is a list of full filepaths to the files containg the required monomer SMILES strings

    Returns:
        set: The set of required monomer SMILES strings
    """

    required = set()
    for file in fp_req:
        with open(file, 'r') as f:
            for smiles in f.readlines():
                required.add(smiles)

    return required


def connect_monomers(monomers):
    """
    Takes a list of monomer SMILES strings and connects them to form a peptide and returns the peptide SMILES string

    Args:
        monomers (list): List of length args.num_monomers of monomer SMILES strings to be connected into a peptide

    Returns:
        string: The SMILES string representation of the resulting peptide
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
                elif atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 0:
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

        try:
            combo.AddBond(monomer_atom, peptide_atom, order=Chem.rdchem.BondType.SINGLE)
            Chem.SanitizeMol(combo)
            peptide = combo
        except:
            print('Error! Could not sanitize mol-------------')
            print('Monomer:', Chem.MolToSmiles(monomer))
            print('Peptide:', Chem.MolToSmiles(peptide))
            print('Combo:', Chem.MolToSmiles(combo))

    return Chem.MolToSmiles(peptide)


def main():
    parser = argparse.ArgumentParser(
        description='Connects specified number of monomers to form a peptide. Takes input file(s) containing '
        'monomer SMILES strings and outputs to a text file the peptides as SMILES strings. Input and output are text '
        'files because this script will run on a super computer')
    parser.add_argument('num_monomers', type=int, help='The number of monomers per peptide')
    parser.add_argument('-n', '--num_pep', dest='num_peptides', type=int,
                        default=None, help='The number of peptides to make; defaults to the length of the cartesian '
                        'product of all monomers')
    parser.add_argument('-i', '--in', dest='in_file', nargs='+',
                        default=['custom_S.txt', 'custom_R.txt', 'natural_D.txt', 'natural_L.txt'],
                        help='The input text file(s) containing monomer SMILES strings')
    parser.add_argument('-o', '--out', dest='out_file', default='length',
                        help='The output text file to write the peptide SMILES strings')
    parser.add_argument('-r', '--req', dest='required', nargs='+',
                        default=['custom_required.txt', 'natural_required.txt'],
                        help='The text file defining monomers that a peptide is required to have at least one of')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='smiles/monomers',
                        help='The filepath to the input files relative to base project directory')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/peptides',
                        help='The filepath to the output file relative to base project directory')

    args = parser.parse_args()

    # format output file based on specified length of peptides
    if args.out_file == 'length':
        args.out_file += str(args.num_monomers)
        args.out_file += '_all.txt' if args.num_peptides is None else '_' + str(args.num_peptides)

    # get full filepath to input, output, and requirement files
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.in_file]
    fp_out = str(base_path / args.fp_out / args.out_file)
    fp_req = [str(base_path / args.fp_in / file) for file in args.required]

    # get all monomers and create all possible "combinations" (actually cartesian product)
    monomers, length = get_monomers(fp_in)
    mono_prod = product(monomers, repeat=args.num_monomers)    # cartesian product

    # get set of required monomers
    required = get_required(fp_req)

    # create peptides
    with open(fp_out, 'w') as f:
        for count, monomers in tqdm(enumerate(mono_prod), total=length ** args.num_monomers):

            # only produce num_peptides peptides
            if args.num_peptides is not None and count >= args.num_peptides:
                break

            # check to see if list of monomers has at least one required monomer
            mono_set = set(monomers)
            if not mono_set.intersection(required):
                continue

            monomers_str = ''
            for monomer in monomers:
                monomers_str += ',' + monomer.rstrip()

            # create peptide and write to SMILES to file
            f.write(connect_monomers(monomers) + monomers_str)
            f.write('\n')


if __name__ == '__main__':
    rand.seed(time())
    main()
