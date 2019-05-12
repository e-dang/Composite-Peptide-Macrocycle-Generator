import argparse
import json
import random as rand
from itertools import chain, product
from pathlib import Path
from time import time

from rdkit import Chem
from tqdm import tqdm


def get_monomers(fp_in):
    """
    Read monomer SMILES strings from all files defined by fp_in into list

    Args:
        fp_in (list): List of full filepaths to all input files

    Returns:
        (list): List of tuples, each tuple contains the monomer SMILES string and the monomer type
        (int): The length of the monomers list
    """

    monomers = []
    for file in fp_in:
        with open(file, 'r') as f:
            for doc in json.load(f):
                monomers.append((doc['monomer'], doc['type']))

    return monomers, len(monomers)


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
            for doc in json.load(f):
                required.add(doc['monomer'])

    return required


def connect_monomers(monomers):
    """
    Takes a list of monomer SMILES strings and connects them to form a peptide and returns the peptide SMILES string

    Args:
        monomers (list): List of length args.num_monomers of monomer SMILES strings to be connected into a peptide

    Returns:
        string: The SMILES string representation of the resulting peptide
        string: The type of monomer on the N-terminus of peptide
    """

    # begin conneting each monomer in n_monomers
    peptide = None
    previous_type = None
    peptide_Nterm = None
    for i, monomer_tup in enumerate(monomers):

        monomer = Chem.MolFromSmiles(monomer_tup[0])

        # start peptide with first monomer
        if i == 0:
            peptide = monomer
            previous_type = monomer_tup[1]
            peptide_Nterm = monomer_tup[1]  # need to know what type of a monomer on N-term
            continue

        # choose peptide matching pattern based on previous monomer backbone type
        if previous_type == 'alpha_amino_acid':
            patt_prev = Chem.MolFromSmarts('NCC(=O)O')
        elif previous_type == 'beta3_amino_acid' or previous_type == 'beta2_amino_acid':
            patt_prev = Chem.MolFromSmarts('NCCC(=O)O')

        # choose monomer matching pattern based on backbone type
        if monomer_tup[1] == 'alpha_amino_acid':
            patt_cur = Chem.MolFromSmarts('NCC(=O)O')
        elif monomer_tup[1] == 'beta3_amino_acid' or monomer_tup[1] == 'beta2_amino_acid':
            patt_cur = Chem.MolFromSmarts('NCCC(=O)O')

        previous_type = monomer_tup[1]
        peptide_match = peptide.GetSubstructMatches(patt_prev)
        monomer_match = monomer.GetSubstructMatches(patt_cur)

        # set atom map number for nitrogen connection point of monomer
        monomer_old_attach_idx = None
        for pairs in monomer_match:
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(monomer, atom_idx)
                if atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) != 0:
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

    return Chem.MolToSmiles(peptide), peptide_Nterm


def main():
    parser = argparse.ArgumentParser(
        description='Connects specified number of monomers to form a peptide. Takes input file(s) containing '
        'monomer SMILES strings and outputs to a json file the peptides as SMILES strings. Input and output are json '
        'files because this script will run on a super computer')
    parser.add_argument('num_monomers', type=int, help='The number of monomers per peptide')
    parser.add_argument('-n', '--num_pep', dest='num_peptides', type=int,
                        default=None, help='The number of peptides to make; defaults to the length of the cartesian '
                        'product of all monomers')
    parser.add_argument('-i', '--in', dest='in_file', nargs='+',
                        default=['custom_CCW.json', 'custom_CW.json', 'natural_D.json',
                                 'natural_L.json'],
                        help='The input json file(s) containing monomer SMILES strings')
    parser.add_argument('-o', '--out', dest='out_file', default='length',
                        help='The output json file to write the peptide SMILES strings')
    parser.add_argument('-r', '--req', dest='required', nargs='+',
                        default=['custom_required.json', 'natural_required.json'],
                        help='The json file defining monomers that a peptide is required to have at least one of')
    parser.add_argument('-fi', '--fin', dest='fp_in', default='smiles/monomers',
                        help='The filepath to the input files relative to base project directory')
    parser.add_argument('-fo', '--fout', dest='fp_out', default='smiles/peptides',
                        help='The filepath to the output file relative to base project directory')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')

    args = parser.parse_args()

    # format output file based on specified length of peptides
    if args.out_file == 'length':
        args.out_file += str(args.num_monomers)
        args.out_file += '_all.json' if args.num_peptides is None else '_' + str(args.num_peptides) + '.json'

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
    collection = []
    count = 0
    skipped = 0
    length = length ** args.num_monomers if args.num_peptides is None else args.num_peptides
    for monomers_tup in tqdm(mono_prod, total=length, disable=args.progress):

        # only produce num_peptides peptides
        if args.num_peptides is not None and count >= args.num_peptides:
            break

        # check to see if list of monomers has at least one required monomer
        monomers = [tup[0] for tup in monomers_tup]
        mono_set = set(monomers)
        if not mono_set.intersection(required):
            skipped += 1
            continue

        # create peptide and accumulate data
        doc = {}
        doc['peptide'], doc['N-term'] = connect_monomers(monomers_tup)
        doc['monomers'] = monomers
        collection.append(doc)
        count += 1

        # prevent list from getting too large
        if len(collection) > 500000000:
            with open(fp_out, 'w') as f:
                json.dump(collection, f)
                collection = []

    # write SMILES to json
    print('# Invalid Seqs:', skipped)
    with open(fp_out, 'w') as f:
        json.dump(collection, f)
        print('Complete!')


if __name__ == '__main__':
    main()
