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


def generate_cross_product(num_monomers, num_peptides, fp_in, mode, fp_out, progress):
    """
    Generate the cross product of all monomers and either write to output json file or return for further processing.

    Args:
        num_monomers (int): The number of monomers per peptide
        num_peptides (int): The number of peptides to create if specified by the user
        fp_in (str): The filepath to the input json file containing all monomers to incorporate into peptides
        mode (int): The execution mode of the program
        fp_out (str): The filepath to the output json file
        progress (bool): Show progress bar

    Returns:
        list: Contains lists of monomers to be connected together to form a peptide
        int: The number of unique monomers that form the peptides
    """

    # get all monomers and create all possible "combinations" (actually cartesian product)
    monomers, length = get_monomers(fp_in)
    mono_prod = product(monomers, repeat=num_monomers)    # cartesian product

    # if mode 0 save cross product results else return results
    if mode == 0:

        collection = []
        length = length ** num_monomers if num_peptides is None else num_peptides
        for count, monomers_tup in tqdm(enumerate(mono_prod), total=length, disable=progress):

            if count >= length:
                break

            collection.append(monomers_tup)

            # prevent list from getting too large
            if len(collection) > 500000000:
                with open(fp_out, 'w') as f:
                    json.dump(collection, f)
                    collection = []

        with open(fp_out, 'w') as f:
            json.dump(collection, f)

    else:
        return mono_prod, length


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


def generate_peptides(mono_prod, num_monomers, num_peptides, length, fp_out, fp_req, progress):
    """
    Connects a set of monomers together to form a peptide and write the results to a json file.

    Args:
        mono_prod (iterable): Contains lists of monomers to connect together to form a peptide
        num_monomers (int): The number of monomers per peptide
        num_peptides (int): The number of peptides to create
        length (int): The number of unique monomers that form the peptides
        fp_out (str): The filepath to write the output json file
        fp_req (str): The filepath to the json file containing the required monomers
        progress (bool): Show progress bar
    """

    # get set of required monomers
    required = get_required(fp_req)

    # create peptides
    collection = []
    count = 0
    skipped = 0
    length = length ** num_monomers if num_peptides is None else num_peptides
    for monomers_tup in tqdm(mono_prod, total=length, disable=progress):

        # only produce num_peptides peptides
        if num_peptides is not None and count >= num_peptides:
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


def main():
    parser = argparse.ArgumentParser(
        description='Connects specified number of monomers to form a peptide. Takes input file(s) containing '
        'monomer SMILES strings and outputs to a json file the peptides as SMILES strings. Input and output are json '
        'files because this script will run on a super computer')
    parser.add_argument('num_monomers', type=int, help='The number of monomers per peptide')
    parser.add_argument('-n', '--num_pep', dest='num_peptides', type=int,
                        default=None, help='The number of peptides to make; defaults to the length of the cartesian '
                        'product of all monomers')
    parser.add_argument('--mode', dest='mode', choices=[0, 1, 2], type=int, default=2, help='Determines if the scripts generates '
                        'cross product of monomers(0), generates peptides(1), or both(2). If 0 or 2, then input file '
                        'needs to be output of generate_monomers.py, if 1 input file needs to be output of '
                        'generate_cross_product().')
    parser.add_argument('-i', '--in', dest='in_file', nargs='+',
                        default=['custom_CCW.json', 'custom_CW.json', 'natural_D.json',
                                 'natural_L.json'],
                        help='The input json file(s) containing monomer SMILES strings if in mode 0 or 2. If in '
                        'mode 1, then this is the path to output file generated by generate_cross_product().')
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
        args.out_file += '_monomers' if args.mode == 0 else ''
        args.out_file += '_all.json' if args.num_peptides is None else '_' + str(args.num_peptides) + '.json'

    # get full filepath to input, output, and requirement files
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file)
             for file in args.in_file] if args.mode in [0, 2] else str(base_path / args.fp_out / args.in_file[0])
    fp_out = str(base_path / args.fp_out / args.out_file)
    fp_req = [str(base_path / args.fp_in / file) for file in args.required]

    # execute program based on mode
    if args.mode == 0:  # generate and write to json file the cross products of all monomers
        generate_cross_product(args.num_monomers, args.num_peptides, fp_in, args.mode, fp_out, args.progress)

    elif args.mode == 1:    # using the file containing cross products of all monomers, generate peptides
        with open(fp_in, 'r') as f:
            mono_prod = json.load(f)
            generate_peptides(mono_prod, args.num_monomers, args.num_peptides, 0, fp_out, fp_req, args.progress)

    else:   # generate both the cross products of monomers and create peptide
        mono_prod, length = generate_cross_product(
            args.num_monomers, args.num_peptides, fp_in, args.mode, fp_out, args.progress)
        generate_peptides(mono_prod, args.num_monomers, args.num_peptides, length, fp_out, fp_req, args.progress)


if __name__ == '__main__':
    main()
