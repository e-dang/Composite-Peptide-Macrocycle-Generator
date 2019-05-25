import argparse
import json
import multiprocessing
from itertools import islice, product
from pathlib import Path
from random import sample

from rdkit import Chem
from tqdm import tqdm

from utils import ranges

NITROGEN_MAP_NUM = 1
CARBON_MAP_NUM = 2

ALPHA = 'NCC(=O)O'
BETA = 'NCCC(=O)O'

CAPACITY = 500000000


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
                monomers.append((doc['monomer'], doc['type'], doc['required']))

    return monomers, len(monomers)


def generate_peptides(num_monomers, num_peptides, fp_in, num_jobs, job_num, random_sample=False):
    """
    Top level function that initiates peptide generation. Breaks up monomer cross product set into num_jobs number of
    chunks and passes the corresponding chunk to parallelize() for further processing.

    Args:
        num_monomers (int): The number of monomers per peptide
        num_peptides (int): The number of peptides to make
        fp_in (list): A list of filepaths to input files containing monomers
        num_jobs (int): The totoal number of jobs used to generate peptides
        job_num (int): The job number or ID that identifies which chunk of monomer cross product set to convert into
            peptides
        random_sample (bool, optional): Generate peptides randomly from monomer cross procut set with replacement.
            Defaults to False.

    Returns:
        generator: A generator of tuples containing the peptides SMILES string, the monomers' SMILES strings, and the
            type of backbone on N-terminus
    """

    global ITERATIONS

    # specified number of sets of monomers to form into peptide
    monomers, total_monomers = get_monomers(fp_in)
    if random_sample:
        mono_prod = [sample(monomers, num_monomers) for i in range(num_peptides)]   # random sampling
    else:
        start, stop = ranges(total_monomers ** num_monomers, num_jobs)[job_num - 1]
        mono_prod = islice(product(monomers, repeat=num_monomers), start, stop)  # create cartesian product

    ITERATIONS = stop - start if num_peptides is None else num_peptides

    return parallelize(mono_prod)


def parallelize(mono_prod):
    """
    Makes parallel calls to connect_monomers() and returns the results in a single generator. Helper function of
        generate_peptides().

    Args:
        mono_prod (iterable): An iterable of lists containing num_monomers number of tuples, each of which monomer
            SMILES strings, the backbone type, and the monomer's required flag

    Yields:
        tuple: The peptide's SMILES string, a list of monomers' SMILES strings, and the backbone type on the N-terminus
    """
    with multiprocessing.Pool() as pool:
        for peptide, monomers, n_term in pool.imap_unordered(connect_monomers, mono_prod):
            yield peptide, monomers, n_term


def connect_monomers(monomer_data):
    """
    Takes a list of monomer SMILES strings and connects them to form a peptide and returns the peptide SMILES string

    Args:
        monomers (list): List of length args.num_monomers of monomer SMILES strings to be connected into a peptide

    Returns:
        string: The SMILES string representation of the resulting peptide
        string: The type of monomer on the N-terminus of peptide
    """

    # check if list of monomers has at least one required monomer
    monomers, types, required = zip(*monomer_data)
    if True not in required:
        return None, None, None

    # begin conneting each monomer in monomers
    peptide = None
    previous_type = None
    n_term = None
    for i, monomer_tup in enumerate(zip(monomers, types)):

        monomer = Chem.MolFromSmiles(monomer_tup[0])

        # start peptide with first monomer
        if i == 0:
            peptide = monomer
            previous_type = monomer_tup[1]
            n_term = monomer_tup[1]  # need to know what type of a monomer on N-term
            continue

        # choose peptide matching pattern based on previous monomer backbone type
        if previous_type == 'alpha_amino_acid':
            patt_prev = Chem.MolFromSmarts(ALPHA)
        elif previous_type in ('beta3_amino_acid', 'beta2_amino_acid'):
            patt_prev = Chem.MolFromSmarts(BETA)

        # choose monomer matching pattern based on backbone type
        if monomer_tup[1] == 'alpha_amino_acid':
            patt_cur = Chem.MolFromSmarts(ALPHA)
        elif monomer_tup[1] in ('beta3_amino_acid', 'beta2_amino_acid'):
            patt_cur = Chem.MolFromSmarts(BETA)

        previous_type = monomer_tup[1]
        peptide_match = peptide.GetSubstructMatches(patt_prev)
        monomer_match = monomer.GetSubstructMatches(patt_cur)

        # set atom map number for nitrogen connection point of monomer
        monomer_old_attach_idx = None
        for pairs in monomer_match:
            for atom_idx in pairs:
                atom = Chem.Mol.GetAtomWithIdx(monomer, atom_idx)
                if atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) != 0:
                    atom.SetAtomMapNum(NITROGEN_MAP_NUM)
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
                    atom.SetAtomMapNum(CARBON_MAP_NUM)
                    pep_old_attach_idx = atom_idx

        # combine mols to enable modification
        combo = Chem.RWMol(Chem.CombineMols(monomer, peptide))

        # get atom indicies on new combo mol and remove old atom map numbers
        peptide_atom = None
        monomer_atom = None
        for atom in combo.GetAtoms():
            if atom.GetAtomMapNum() == NITROGEN_MAP_NUM:
                monomer_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)
                Chem.Mol.GetAtomWithIdx(monomer, monomer_old_attach_idx).SetAtomMapNum(0)
            elif atom.GetAtomMapNum() == CARBON_MAP_NUM:
                peptide_atom = atom.GetIdx()
                atom.SetAtomMapNum(0)
                Chem.Mol.GetAtomWithIdx(peptide, pep_old_attach_idx).SetAtomMapNum(0)

        try:
            combo.AddBond(monomer_atom, peptide_atom, order=Chem.rdchem.BondType.SINGLE)
            Chem.SanitizeMol(combo)
            peptide = combo
        except ValueError:
            print('Error!')
            print('Monomer:', Chem.MolToSmiles(monomer))
            print('Peptide:', Chem.MolToSmiles(peptide))
            print('Combo:', Chem.MolToSmiles(combo) + '\n')

    return Chem.MolToSmiles(peptide), monomers, n_term


def write_data(peptide_generator, num_peptides, job_num, fp_out, show_progress=False):
    """
    Parses data from peptide generator object and writes it to file. File names are based on num_peptides,
    template type, job_num, and dump_count. This function dumps data into a file when either num_peptides is reached,
    generator object has been depleted, or length of list containing parsed data reaches CAPACITY .

    Args:
        peptide_generator (generator): Generator object containing peptide SMILES strings, a list of monomers' SMILES
            strings, and backbone type on N-terminus
        num_peptides (int): The number of peptides being created
        job_num (int): The job number or ID
        fp_out (str): The filepath to the output file, where the file name is incomplete (it is finished during
            execution of this function)
        show_progress (bool, optional): Show progress bar. Defaults to False.

    Returns:
        Bool: True if function completes
    """

    def write_json(collection, fp_out, dump_count):
        with open(fp_out + str(dump_count) + '.json', 'w') as f:
            json.dump(collection, f)

    collection = []
    dump_count = 1
    # for data in peptide_generator:
    for count, data in enumerate(tqdm(peptide_generator, total=ITERATIONS, desc='Peptides: ', disable=show_progress)):

        peptide, monomers, n_term = data[0], data[1], data[2]
        # result from peptide without any required monomers
        if None in (peptide, monomers, n_term):
            continue

        # combine data
        doc = {}
        doc['peptide'], doc['monomers'], doc['N-term'] = peptide, monomers, n_term
        collection.append(doc)

        if num_peptides is not None and len(collection) >= num_peptides:  # early break
            break
        elif len(collection) >= CAPACITY or count > ITERATIONS:  # early dump due to large list
            write_json(collection, fp_out, dump_count)
            collection = []
            dump_count += 1
            print('Dump:', dump_count, 'Job Num:', job_num)

    # write SMILES to json
    write_json(collection, fp_out, dump_count)
    print('Complete! Job Num:', job_num)

    return True


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
                                 'natural_L.json', 'modified_prolines_CCW.json', 'modified_prolines_CW.json'],
                        help='The input json file(s) containing monomer SMILES strings')
    parser.add_argument('-o', '--out', dest='out_file', default='length',
                        help='The output json file to write the peptide SMILES strings')
    parser.add_argument('--fp_in', dest='fp_in', default='smiles/monomers',
                        help='The filepath to the input files relative to base project directory')
    parser.add_argument('--fp_out', dest='fp_out', default='smiles/peptides',
                        help='The filepath to the output file relative to base project directory')
    parser.add_argument('--num_jobs', dest='num_jobs', default=1, type=int,
                        help='The number of jobs to run on a job array.')
    parser.add_argument('--job_num', dest='job_num', default=1, type=int, help='The job number or ID.')
    parser.add_argument('--random_sample', dest='random_sample', action='store_true',
                        help='Generate peptides randomly (defaults to False)')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar (defaults to False)')

    args = parser.parse_args()

    # format output file based on specified length of peptides
    if args.out_file == 'length':
        args.out_file += str(args.num_monomers)
        args.out_file += '_all' if args.num_peptides is None else '_' + str(args.num_peptides)
        args.out_file += '_' + str(args.job_num) + '_'

    # get full filepath to input, output, and requirement files
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.in_file]
    fp_out = str(base_path / args.fp_out / args.out_file)

    # generate peptides and write to file
    peptide_generator = generate_peptides(args.num_monomers, args.num_peptides, fp_in, args.num_jobs, args.job_num,
                                          args.random_sample)
    write_data(peptide_generator, args.num_peptides, args.job_num, fp_out, args.progress)


if __name__ == '__main__':
    main()
