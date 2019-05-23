import argparse
import json
from itertools import chain
from pathlib import Path

from pymongo.errors import DuplicateKeyError
from rdkit import Chem
from tqdm import tqdm

from utils import Database, read_mols

RXN_MAP_NUM = 2  # atom map number for atom participating in EAS reaction
ATT_MAP_NUM = 4  # atom map number for attachment point of side chain to peptide backbone
ETHYL = 1   # modifcation key
PROPYL = 2  # modifcation key


def get_side_chains(fp_in):
    """
    Retrieve side chain structures from all files in fp_in and compress them into generator object

    Args:
        fp_in (list): List containing all absolute filepath strings to the input file(s)

    Returns:
        list: A list of tuples containing a side chain as rdkit Mols and the modification codes performed on
        the side chain
    """

    side_chains = []
    for file in fp_in:
        with open(file, 'r') as f:
            for doc in json.load(f):
                side_chains.append((Chem.MolFromSmiles(doc['side_chain']), doc['modification']))

    return side_chains


def set_atom_map_nums(side_chain):
    """
    Enumerate all possible regioisomer locations on each sidechain for EAS reactions and all possible attachment points
    of side chain to peptide backbone. Possible regioisomer locations are determined by:
        Carbon - has at least 1 hydrogen, is aromatic, and atom is not designated as peptide backbone attachment point
        Nitrogen - has at least 1 hydrogen
        Oxygen - is an alcohol/can't be an ether

    Args:
        side_chain (rdkit Mol): The rdkit Mol representation of the side chain

    Returns:
        dict: Keys are SMILES strings, and values are reacting atom index. Each SMILES string is a different regioisomer
            or has a different peptide backbone attachment point
    """

    smiles_idx = {}

    # assign atom map number for attachment and reacting atom of side chain
    patt = Chem.MolFromSmarts('[CH3]*')  # methyl carbon
    matches = side_chain.GetSubstructMatches(patt, useChirality=False)
    for pairs in matches:
        old_atom_idx = None
        for atom_idx in pairs:
            atom = Chem.Mol.GetAtomWithIdx(side_chain, atom_idx)
            if atom.GetSymbol() == 'C' and Chem.Atom.GetTotalNumHs(atom) == 3:
                atom.SetAtomMapNum(ATT_MAP_NUM)
                old_atom_idx = atom_idx

        # set atom map numbers for every possible EAS reacting atom
        for atom in side_chain.GetAtoms():
            atom_idx = None
            valid = False
            if atom.GetSymbol() == 'C' and atom.GetAtomMapNum() != 4 and Chem.Atom.GetTotalNumHs(atom) != 0 and Chem.Atom.GetIsAromatic(atom):
                atom.SetAtomMapNum(RXN_MAP_NUM)
                atom_idx = atom.GetIdx()
                valid = True
            elif atom.GetSymbol() == 'N' and Chem.Atom.GetTotalNumHs(atom) != 0:
                atom.SetAtomMapNum(RXN_MAP_NUM)
                atom_idx = atom.GetIdx()
                valid = True
            elif atom.GetSymbol() == 'O' and Chem.Atom.GetTotalNumHs(atom) == 1 or atom.GetSymbol() == 'S' and Chem.Atom.GetTotalNumHs(atom) == 1:
                atom.SetAtomMapNum(RXN_MAP_NUM)
                atom_idx = atom.GetIdx()
                valid = True

            if valid:

                # format SMILES string
                mol_str = Chem.MolToSmiles(side_chain)
                ind = mol_str.find(f'[CH3:{ATT_MAP_NUM}]') + 7  # length of atom map string
                mol_str = mol_str[:ind] + f'[*:{ATT_MAP_NUM}]' + mol_str[ind:]
                mol_str = mol_str.replace(f'[CH3:{ATT_MAP_NUM}]', 'C')

                atom.SetAtomMapNum(0)   # reset reacting atom map number

                try:
                    Chem.MolFromSmiles(mol_str)
                    smiles_idx[mol_str] = atom_idx
                except:
                    print('Error: can not convert SMILES string to mol:')
                    print('Side Chain SMILES - ' + Chem.MolFromSmiles(side_chain) + '\nAtom Mapped SMILES - ' + mol_str)

        Chem.Mol.GetAtomWithIdx(side_chain, old_atom_idx).SetAtomMapNum(0)  # reset attachment atom map number

    return smiles_idx


def main():
    global COUNT
    parser = argparse.ArgumentParser(description='Enumerates all reacting atoms (regioisomers) of each side chain using'
                                     ' atom map numbers. Atom map number 2 corresponds to atom that participates in the electrophilic aromatic '
                                     'substitution reaction that closes the macrocycle ring. Atom map number 4 (wildcard atom) corresponds to the '
                                     'connection point between the side chain and the peptide backbone. Stores resulting SMILES strings in MongoDB database by defualt.')
    parser.add_argument('-i', '--in', dest='fin', nargs='+', default=['side_chains.json'],
                        help='The input json file(s) containing side chain structures')
    parser.add_argument('-fi', '--fp_in', dest='fp_in', default='smiles/pre_monomer',
                        help='The filepath to the input files relative to base project directory')
    parser.add_argument('-d', '--db', dest='database', default='rxn_templates',
                        help='The mongoDB database to connect to')
    parser.add_argument('-hn', '--host', dest='host', default='localhost',
                        help='The host MongoDB server to connect to')
    parser.add_argument('-p', '--port', dest='port', type=int, default=27017,
                        help='The port on host server to connect to')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')

    args = parser.parse_args()

    # set up full filepath to input file
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.fin]

    # establish database connection
    db = Database(host=args.host, port=args.port, db=args.database)

    # enumerate all regioisomers with atom map numbers and write to output
    for side_chain in tqdm(get_side_chains(fp_in), disable=args.progress):

        # skip enumeration for ethyl, propyl attachment points, only need to do methyl attachment point
        if ETHYL in side_chain[1] or PROPYL in side_chain[1]:
            continue

        for smiles, atom_idx in set_atom_map_nums(side_chain[0]).items():
            try:
                db.insert_sidechain(Chem.MolToSmiles(side_chain[0]), smiles, ATT_MAP_NUM,
                                    RXN_MAP_NUM, atom_idx)
            except DuplicateKeyError:
                print('Duplicate:', smiles)
                continue


if __name__ == '__main__':
    main()
