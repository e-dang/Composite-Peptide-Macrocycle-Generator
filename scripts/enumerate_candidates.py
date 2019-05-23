import argparse
import json
from pathlib import Path
from copy import deepcopy

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from tqdm import tqdm

from utils import Database

# TODO: Can make enumeration process quicker by only searching for reaction templates that correspond to side chains in
# the peptide. Use monomer information in template_peptide documents to traceback side_chain structure and query db to
# only find reaction templates that contain that side chain.


def generate_candidates(reactant, rxn_docs):
    """
    Applys all reaction templates to a reactant and collects all unique products

    Args:
        reactant (str): The reactant SMILES string
        rxn_docs (pymongo cursor): The collection of reaction documents containing reaction template SMARTS strings
        verbose (bool, optional): Print all results to console. Defaults to False.
        show (bool, optional): Display image of all products. Defaults to False.

    Returns:
        list: A list of SMARTS strings of all unique products
        list: A list of side chains that reacted in order to form the products
        list: A list of reacting atom indices in the side_chains that formed products
    """

    # print(reactant)
    reactant = Chem.MolFromSmiles(reactant)

    # apply reactions
    reaction_info = {}
    for doc in rxn_docs:
        # print(doc)
        rxn = AllChem.ReactionFromSmarts(doc['reaction_smarts'])
        prod = rxn.RunReactants((reactant,))
        # print(prod)

        # get unique products
        for mol in prod:
            Chem.SanitizeMol(mol[0])
            smiles = Chem.MolToSmiles(mol[0])
            # print(smiles)
            if smiles not in reaction_info.keys():
                reaction_info[smiles] = (doc['side_chain'], doc['atom_idx'])

    # print(reaction_info)
    unique_prods, sc_and_idx = zip(*[(key, value) for key, value in reaction_info.items()])
    reacting_side_chains, atom_idxs = zip(*list(sc_and_idx))
    # print(reacting_side_chains)
    # print(atom_idxs)

    return unique_prods, reacting_side_chains, atom_idxs


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', '--fin', dest='in_file', nargs='+', default=['length3_temp1.json'],
                        help='The input json file(s) containing the SMILES strings of template-peptide connected molecule')
    parser.add_argument('-fp', '--fp_in', dest='fp_in', default='smiles/template_peptide/c_term/',
                        help='The filepath to the input files relative to base project directory')
    parser.add_argument('-d', '--db', dest='database', default='rxn_templates',
                        help='The mongoDB database to connect to')
    parser.add_argument('-hn', '--host', dest='host', default='localhost',
                        help='The host MongoDB server to connect to')
    parser.add_argument('-p', '--port', dest='port', type=int, default=27017,
                        help='The port on host server to connect to')
    parser.add_argument('-s', '--store', dest='store', action='store_false', help='Toggle to not store results')
    parser.add_argument('--show_progress', dest='progress', action='store_false',
                        help='Show progress bar. Defaults to False')

    args = parser.parse_args()

    # generate filepath to the input file
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.in_file]

    # establish database connection and retrieve reaction documents
    db = Database(host=args.host, port=args.port, db=args.database)
    rxn_docs = db.find_all('reactions')

    # for each reactant generate candidates
    db.db = db.client['molecules']
    for fp in fp_in:
        with open(fp, 'r') as f:
            for doc in tqdm(json.load(f), disable=args.progress):

                # extract data
                reactant = doc['template-peptide']
                peptide = doc['peptide']
                template = doc['template']
                monomers = doc['monomers']

                products, reacting_side_chains, atom_idx = generate_candidates(reactant, deepcopy(rxn_docs))

                if args.store:
                    db.insert_candidates(reactant, products, len(products), peptide,
                                         template, monomers, reacting_side_chains, atom_idx)


if __name__ == '__main__':
    main()
