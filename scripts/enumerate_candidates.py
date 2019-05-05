import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from utils import Database


def generate_candidates(reactant, rxn_docs, verbose=False, show=False):
    """
    Applys all reaction templates to a reactant and collects all unique products

    Args:
        reactant (str): The reactant SMILES string
        rxn_docs (pymongo cursor): The collection of reaction documents containing reaction template SMARTS strings
        verbose (bool, optional): Print all results to console. Defaults to False.
        show (bool, optional): Display image of all products. Defaults to False.

    Returns:
        (list, set): Products is a list of SMARTS strings of all unique products, and reacting_side_chains is a set of
        unique side chains that reacted in order to form the products
    """

    reactant = Chem.MolFromSmiles(reactant)

    # apply reactions
    unique_prod = {}
    reacting_side_chains = set()
    for doc in rxn_docs:
        rxn = AllChem.ReactionFromSmarts(doc['reaction_smarts'])
        prod = rxn.RunReactants((reactant,))

        # get unique products
        for mol in prod:
            reacting_side_chains.add(doc['sc_smiles'])
            Chem.SanitizeMol(mol[0])
            smiles = Chem.MolToSmiles(mol[0])
            unique_prod[smiles] = mol[0]

    products = sorted(unique_prod.keys())

    if verbose:
        for mol in unique_prod.values():
            print(mol)

    if show:
        for mol in unique_prod.values():
            Draw.MolToImage(mol).show()

    return products, reacting_side_chains


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', '--fin', dest='in_file', nargs='+', default=['length3_all.txt'],
                        help='The input text file(s) containing the SMILES strings of template-peptide connected molecule')
    parser.add_argument('-fp', '--fp_in', dest='fp_in', default='/smiles/template_peptide/c_term/',
                        help='The filepath to the input files relative to base project directory')
    parser.add_argument('-d', '--db', dest='database', default='rxn_templates',
                        help='The mongoDB database to connect to')
    parser.add_argument('-hn', '--host', dest='host', default='localhost',
                        help='The host MongoDB server to connect to')
    parser.add_argument('-p', '--port', dest='port', type=int, default=27017,
                        help='The port on host server to connect to')
    parser.add_argument('-s', '--store', dest='store', action='store_false', help='Toggle to not store results')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='Toggle to print reaction SMARTS to console')
    parser.add_argument('-sh', '-show', dest='show', action='store_true',
                        help='Toggle to show reaction SMARTS image; Do not toggle if generating a lot of reaction '
                        'templates')

    args = parser.parse_args()

    # generate filepath to the input file
    base_path = Path(__file__).resolve().parents[1]
    fp_in = [str(base_path / args.fp_in / file) for file in args.in_file]

    # establish database connection and retrieve reaction documents
    db = Database(host=args.host, port=args.port, db=args.db)
    rxn_docs = db.find_all('reactions')

    # for each reactant generate candidates
    with open(fp_in, 'r') as f:
        for line in f.readlines():

            # extract data
            line = line.split(',')
            reactant = line[0]
            peptide = line[1]
            template = line[2]
            monomers = line[3:]

            products, reacting_side_chains = generate_candidates(
                reactant, rxn_docs, verbose=args.verbose, show=args.show)

            if args.store:
                db.insert_candidates(reactant, products, len(products), peptide,
                                     template, monomers, list(reacting_side_chains))


if __name__ == '__main__':
    main()
