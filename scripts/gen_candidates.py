from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from utils import Database

def generate_candidates(reactant, template=None, sidechain=None, store=False, verbose=False):
    """
    """

    # check if reactant is mol, if not convert
    if isinstance(reactant, str):
        reactant = Chem.MolFromSmiles(reactant)

    # get reaction smarts
    db = Database()
    rxn_docs = db.find_all('reactions')

    # apply reactions
    unique_prod = {}
    for doc in rxn_docs:
        rxn = AllChem.ReactionFromSmarts(doc['smarts'])
        prod = rxn.RunReactants((reactant,))

        for mol in prod:
            Chem.SanitizeMol(mol[0])
            smiles = Chem.MolToSmiles(mol[0])
            unique_prod[smiles] = mol[0]

    products = sorted(unique_prod.keys())

    if verbose:
        for mol in unique_prod.values():
            Draw.MolToImage(mol).show()

    if store:
        db.insert_candidates(Chem.MolToSmiles(reactant), products, len(products), template, sidechain)

    return products

with open('/Users/ericdang/Documents/UCLA_Research/smiles/template_peptide/demo_2.txt', 'r') as f:
    for line in f.readlines():
        generate_candidates(line, store=True)
