from rdkit import Chem
from rdkit.Chem import AllChem
from pymongo import MongoClient
import os
from rdkit.Chem import Draw
from pathlib import Path


class Database():

    def __init__(self, client=None, db='rxn_templates', verbose=False):

        self.client = MongoClient('localhost', 27017) if client is None else client
        self.db = self.client[db]

        if verbose:
            print(db.collection_names())

    def __del__(self):
        self.client.close()

    def insert(self, collection, document):
        """
        """

        col = self.db[collection]

        if isinstance(document, list):
            result = col.insert_many(document)
        else:
            result = col.insert_one(document)

        return result

    def find(self, collection, query):
        """
        """

        return self.db[collection].find(query)

    def find_all(self, collection):
        """
        """

        return self.db[collection].find({})

    def insert_sidechain(self, smarts, smiles, chain_ind, rxn_ind):
        """
        Insert a new side_chain document into the database's side_chains collection
        """

        return self.db['side_chains'].insert_one({'smarts': smarts, 'smiles': smiles, 'chain_ind': chain_ind, 'rxn_ind': rxn_ind})

    def insert_reaction(self, template, side_chain, smarts):
        """
        Insert a new reaction template document into the database's reaction collection
        """

        return self.db['reactions'].insert_one({'template': template, 'side_chain': side_chain, 'smarts': smarts})

    def insert_candidates(self, reactant, products, num_products, template, sidechain):
        """
        Insert a new reaction template document into the database's candidates collection
        """

        return self.db['candidates'].insert_one({'reactant': reactant, 'products': products, 'num_products': num_products, 'template': template, 'sidechain': sidechain})


def merge(template, side_chain):
    """
    Takes template and side_chain mongo docs, merges the SMILES strings together at the designated reacting site, and returns the merged SMILES string.
    """

    # convert smiles strings to mols
    temp = Chem.MolFromSmiles(template['smiles'])
    sc = Chem.MolFromSmiles(side_chain['smiles'])

    # get atom map number of reacting site
    temp_map_num = int(template['rxn_ind'])
    sc_map_num = int(side_chain['rxn_ind'])

    # remove substruct and combine mols
    temp = Chem.DeleteSubstructs(temp, Chem.MolFromSmiles(template['substruct']))
    combo = Chem.RWMol(Chem.CombineMols(temp, sc))

    # get reacting atom indicies
    temp_atom = None
    sc_atom = None
    for atom in combo.GetAtoms():
        if atom.GetAtomMapNum() == temp_map_num:
            temp_atom = atom.GetIdx()
        elif atom.GetAtomMapNum() == sc_map_num:
            sc_atom = atom.GetIdx()

    # check if reacting atom is a nitrogen and if so remove all hydrogens
    if combo.GetAtomWithIdx(sc_atom).GetSymbol() == 'N':
        combo.GetAtomWithIdx(sc_atom).SetNumExplicitHs(0)

    # create bond
    combo.AddBond(temp_atom, sc_atom, order=Chem.rdchem.BondType.SINGLE)
    Chem.SanitizeMol(combo)

    return Chem.MolToSmiles(combo)


def generate_rxn_temp(template, side_chain, store=False, verbose=False):
    """
    Takes template and side_chain mongo docs and combines them to create and return a SMARTS reaction template.
    """

    # get smiles strings
    temp = template['smiles']
    sc = side_chain['smiles']
    prod = merge(template, side_chain)

    # combine smiles strings
    rxn = '(' + temp + '.' + sc + ')>>' + prod

    if verbose:
        Draw.ReactionToImage(AllChem.ReactionFromSmarts(rxn), subImgSize=(500, 500)).show()

    if store:
        new_rxn = rxn.replace('\\', '\\\\')
        db = Database()
        db.insert('reactions',
                  {"template": template['name'], "side_chain": side_chain['name'], "smiles": new_rxn})

    return rxn


def generate_candidates(reactant, store=False, verbose=True):
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
        rxn = AllChem.ReactionFromSmarts(doc['smiles'])
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
        db.insert('candidates', {'reactant': Chem.MolToSmiles(reactant),
                                 'products': products, 'number_products': len(products)})

    return products


def read_mols(filepath=None, verbose=False):

    # set default
    if filepath is None:
        filepath = str(Path(__file__).resolve().parents[1] / 'chemdraw/test_rxn.sdf')

    mols = Chem.SDMolSupplier(filepath)

    if verbose:
        for mol in mols:
            print(Chem.MolToSmiles(mol))
            Draw.MolToImage(mol).show()

    return mols

# def draw_mols(mols=None, file=None, verbose=False):
#     """
#     """
#     if mols is None and file is None:
#         print('Need mol smiles string or text file to read from')
#         exit()
#
#     with open()


def write_mols(mols, file):
    """
    Write a list of molecules to an sdf file.
    """

    fp = '/Users/ericdang/Documents/UCLA_Research/chemdraw/' + file
    writer = Chem.SDWriter(fp)

    for mol in mols:

        # check if mol is a molecule or smiles string
        if isinstance(mol, str):
            mol = Chem.MolFromSmiles(mol)

        writer.write(mol)
