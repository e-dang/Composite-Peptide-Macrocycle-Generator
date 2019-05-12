import json
from itertools import chain
from pathlib import Path

from pymongo import MongoClient
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


class Database():
    """
    A class to establish a connection to the MongoDB database

    Returns:
        Database: An instance of self
    """

    def __init__(self, host='localhost', port=27017, client=None, db='rxn_templates'):
        """
        Constructor - initializes database connection

        Args:
            host (str, optional): The server host name. Defaults to 'localhost'.
            port (int, optional): the port number. Defaults to 27017.
            client (pymongo mongoclient, optional): A preinitialized pymongo client. Defaults to None.
            db (str, optional): The database to connect to. Defaults to 'rxn_templates'.
            verbose (bool, optional): Prints all collections within database. Defaults to False.
        """

        self.client = MongoClient(host, port) if client is None else client
        self.db = self.client[db]

    def __del__(self):
        """
        Destructor - properly close connection
        """
        self.client.close()

    def insert(self, collection, document):
        """
        Insert a document into a collection

        Args:
            collection (str): The collection to insert the document into
            document (dict): Dictionary containing the data in attribute: value format

        Returns:
            bool: True if successful
        """

        col = self.db[collection]

        if isinstance(document, list):
            result = col.insert_many(document)
        else:
            result = col.insert_one(document)

        return result

    def find(self, collection, query, projection=None):
        """
        Fetch a set of documents based on query

        Args:
            collection (str): The collection to search
            query (dict): Dictionary containing the query in attribute: value format

        Returns:
            pymongo cursor: The results of the query
        """

        return self.db[collection].find(query) if projection is None else self.db[collection].find(query, projection)

    def find_all(self, collection):
        """
        Fetches all documents in a collection

        Args:
            collection (str): The collection to query

        Returns:
            pymongo cursor: The results of the query
        """

        return self.db[collection].find({})

    def insert_sidechain(self, smiles, atom_mapped_smiles, chain_map_num, rxn_map_num, atom_idx, collection='side_chains'):
        """
        Insert a new side_chain document into the database's side_chains collection

        Args:
            smiles (str): The side chain's SMILES string
            atom_mapped_smiles (str): The side chain's atom mapped SMILES string to be used for generating reaction templates
            chain_map_num (int): The atom map number of the atom connecting to the peptide backbone
            rxn_map_num (int): The atom map number of the atom reacting in the reaction template
            atom_ind (int): The atom index of the reacting atom
            collection (str, optional): A collection name to insert into. Defaults to 'side_chains'.

        Returns:
            bool: True if successful
        """

        # check if connected to correct database
        if self.db.name != 'rxn_templates':
            print('Not connected to rxn_templates database.')
            return

        return self.db[collection].insert_one({'smiles': smiles, 'atom_mapped_smiles': atom_mapped_smiles,
                                               'chain_map_num': chain_map_num, 'rxn_map_num': rxn_map_num, 'atom_idx': atom_idx})

    def insert_reaction(self, template, side_chain, reaction_smarts, collection='reactions'):
        """
        Insert a new reaction template document into the database's reaction collection

        Args:
            temp_name (str): The name of the template
            template (str): The template's SMILES string
            side_chain (str): The side chain's SMILES string
            reaction_smarts (str): The reaction SMARTS string
            collection (str): The name of the collection to insert into. Defaults to 'reactions'.

        Returns:
            bool: True if successful
        """

        # check if connected to correct database
        if self.db.name != 'rxn_templates':
            print('Not connected to rxn_templates database.')
            return

        return self.db[collection].insert_one({'template': template, 'side_chain': side_chain,
                                               'reaction_smarts': reaction_smarts})

    def insert_candidates(self, reactant, products, num_products, template, peptide, monomers,
                          reacting_side_chains, collection='candidates'):
        """
        Insert the result of applying a reaction template to a reactant into the database's candidates collection

        Args:
            reactant (str): The reactant SMILES string
            products (list): A list of all product SMILES strings
            num_products (int): The number of total products enumerated
            temp_name (str): The name of the template in the reactant
            template (str): The template's SMILES string
            peptide (str): The peptide's SMILES string
            monomers (list): A list of the monomers that compose the peptide as SMILES strings
            collection (str, optional): The collection to insert into. Defaults to 'candidates'.

        Returns:
            bool: True if successful
        """

        # check if connected to correct database
        if self.db.name != 'rxn_templates':
            print('Not connected to rxn_templates database.')
            return

        return self.db[collection].insert_one({'reactant': reactant, 'products': products, 'num_products': num_products,
                                               'peptide': peptide, 'template': template,
                                               'monomers': monomers, 'reacting_side_chains': reacting_side_chains})


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


def read_multiple_sdf(filepaths, verbose=False):

    mols = []
    for file in filepaths:
        mols = chain(mols, read_mols(file))

    return mols


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


def create_monomer_requirements(fp='smiles/monomers/required.json'):
    """
    Takes all required monomers from the monomers collection in the molecules database and writes them to a json file.

    Args:
        fp (str, optional): The filepath to the output json file. Defaults to '/smiles/monomers/required.json'.
    """

    # get all required monomers from database
    db = Database(db='molecules')
    required = db.find('monomers', {'required': True}, {'_id': 0})
    collection = []
    for monomer in required:
        collection.append(monomer)
    # print(collection)
    # write monomers to file
    fp = str(Path(__file__).resolve().parents[1] / fp)
    with open(fp, 'w') as f:
        json.dump(collection, f)
