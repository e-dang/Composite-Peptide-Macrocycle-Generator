from pathlib import Path

from pymongo import MongoClient
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from itertools import chain


class Database():
    """
    A class to establish a connection to the MongoDB database

    Returns:
        Database: An instance of self
    """

    def __init__(self, host='localhost', port=27017, client=None, db='rxn_templates', verbose=False):
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

        if verbose:
            print(db.collection_names())

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

    def find(self, collection, query):
        """
        Fetch a set of documents based on query

        Args:
            collection (str): The collection to search
            query (dict): Dictionary containing the query in attribute: value format

        Returns:
            pymongo cursor: The results of the query
        """

        return self.db[collection].find(query)

    def find_all(self, collection):
        """
        Fetches all documents in a collection

        Args:
            collection (str): The collection to query

        Returns:
            pymongo cursor: The results of the query
        """

        return self.db[collection].find({})

    def insert_sidechain(self, smiles, atom_mapped_smiles, chain_ind, rxn_ind, collection='side_chains'):
        """
        Insert a new side_chain document into the database's side_chains collection

        Args:
            smiles (str): The side chain's SMILES string
            atom_mapped_smiles (str): The side chain's atom mapped SMILES string to be used for generating reaction templates
            chain_ind (int): The atom map number of the atom connecting to the peptide backbone
            rxn_ind (int): The atom map number of the atom reacting in the reaction template
            collection (str, optional): A collection name to insert into. Defaults to 'side_chains'.

        Returns:
            bool: True if successful
        """

        return self.db[collection].insert_one({'smiles': smiles, 'atom_mapped_smiles': atom_mapped_smiles,
                                               'chain_ind': chain_ind, 'rxn_ind': rxn_ind})

    def insert_reaction(self, temp_smiles, sc_smiles, reaction_smarts, collection='reactions'):
        """
        Insert a new reaction template document into the database's reaction collection

        Args:
            temp_name (str): The name of the template
            temp_smiles (str): The template's SMILES string
            sc_smiles (str): The side chain's SMILES string
            reaction_smarts (str): The reaction SMARTS string
            collection (str): The name of the collection to insert into. Defaults to 'reactions'.

        Returns:
            bool: True if successful
        """

        return self.db[collection].insert_one({'temp_smiles': temp_smiles, 'sc_smiles': sc_smiles,
                                               'reaction_smarts': reaction_smarts})

    def insert_candidates(self, reactant, products, num_products, temp_smiles, pep_smiles, monomers,
                          reacting_side_chains, collection='candidates'):
        """
        Insert the result of applying a reaction template to a reactant into the database's candidates collection

        Args:
            reactant (str): The reactant SMILES string
            products (list): A list of all product SMILES strings
            num_products (int): The number of total products enumerated
            temp_name (str): The name of the template in the reactant
            temp_smiles (str): The template's SMILES string
            sc_smiles (str): The peptide's SMILES string
            sc_monomers (list): A list of the monomers that compose the peptide as SMILES strings
            collection (str, optional): The collection to insert into. Defaults to 'candidates'.

        Returns:
            bool: True if successful
        """

        return self.db[collection].insert_one({'reactant': reactant, 'products': products, 'num_products': num_products,
                                               'pep_smiles': pep_smiles, 'temp_smiles': temp_smiles,
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
