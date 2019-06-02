from pymongo import MongoClient


class Database():
    """
    A class to establish a connection a MongoDB database

    Returns:
        Database: An instance of self
    """

    def __init__(self, database, host='localhost', port=27017, client=None):
        """
        Constructor - initializes database connection

        Args:
            database (str): The database to connect to.
            host (str, optional): The server host name. Defaults to 'localhost'.
            port (int, optional): the port number. Defaults to 27017.
            client (pymongo mongoclient, optional): A preinitialized pymongo client. Defaults to None.
        """

        self.client = MongoClient(host, port) if client is None else client
        self.db = self.client[database]

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
            projection (dict): Determines which fields to retieve

        Returns:
            pymongo cursor: The results of the query
        """

        return self.db[collection].find(query) if projection is None else self.db[collection].find(query, projection)

    def find_all(self, collection, projection=None):
        """
        Fetches all documents in a collection

        Args:
            collection (str): The collection to query

        Returns:
            pymongo cursor: The results of the query
        """

        return self.db[collection].find({}) if projection is None else self.db[collection].find({}, projection)


class RxnDatabase(Database):

    def __init__(self, database='rxn_templates', host='localhost', port=27017):
        super().__init__(database, host, port)

    def insert_sidechain(self, smiles, atom_mapped_smiles, chain_map_num, rxn_map_num, atom_idx,
                         collection='side_chains'):
        """
        Insert a new side_chain document into the database's side_chains collection

        Args:
            smiles (str): The side chain's SMILES string
            atom_mapped_smiles (str): The side chain's atom mapped SMILES string to be used for generating reaction
                templates
            chain_map_num (int): The atom map number of the atom connecting to the peptide backbone
            rxn_map_num (int): The atom map number of the atom reacting in the reaction template
            atom_idx (int): The atom index of the reacting atom (need for regioSQM filter)
            collection (str, optional): A collection name to insert into. Defaults to 'side_chains'.

        Returns:
            bool: True if successful
        """

        return self.db[collection].insert_one({'smiles': smiles, 'atom_mapped_smiles': atom_mapped_smiles,
                                               'chain_map_num': chain_map_num, 'rxn_map_num': rxn_map_num,
                                               'atom_idx': atom_idx})

    def insert_reaction(self, reaction_smarts, template, side_chain, atom_idx, collection='reactions'):
        """
        Insert a new reaction template document into the database's reaction collection

        Args:
            reaction_smarts (str): The reaction SMARTS string
            temp_name (str): The name of the template
            template (str): The template's SMILES string
            side_chain (str): The side chain's SMILES string
            atom_idx (int): The atom index of the reacting side chain atom
            collection (str): The name of the collection to insert into. Defaults to 'reactions'.

        Returns:
            bool: True if successful
        """

        return self.db[collection].insert_one({'reaction_smarts': reaction_smarts, 'template': template,
                                               'side_chain': side_chain, 'atom_idx': atom_idx})


class MolDatabase(Database):

    def __init__(self, database='molecules', host='localhost', port=27017):
        super().__init__(database, host, port)

    def insert_candidates(self, reactant, products, num_products, template, peptide, monomers,
                          reacting_side_chains, atom_idx, collection='candidates'):
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
            atom_idx (list): A list of indices for each reacting atom in the side chains of the reaction templates that
                produced the corresponding candidates
            collection (str, optional): The collection to insert into. Defaults to 'candidates'.

        Returns:
            bool: True if successful
        """

        return self.db[collection].insert_one({'reactant': reactant, 'products': products, 'num_products': num_products,
                                               'peptide': peptide, 'template': template, 'monomers': monomers,
                                               'reacting_side_chains': reacting_side_chains, 'atom_idx': atom_idx})

    def insert_heterocycles(self, data, collection='heterocycles'):

        return self.db[collection].insert_many(data)

    def insert_filtered_candidates(self, reactant, filter_type, reacting_side_chains, collection='filtered_candidates'):

        return self.db[collection].insert_one({'reactant': reactant, 'filter': filter_type,
                                               'reacting_side_chains': reacting_side_chains})

    def insert_conformers(self, candidate, binary, convergences, energies, collection='conformers'):

        return self.db[collection].insert_one({'candidate': candidate, 'binary': binary,
                                               'num_conformers': len(convergences), 'convergences': convergences,
                                               'energies': energies})
