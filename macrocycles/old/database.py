"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

from collections import OrderedDict
from logging import INFO

from pymongo import ASCENDING, MongoClient
from pymongo.errors import (BulkWriteError, CollectionInvalid,
                            DuplicateKeyError, OperationFailure, ConnectionFailure, InvalidName)

from macrocycles.config import COLLECTIONS, MONGO_SETTINGS, VALIDATORS
from macrocycles.utils.utils import create_logger

LOGGER = create_logger(name=__name__, level=INFO)


class MongoDataBase():
    """
    A class to establish a connection the Mongo database.
    """

    def __init__(self, settings=MONGO_SETTINGS, logger=LOGGER, client=None):
        """
        Constructor - initializes database connection

        Args:
            database (str, optional): The database to connect to. Defaults to 'macrocycles'.
            host (str, optional): The server host name. Defaults to 'localhost'.
            port (int, optional): the port number. Defaults to 27017.
            client (pymongo mongoclient, optional): A pre-initialized pymongo client. Defaults to None.
        """

        try:
            self.logger = logger
            self.client = MongoClient(settings.host, settings.port) if client is None else client
            self.database = self.client[settings.database]
        except (ConnectionFailure, InvalidName, TypeError):
            self.logger.exception(
                f'Settings: host = {settings.host}, port = {settings.port}, database = {settings.database}')
        else:
            self.logger.info(
                f'Established connection to database \'{settings.database}\' with {self.client}')

    def __enter__(self):
        """
        Create object in context manager.
        """

        self.logger.info('Creating instance of MongoDataBase in context manager')
        return self

    def __del__(self):
        """
        Close connection through garbage collection.
        """

        self.logger.info(f'Closing connection with {self.client} through garbage collection')
        self.client.close()

    def __exit__(self, e_type, e_val, traceback):
        """
        Close connection in context manager.
        """

        self.logger.info(f'Closing connection with {self.client} through context manager')
        self.client.close()

    def setup(self, validation_level='moderate', clear=False):
        """
        Set up the database scheme specified in the MongoDataBase section of config.py.

        Args:
            validation_level (str, optional): Set the validation level of each collection. Defaults to 'moderate'.
            clear (bool, optional): If True, clears collections in the database with the same names listed in
                config.COLLECTIONS. Defaults to False.

        Returns:
            bool: True if successful.
        """

        # clear existing collection in database
        if clear:
            for collection in COLLECTIONS:
                self.database[collection].drop()
                self.logger.info(f'Dropped collection \'{collection}\' from database \'{self.database.name}\'')

        # create collections and add validators
        for collection, validator in zip(COLLECTIONS, VALIDATORS):
            try:
                query = OrderedDict([('collMod', collection),
                                     ('validator', validator),
                                     ('validationLevel', validation_level)])
                self.database.create_collection(collection)
                self.database.command(query)
                self.database[collection].create_index([('ID', ASCENDING)], unique=True)
                self.logger.info(
                    f'Created collection \'{collection}\' in database \'{self.database.name}\' and applied validation schema {validator}')
            except (CollectionInvalid, OperationFailure):
                self.logger.exception(f'Variables: collection = {collection}, schema = {validator}')
                break
        else:
            # intialize counter
            self.database[COLLECTIONS.counter].insert_one({'type': 'parent_side_chain', 'count': 0, 'prefix': ''})

            self.logger.info('Successfully completed setup()!')
            return True

        return False

    def insert(self, collection, docs, create_id=False):
        """
        Insert a document into a collection.

        Args:
            collection (str): The collection to insert the document into.
            document (iterable / dict): Iterable of dictionaries or single dictionary containing the data.

        Returns:
            bool: True if successful
        """

        # convert to list
        if not isinstance(docs, list):
            docs = [docs]

        # give ID if needed
        if create_id:
            docs = self.assign_id(docs)

        return self.database[collection].insert_many(docs)

    def assign_id(self, docs):
        """
        Assigns new IDs to each document in docs.

        Args:
            docs (iterable): An iterable containing the documents to be modified as dictionaries.

        Returns:
            iterable: An iterable containing the modified documents.
        """
        # get latest count and prefixes
        counter = self.find(COLLECTIONS.counter, {'type': 'parent_side_chain'}, {'count': 1, 'prefix': 1})[0]
        count, prefix = counter['count'], counter['prefix']

        # assign ID
        for doc in docs:
            doc['ID'], count, prefix = generate_id(count, prefix)

        # update count and prefix
        self.database[COLLECTIONS.counter].find_one_and_update({'type': 'parent_side_chain'},
                                                               {'$set': {'count': count, 'prefix': prefix}})
        self.logger.info(
            f'Successfully created {len(docs)} new IDs and updated the count and prefix to {count}, \'{prefix}\'')

        return docs

    # def insert_parent_side_chains(self, docs, collection=COLLECTIONS.mols, ordered=False):
    #     """
    #     Inserts the

    #     Args:
    #         docs ([type]): [description]
    #         collection ([type], optional): [description]. Defaults to COLLECTIONS.mols.
    #         ordered (bool, optional): [description]. Defaults to False.

    #     Returns:
    #         [type]: [description]
    #     """

    #     # convert docs to list
    #     if not isinstance(docs, list):
    #         docs = [docs]

    #     # assign ID
    #     counter = self.find(COLLECTIONS.counter, {'type': 'parent_side_chain'}, {'count': 1, 'prefix': 1})[0]
    #     count = counter['count']
    #     prefix = counter['prefix']
    #     for doc in docs:
    #         doc['ID'], count, prefix = create_id(count, prefix)

    #     # insert documents
    #     try:
    #         self.database[collection].insert_many(docs, ordered=ordered)
    #         self.database[COLLECTIONS.counter].find_one_and_update({'type': 'parent_side_chain'},
    #                                                                {'$set': {'count': count, 'prefix': prefix}})
    #     except (DuplicateKeyError, ValueError):
    #         self.logger.exception(
    #             f'Error inserting into collection \'{collection}\' in database \'{self.database.name}\'')
    #     except BulkWriteError as err:
    #         self.logger.exception(f'{err.details}')
    #     else:
    #         self.logger.info(
    #             f'Successfully inserted {len(docs)} parent_side_chains into collection \'{collection}\' in database \'{self.database.name}\'')
    #         return True

    #     return False

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

        return self.database[collection].find(query) if projection is None else self.database[collection].find(query, projection)

    def find_all(self, collection, projection=None):
        """
        Fetches all documents in a collection

        Args:
            collection (str): The collection to query

        Returns:
            pymongo cursor: The results of the query
        """

        return self.database[collection].find({}) if projection is None else self.database[collection].find({}, projection)


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
        print('__init__')

    def __enter__(self):
        """
        Open connection in context manager.
        """
        print('__enter__')
        return self

    def __del__(self):
        """
        Close connection
        """
        print('connection closing from __del__')
        self.client.close()

    def __exit__(self, e_type, e_val, traceback):
        """
        Close connection in context manager
        """
        print('connection closing from __exit__')
        self.client.close()

    def insert(self, collection, document):
        """
        Insert a document into a collection

        Args:
            collection (str): The collection to insert the document into
            document (dict, list): List or Dictionary containing the data in attribute: value format

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

    def insert_side_chains(self, data, collection='side_chains'):

        return self.db[collection].insert(data)

    def insert_filtered_candidates(self, reactant, filter_type, reacting_side_chains, collection='filtered_candidates'):

        return self.db[collection].insert_one({'reactant': reactant, 'filter': filter_type,
                                               'reacting_side_chains': reacting_side_chains})

    def insert_conformers(self, candidate, binary, convergences, energies, collection='conformers'):

        return self.db[collection].insert_one({'candidate': candidate, 'binary': binary,
                                               'num_conformers': len(convergences), 'convergences': convergences,
                                               'energies': energies})


def generate_id(count, prefix):
    aa_codes = 'A R N D C G Q E H I L K M F P S T W Y V'.split(' ')  # taken by natural amino acids
    alphabet = 'A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'.split(' ')

    # find next un taken ID
    code = 'A'
    while code in aa_codes:
        try:
            code = prefix + alphabet[count]
        except IndexError:
            count = 0
            prefix = set_prefix(prefix)
        else:
            count += 1

    return code, count, prefix


def set_prefix(prefix):

    # initial prefix assignment
    if prefix == '':
        prefix = 'A'
        return prefix

    # increment last letter of prefix and recursively wrap if necessary
    ending = ord(prefix[-1]) + 1
    if ending > ord('Z'):
        prefix = set_prefix(prefix[:-1]) if len(prefix) > 1 else 'A'
        prefix = prefix + 'A'
        return prefix

    # no recursive wrapping needed
    prefix = prefix[:-1] + str(chr(ending))
    return prefix
