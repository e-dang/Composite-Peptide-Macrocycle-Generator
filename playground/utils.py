"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

import inspect
import json
import logging
import logging.handlers
import os
from collections import OrderedDict, namedtuple
from itertools import chain
from pathlib import Path
from copy import deepcopy
from bson import json_util, errors

from pymongo import ASCENDING, MongoClient
from pymongo.errors import (BulkWriteError, CollectionInvalid,
                            ConnectionFailure, DuplicateKeyError, InvalidName,
                            OperationFailure)
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import Draw
from macrocycles.config import (COLLECTIONS, DATA_DIR, MONGO_SETTINGS, BB_MAP_NUM,
                                PROJECT_DIR, LOG_DIR, VALIDATORS, INDEX, DI_INPUT_DIR)
from macrocycles.exceptions import (SavingMongoError, SavingSQLError,
                                    WritingJsonError, WritingTxtError, InvalidSmilesString)

########################################################################################################################
########################################################################################################################
########################################################################################################################
Flags = namedtuple('Flags', ['json_flag', 'txt_flag', 'mongo_flag', 'sql_flag'])
IOPaths = namedtuple('IOPaths', ['input', 'output'])
MongoParams = namedtuple('MongoParameters', ['input_cols', 'input_types', 'output_col'])


class CustomFormatter(logging.Formatter):
    """
    Formatter for logging. Customizes the exception traceback format.
    """

    def format(self, record):
        """
        Overloaded method for applying custom formats to log records.

        Args:
            record (str): The log record.

        Returns:
            str: The formatted log record.
        """

        result = super(CustomFormatter, self).format(record)
        if record.exc_text:
            result = result.replace('\n', '\n\t')
        return result

    def formatException(self, exc_info):
        """
        Overloaded method for applying custom formats to exception log records.

        Args:
            exc_info (str): The Exception log record.

        Returns:
            str: The exception log record.
        """

        result = super(CustomFormatter, self).formatException(exc_info)
        return result


def create_logger(name, level, path=None):
    """
    Creates a logger with a file handler and CustomFormatter.

    Args:
        name (str): The name of the logger.
        level (int): The level of logs to be recorded by the logger.
        path (str, optional): The path to the log file. Defaults to the LOG_DIR/<caller's name>.log.

    Returns:
        Logger: The logger.
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # set path to log file
    if path is None:
        file = inspect.stack()[-1][1].strip('.py')
        path = os.path.join(LOG_DIR, file + '.log')

    # add file handler
    file_handler = logging.handlers.RotatingFileHandler(path, maxBytes=1000000, backupCount=5)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = CustomFormatter(
        '[%(asctime)-19s] [%(levelname)-8s] [%(name)s - %(funcName)s - %(lineno)d] -- %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # add console handler
    # console_handler = logging.StreamHandler(sys.stdout)
    # console_handler.setLevel(logging.ERROR)
    # console_formatter = logging.Formatter('%(levelname)s %(name)s - %(message)s')
    # console_handler.setFormatter(console_formatter)
    # logger.addHandler(console_handler)

    return logger


LOGGER = create_logger(name=__name__, level=logging.INFO)

########################################################################################################################
########################################################################################################################
########################################################################################################################


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

        self.client.close()
        self.logger.info(f'Closing connection with {self.client} through garbage collection')

    def __exit__(self, e_type, e_val, traceback):
        """
        Close connection in context manager.
        """

        self.logger.info(f'Closing connection with {self.client} through context manager')
        self.client.close()

    def __getitem__(self, collection):
        """
        Overloaded [] operator. Allows for accessing collections from class instances rather than dereferencing the
        database attribute.

        Args:
            collection (str): The collection to access in the database.

        Returns:
            Collection: An instance of the MongoDB collection.
        """

        try:
            return self.database[collection]
        except (TypeError, InvalidName):
            self.logger.exception(f'Invalid access attempt ([] operator) on {self}')
            raise

    def __repr__(self):
        """
        String representation of the class.
        """

        return f'database \'{self.database.name}\' using \'{self.client}\''

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
                self[collection].drop()
                self.logger.info(f'Dropped collection \'{collection}\' from database \'{self.database.name}\'')

        # create collections, add validators and indexes
        for collection, validator, index in zip(COLLECTIONS, VALIDATORS, INDEX):
            try:
                query = OrderedDict([('collMod', collection),
                                     ('validator', validator),
                                     ('validationLevel', validation_level)])
                self.database.create_collection(collection)
                self.database.command(query)

                if len(index) > 1:
                    for ind in index:
                        self[collection].create_index(ind, unique=True)
                else:
                    self[collection].create_index(index, unique=True)

            except (CollectionInvalid, OperationFailure):
                self.logger.exception(f'Variables: collection = {collection}, schema = {validator}')
                break
            else:
                self.logger.info(
                    f'Created collection \'{collection}\' in database \'{self.database.name}\' and applied index '
                    f'{index} to validation schema {validator}')
        else:
            # intialize counter
            self[COLLECTIONS.counter].insert_one({'type': 'parent_side_chain', 'count': 0, 'prefix': ''})

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

        return self[collection].insert_many(docs)

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
            unique_id, count, prefix = generate_id(count, prefix)

            # new monomers will have 'm' in ID field
            if doc['ID'] is not None:
                doc['ID'] += unique_id
            else:
                doc['ID'] = unique_id
            # print('unique id', unique_id)

        # update count and prefix
        self[COLLECTIONS.counter].find_one_and_update({'type': 'parent_side_chain'},
                                                      {'$set': {'count': count, 'prefix': prefix}})
        self.logger.info(
            f'Successfully created {len(docs)} new IDs and updated the count and prefix to {count}, \'{prefix}\'')

        return docs

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

        return self[collection].find(query) if projection is None else \
            self[collection].find(query, projection)


def generate_id(count, prefix):
    """
    Method for generating new IDs for parent_side_chains.

    Args:
        count (int): The count stored in the MongoDataBase, which is used for determining the next letter in the ID.
        prefix (str): The prefix to which the next letter in the ID will be appended to.

    Returns:
        str: The new ID.
    """

    aa_codes = 'A R N D C G Q E H I L K M F P S T W Y V'.split(' ')
    alphabet = 'A B C D E F G H I J K L M N O P Q R S T U V W X Y Z'.split(' ')
    code = 'A'

    # find next untaken ID
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
    """
    Recursive method for rotating the prefix letter once they reach 'Z'. For example, a prefix 'ZZ' will turn into
    'AAA'.

    Args:
        prefix (str): The prefix to be rotated.

    Returns:
        str: The rotated prefix.
    """

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

########################################################################################################################
########################################################################################################################
########################################################################################################################


class Base():
    """
    Class from which all other classes in package will inherit from. Handles data I/O.

    Attributes:
        project_dir: The filepath to the root project directory.
        data_dir: The filepath to the data directory.
        fp_in: The filepath(s) to the input file(s).
        fp_out: The filepath to the output file.
        mongo_db: A connection to the MongoDB where the result_data will be stored.
        mongo_params: A MongoParams namedtuple containing the collection name(s) and value(s) held within the
            'type' field of the documents to be retrieved, as well as the output collection name.
        sql_db: A connection to the SQL database where the result_data will be stored.
        logger: The logger of the child classes' module.
        result_data: A list to contain the result_data.
        input_flags / output_flags: Namedtuples containing the following:
            json_flag: If true, data will be loaded/saved in this format.
            txt_flag: If true, data will be loaded/saved in this format.
            mongo_flag: If true, data will be loaded/saved in this format.
            sql_flag: If true, data will be loaded/saved in this format.
                **Note: Data can only be loaded from one format**
    """

    def __init__(self, paths, mongo_params, logger, input_flags, output_flags):
        """
        Constructor.

        Args:
            paths (IOPaths): A namedtuple containing the I/O paths relative to DATA_DIR.
            mongo_params (MongoParams): A namedtuple containing the collection name(s) and value(s) held within the
                    'type' field of the documents to be retrieved, as well as the output collection name.
            logger (Logger): The logger to be used by this class.
            output_flags (Flags, optional): A namedtuple specifying output formats.
            input_flags (Flags, optional): A namedtuple specifying the input format.
        """

        # I/O
        self.project_dir = PROJECT_DIR
        self.data_dir = DATA_DIR
        self.fp_in = [os.path.join(DATA_DIR, path) for path in paths.input]
        self.fp_out = os.path.join(DATA_DIR, paths.output)
        self.mongo_db = MongoDataBase(logger=logger) if mongo_params is not None else None
        self.mongo_params = mongo_params
        self.sql_db = None
        self.logger = logger

        # data
        self.result_data = []

        # flags
        self.input_flags = input_flags
        self.output_flags = output_flags

    def save_data(self):
        """
        Top level function for saving the data stored in self.result_data to all specified data formats. Calls helper
        function self.save_all_formats().

        Returns:
            bool: True if successful.
        """

        try:
            self.save_all_formats()
        except (WritingJsonError, WritingTxtError, SavingMongoError, SavingSQLError):
            return False
        else:
            if True not in self.output_flags:
                self.logger.warning('No data saved. Please set a data flag specifying which data format to save to.')
                return False
            return True

    def save_all_formats(self):
        """
        Helper function of self.save_data(). Sequentially calls the following functions: self.write_json(),
        self.write_txt(), self.save_mongo_db(), and self.save_sql_db(), in order to specifically diagnose any errors
        raised during exectution of each of these functions.

        Raises:
            WritingJsonError: Failure to write data in json format.
            WritingTxtError: Failure to write data in txt format.
            SavingMongoError: Failure to save data to the Mongo database.
            SavingSQLError: Failure to save data to the SQL database.
        """

        if self.output_flags.json_flag:
            try:
                self.write_json()
            except TypeError:
                self.logger.exception(
                    f'Failed to write result_data to a json file. self.result_data = {self.result_data[:10]}')
                raise WritingJsonError
            else:
                self.logger.info(
                    f'Successfully wrote {len(self.result_data)} data points to file {self.fp_out}.json!')

        if self.output_flags.txt_flag:
            try:
                self.write_txt()
            except TypeError:
                self.logger.exception(
                    f'Failed to write result_data to a txt file. self.result_data = {self.result_data[:10]}')
                raise WritingTxtError
            else:
                self.logger.info(
                    f'Successfully wrote {len(self.result_data)} data points to file {self.fp_out}.txt!')

        if self.output_flags.mongo_flag:
            try:
                self.save_mongo_db()
            except (DuplicateKeyError, ValueError, TypeError, errors.InvalidDocument):
                self.logger.exception(
                    f'Failed to save result_data to the Mongo database. self.result_data = {self.result_data[:10]}')
                raise SavingMongoError
            except BulkWriteError as err:
                self.logger.exception(
                    f'Failed to save result_data to the Mongo database. self.result_data = {self.result_data[:10]}'
                    f'\n{err.details}')
                raise SavingMongoError
            else:
                self.logger.info(f'Successfully saved {len(self.result_data)} data points to the collection '
                                 f'{self.mongo_params.output_col} on {self.mongo_db}')

        if self.output_flags.sql_flag:
            try:
                self.save_sql_db()
            except Exception:
                self.logger.exception(
                    f'Failed to save result_data to the SQL database. self.result_data = {self.result_data[:10]}')
                raise SavingSQLError
            else:
                self.logger.info(f'Successfully saved {len(self.result_data)} data points to the SQL database!')

    def write_json(self):
        """
        Writes data stored in self.result_data to a json file specified by self.fp_out.
        """

        with open(self.fp_out + '.json', 'w') as f:
            json.dump(json.loads(json_util.dumps(self.result_data)), f)

    def write_txt(self):
        """
        Writes data stored in self.result_data to a txt file specified by self.fp_out.
        """

        with open(self.fp_out + '.txt', 'w') as f:

            # write data properties
            f.write(','.join(self.result_data[0].keys()) + '\n')

            # write data
            for data_point in self.result_data:
                data_string = [str(field) for field in data_point.values()]
                f.write(','.join(data_string) + '\n')

    def save_mongo_db(self, create_id=False):
        """
        Saves data stored in self.result_data to the collection defined by self.mongo_params.output_col in the data base
        self.mongo_db.

        Returns:
            bool: True if successful.
        """

        return self.mongo_db.insert(self.mongo_params.output_col, self.result_data, create_id=create_id)

    def save_sql_db(self):
        """
        Saves data stored in self.result_data to the SQL database specified by self.sql_db.

        Returns:
            bool: True if successful.
        """

        return self.sql_db.insert(self.result_data)

    def load_data(self):
        """
        Method to be overloaded by derived classes. Loads data from format specified by self.input_flags in generator
        form. Derived class should unpack the returned list in the overloaded version of this method according to the
        order of files stored in self.fp_in.

        Returns:
            list: A list of generators, where each generator corresponds to data stored in one file.
        """

        if sum(self.input_flags) > 1:
            self.logger.warning(f'No data loaded. Multiple input flags were set, please specify one. self.input_flags: '
                                f'{self.input_flags}')
        elif True not in self.input_flags:
            self.logger.warning(f'No data loaded. No input flags were set, please specify one. self.input_flags: '
                                f'{self.input_flags}')

        data = []
        if self.input_flags.json_flag:
            try:
                for filepath in self.fp_in:
                    data.append(load_json(filepath))
            except (OSError, json.JSONDecodeError):
                self.logger.exception(f'Failed to load the json file(s) {self.fp_in}')
            else:
                self.logger.info(f'Successfully loaded the json file(s) {self.fp_in}')

        if self.input_flags.txt_flag:
            try:
                for filepath in self.fp_in:
                    data.append(load_txt(filepath))
            except OSError:
                self.logger.exception(f'Failed to load the txt file(s) {self.fp_in}')
            else:
                self.logger.info(f'Successfully loaded the txt file(s) {self.fp_in}')

        if self.input_flags.mongo_flag:
            try:
                for input_col, doc_type in zip(*[self.mongo_params.input_cols, self.mongo_params.input_types]):
                    data.append(self.load_mongo_db(input_col, doc_type))
                counts = [cursor.count() for cursor in data]
            except TypeError:
                self.logger.exception(f'Failed to load the data in collection {self.mongo_params.input_cols} '
                                      f'on {self.mongo_db}')
            else:
                self.logger.info(f'Successfully loaded {counts} data points from collection '
                                 f'{self.mongo_params.input_cols} on {self.mongo_db}')

        if self.input_flags.sql_flag:
            try:
                self.load_sql_db()
            except TypeError:
                self.logger.exception(f'Failed to load the data in SQL database')
            else:
                self.logger.info(f'Successfully loaded data from SQL database')

        return data

    def load_mongo_db(self, input_col, doc_type):
        """
        Queries the collection in the database defined by self.mongo_db and retrieves the documents whose 'type' field
        equals taht of doc_type.

        Args:
            input_col (str): The name of the collection.
            doc_type (str): The value of the 'type' field of the documents to retrieve.

        Returns:
            pymongo cursor: The cursor containing the data.
        """

        return self.mongo_db[input_col].find({'type': doc_type})

    def load_sql_db(self):
        """
        Loads the data in the SQL database.
        """

        raise OSError

    @staticmethod
    def merge(mol1, mol2, map_num1, map_num2, stereo=None, clear_map_nums=True):
        """
        Static method for merging two molecules at the specified atoms, and updating the hydrogen counts as needed. Atom
        map numbers should not be the same for both molecules.

        Args:
            mol1 (rdkit Mol): The molecule to be combined with mol2.
            mol2 (rdkit Mol): The molecule to be combined with mol1.
            map_num1 (int): The atom map number of the atom on mol1 that will form a bond with mol2.
            map_num2 (int): The atom map number of the atom on mol2 that will form a bond with mol1.

        Returns:
            rdkit Mol: The resulting molecule from combining mol1 and mol2 at the specified atoms.
        """

        def update_hydrogen_counts(atom, clear_map_nums):
            """
            Inner method for clearing the atom map number and updating hydrogen counts.

            Args:
                atom (rdkit Atom): The atom from a molecule that is going to form a bond with an atom from another
                    molecule.

            Returns:
                rdkit Atom: The reformatted atom.
            """

            if clear_map_nums:
                atom.SetAtomMapNum(0)

            if atom.GetSymbol() in ['N', 'O', 'S']:
                atom.SetNumExplicitHs(0)
            elif atom.GetSymbol() == 'C' and atom.GetNumExplicitHs() != 0:
                atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)

            return atom

        # find atoms that will form a bond together and update hydrogen counts
        combo = Chem.RWMol(Chem.CombineMols(mol1, mol2))
        Chem.SanitizeMol(combo)
        for atom in combo.GetAtoms():
            if atom.GetAtomMapNum() == map_num1:
                mol1_atom = update_hydrogen_counts(atom, clear_map_nums)
            elif atom.GetAtomMapNum() == map_num2:
                mol2_atom = update_hydrogen_counts(atom, clear_map_nums)

        # create bond, remove hydrogens, and sanitize
        combo.AddBond(mol1_atom.GetIdx(), mol2_atom.GetIdx(), order=Chem.rdchem.BondType.SINGLE)
        Chem.rdmolops.RemoveHs(combo)
        Chem.SanitizeMol(combo)

        # add stereochemistry as specified
        stereo_center = mol1_atom if mol1_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and \
            mol1_atom.GetTotalNumHs() != 2 else mol2_atom
        if stereo == 'CCW':
            stereo_center.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif stereo == 'CW':
            stereo_center.SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)

        return combo

########################################################################################################################
########################################################################################################################
########################################################################################################################


class DataInitializer(Base):
    """
    A class for converting data stored in chemdraw files (.sdf) into the proper formats to be operated on. The data is
    then stored in the specified file or database. The data stored in chemdraw files should be the starting point for
    generating macrocyles, such as the templates, parent side chains, backbones, and monomers that won't be derived from
    combinations of side chains and backbones such as natural amino acids and modified prolines. Inherits from Base.
    """

    def __init__(self, f_in, f_out, logger=LOGGER, output_flags=Flags(True, False, True, False)):
        """
        Constructor.

        Args:
            f_in (str): The input filepath relative to DI_INPUT_DIR. Should only be a single file.
            f_out (str): The output filepath relative to DATA_DIR. Should only be a single file.
            logger (Logger, optional): The logger. Defaults to LOGGER.
            output_flags (Flags, optional): A namedtuple containing the flags indicating which format to output data to.
                Defaults to Flags(True, False, True, False).
        """

        # conver to list
        if not isinstance(f_in, list):
            f_in = [f_in]

        # I/O
        f_in = [os.path.join(DI_INPUT_DIR, file) for file in f_in]
        input_flags = Flags(False, False, False, False)
        super().__init__(IOPaths(f_in, f_out), MongoParams(None, None, COLLECTIONS.mols), logger, input_flags,
                         output_flags)

    def load_files(self, template_doc):
        """
        Helper function that reads in all molecules in the sdf file defined by self.fp_in, stores the mol's SMILES
        string and kekule SMILES string into the template document, and stores the result in self.result_data.

        Args:
            template_doc (OrderedDict): The dictionary that will contain the SMILES strings.

        Returns:
            bool: True if successful
        """

        try:
            for mol in read_mols(self.fp_in[0]):
                doc = deepcopy(template_doc)
                Chem.Kekulize(mol)
                doc['smiles'] = Chem.MolToSmiles(mol)
                doc['kekule'] = Chem.MolToSmiles(mol, kekuleSmiles=True)
                self.result_data.append(doc)
        except (OSError, Exception):
            self.result_data = []
            self.logger.exception(
                f'Variables: self.fp_in = {self.fp_in}, mol = {Chem.MolToSmiles(mol)}')
            return False

        return True

    def load_parent_side_chains(self, group=None):
        """
        Top level function that sets up a template document for parent side chains to be stored in, and passes it to
        self.load_files() for filling. If call to self.load_files() is successful then makes call to self.save_data().

        Args:
            group (str, optional): The group name assigned to each molecule in file. Defaults to the last word in the
                file name when split by '_' character.

        Returns:
            bool: True if successful.
        """

        group = self.fp_in[0].split('_')[-1].split('.')[0] if group is None else group
        template_doc = OrderedDict([('ID', None), ('type', 'parent_side_chain'),
                                    ('smiles', None), ('kekule', None), ('group', group)])
        if self.load_files(template_doc):
            return self.save_data(create_id=True)

        return False

    def load_monomers(self, group=None, required=False, backbone=None):
        """
        Top level function that sets up template document for monomers (that are not derived from modified side chains
        and backbones) to be stored in and passes it to self.load_files() for filling. If call to self.load_files() is
        successful then makes call to self.save_data().

        Args:
            group (str, optional): The group name assigned to each molecule in file. Defaults to the last word in the
                file name when split by '_' character.
            required (bool, optional): The value of 'required' stored in the documents. This value determines if the
                peptides these monomers are present in are valid if no other required monomers are present. Defaults to
                False.
            backbone (str, optional): The backbone type these monomers are made of. Defaults to None.

        Returns:
            bool: True if successful.
        """

        group = self.fp_in[0].split('_')[-1].split('.')[0] if group is None else group
        template_doc = OrderedDict([('ID', 'm'), ('type', 'monomer'), ('smiles', None), ('kekule', None),
                                    ('backbone', backbone), ('side_chain', None), ('group', group), ('required', required)])
        if self.load_files(template_doc):
            return self.save_data(create_id=True)

        return False

    def load_templates(self, identifier=['t1', 't2', 't3']):
        """
        Top level function that sets up template documents for the template molecules to be stored in and passes it to
        self.load_files() for filling. If call to self.load_files() is successful then makes call to self.save_data().
        Passed in identifier list must be in same order as molecules are read from the sdf file. Defaults indentifiers
        are as follows:
            - 't1' = template1
            - 't2' = template2
            - 't3' = template3

        Args:
            ID (list, optional): A list of strings to be used as IDs. Defaults to ['t1', 't2', 't3'].

        Returns:
            bool: True if successful.
        """

        template_doc = OrderedDict([('ID', None), ('type', 'template'), ('smiles', None), ('kekule', None)])
        if self.load_files(template_doc):
            for doc, id_val in zip(self.result_data, identifier):
                doc['ID'] = id_val
            return self.save_data(create_id=False)

        return False

    def set_template_atom_maps(self, doc):
        pass

    def load_backbones(self, identifier=['a', 'b', 'c']):
        """
        Top level function that sets up template documents for the backbone molecules to be stored in and passes it to
        self.load_files() for filling. If call to self.load_files() is successful then makes call to self.save_data().
        Passed in identifier list must be in same order as molecules are read from the sdf file. Default identifiers are
        as follows:
            - 'a' = 'alpha amino acid'
            - 'b' = 'beta2 amino acid'
            - 'c' = 'beta3 amino acid'

        Args:
            identifier (list, optional): A list of strings to be used as IDs. Defaults to ['a', 'b', 'c'].

        Returns:
            bool: True if successful.
        """

        template_doc = OrderedDict([('ID', None), ('type', 'backbone'), ('smiles', None), ('kekule', None)])
        if self.load_files(template_doc):
            for doc, id_val in zip(self.result_data, identifier):
                doc['ID'] = id_val
                doc = self.set_backbone_atom_maps(doc)
            return self.save_data(create_id=False)

        return False

    def set_backbone_atom_maps(self, doc):
        """
        Method for setting the atom map numbers on the backbone molecules. Atom map numbers are used to determine where
        on the backbone to attach the side chain.

        Args:
            doc (dict): The dictionary containing the associated backbone data.

        Returns:
            dict: The dictionary containing the atom mapped backbone SMILES string.
        """

        mol = Chem.MolFromSmiles(doc['smiles'])
        for atom in mol.GetAtoms():
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            if doc['ID'] in ('a', 'b') and (neighbors in (['C', 'N'], ['N', 'C'])):
                atom.SetAtomMapNum(BB_MAP_NUM)
                doc['smiles'] = Chem.MolToSmiles(mol)
                break
            elif doc['ID'] == 'c' and neighbors == ['C', 'C']:
                atom.SetAtomMapNum(BB_MAP_NUM)
                doc['smiles'] = Chem.MolToSmiles(mol)
                break

        return doc

    def save_data(self, create_id=False):
        """
        Saves data stored in self.result_data to the specified formats/locations defined by self.output_flags.

        Args:
            create_id (bool, optional): Determines whether to create new ID values for the data. Should only be used
                when inserting data into the database that is not derived. Defaults to False.

        Returns:
            bool: True if successful.
        """

        # save to Mongo database first to create new IDs
        try:
            self.save_mongo_db(create_id=create_id)
        except (DuplicateKeyError, ValueError, TypeError):
            self.logger.exception('Failed to save result_data to the Mongo database')
            return False
        except BulkWriteError as err:
            self.logger.exception(f'Failed to save result_data to the Mongo database\n{err.details}')
            return False
        else:
            self.logger.info(f'Successfully saved {len(self.result_data)} data points to the Mongo database!')

            # get doc_type and reset flags
            doc_type = self.result_data[0]['type']
            flag = self.output_flags
            self.output_flags = Flags(flag.json_flag, flag.txt_flag, False, flag.sql_flag)

            # load data just written with assigned IDs, and write to remaining formats
            self.result_data = list(self.load_mongo_db(COLLECTIONS.mols, doc_type))
            return super().save_data()

########################################################################################################################
########################################################################################################################
########################################################################################################################


def load_json(filepath):
    """
    Loads the data defined by the file in filepath and returns a generator object containing the data.

    Args:
        filepath (str): The filepath to the json file.

    Yields:
        (dict): A single data point.
    """

    with open(filepath, 'r') as f:
        for doc in json_util.loads(json_util.dumps(json.load(f))):  # need to convert _id back to ObjectID()
            yield doc


def load_txt(filepath):
    """
    Loads the data defined by the file in filepath and returns a generator object containing the data.

    Args:
        filepath (str): The filepath to the txt file.

    Yields:
        (str): A single data point.
    """

    with open(filepath, 'r') as f:
        while True:
            line = f.readline().rstrip('\n')
            if not line:
                break

            yield line


def read_mols(filepath=None, verbose=False):
    """
    Reads in the molecules defined by the sdf file in filepath.

    Args:
        filepath (str, optional): The filepath to the sdf file to be read. Defaults to None.
        verbose (bool, optional): If True, prints the molecules' SMILES strings to the console. Defaults to False.

    Returns:
        iterable: An iterable containing rdkit Mols.
    """

    # set default
    if filepath is None:
        filepath = os.path.join(DATA_DIR, 'chemdraw', 'test_rxn.sdf')

    mols = Chem.SDMolSupplier(filepath)

    # print smiles
    if verbose:
        for mol in mols:
            print(Chem.MolToSmiles(mol))

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

    fp = '/Users/ericdang/Desktop/' + file
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

    # write monomers to file
    fp = str(Path(__file__).resolve().parents[1] / fp)
    with open(fp, 'w') as f:
        json.dump(collection, f)


def ranges(total, chunks):
    step = total / chunks
    return [(round(step*i), round(step*(i+1))) for i in range(chunks)]


def kekulize_smiles(mols, smiles):

    if isinstance(smiles, str):
        smiles = [smiles]

    mols = [Chem.MolToSmiles(mol) for mol in mols]
    smiles.extend(mols)

    mols = [Chem.MolFromSmiles(mol) for mol in smiles]
    [Chem.Kekulize(mol) for mol in mols]
    [print(Chem.MolToSmiles(mol, kekuleSmiles=True)) for mol in mols]
    return True


def conformers_to_pdb(candidates=None):

    if candidates is None:
        cursor = Database(db='molecules').find_all('conformers')
    else:
        cursor = Database(db='molecules').find('conformers', {'candidate': candidates, 'num_conformers': 500})
    # else:
        # if isinstance(candidates, str):
        #     candidates = [candidates]
        #     for candidate in candidates:
        #     confs = Database(db='molecules').find({'candidate': candidate}, {'conformers': 1})

    for ind, mols in enumerate(cursor):
        mol = Chem.Mol(mols['binary'])
        # for conf in mols['conformers']:
        #     print([re.split(r'Energy: \d+\.\d+', conf)[1]])
        #     mol = Chem.MolFromMolBlock(re.split(r'Energy: \d+\.\d+\s', conf)[1])
        #     print(Chem.MolToSmiles(mol))
        #     exit()
        # [mol.AddConformer(Chem.MolFromMolBlock(conf)) for conf in mols['conformers']]
        rdmolfiles.MolToPDBFile(mol, '500.pdb')

    return True


def set_flags(inputs, outputs):
    input_flags = [False, False, False, False]
    output_flags = [False, False, False, False]

    # set input flags
    for i, val in enumerate(['json', 'txt', 'mongo', 'sql']):
        if val in inputs:
            input_flags[i] = True

        if val in outputs:
            output_flags[i] = True

    input_flags = Flags(input_flags[0], input_flags[1], input_flags[2], input_flags[3])
    output_flags = Flags(output_flags[0], output_flags[1], output_flags[2], output_flags[3])
    return input_flags, output_flags


def test_valid_smiles(smiles):
    """
    Test for a valid SMILES string. Rdkit does not provide a way to catch the error associated with trying to convert an
    invalid SMILES string to a Mol in python, thus this function is necessary to do so.

    Args:
        smiles (str): The SMILES string to be tested.

    Raises:
        InvalidSmilesString: The exception that is raised if unable to convert the SMILES string to a Mol.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSmilesString(f'Could not convert SMILES string \'{smiles}\' to rdkit Mol')
