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

from pymongo import ASCENDING, MongoClient
from pymongo.errors import (BulkWriteError, CollectionInvalid,
                            ConnectionFailure, DuplicateKeyError, InvalidName,
                            OperationFailure)
from rdkit import Chem
from rdkit.Chem import rdmolfiles

from macrocycles.config import (COLLECTIONS, DATA_DIR, MONGO_SETTINGS,
                                PROJECT_DIR, LOG_DIR, VALIDATORS, INDEX)
from macrocycles.exceptions import (SavingMongoError, SavingSQLError,
                                    WritingJsonError, WritingTxtError)

########################################################################################################################
########################################################################################################################
########################################################################################################################
Flags = namedtuple('Flags', ['json_flag', 'txt_flag', 'mongo_flag', 'sql_flag'])
Smiles = namedtuple('Smiles', ['smiles', 'kekule'])
IOPaths = namedtuple('IOPaths', ['input', 'output'])
Mongo = namedtuple('MongoAccess', ['collection', 'doc_type'])


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
        '[%(asctime)-19s] [%(levelname)-8s] [%(name)s] -- %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
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

        self.logger.info(f'Closing connection with {self.client} through garbage collection')
        self.client.close()

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
            doc['ID'], count, prefix = generate_id(count, prefix)

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
        fp_in: The filepath to the input file.
        fp_out: The filepath to the output file.
        mongo_db: A connection to the MongoDB where the result_data will be stored.
        collection: The collection in the MongoDB that result_data will be inserted into.
        doc_type: The value held within the field 'type' in the documents to be loaded.
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

    def __init__(self, paths, mongo_setup, logger, output_flags=Flags(False, False, False, False),
                 input_flags=Flags(False, False, True, False)):
        """
        Constructor.

        Args:
            paths (IOPaths): A namedtuple containing the I/O paths relative to DATA_DIR.
            mongo_setup (Mongo): A namedtuple containing the collection and doc_type names.
            logger (Logger): The logger to be used by this class.
            output_flags (Flags, optional): A namedtuple specifying output formats.
                Defaults to Flags(False, False, False, False).
            input_flags (Flags, optional): A namedtuple specifying the input format.
                Defaults to Flags(False, False, True, False).
        """
        # I/O
        self.project_dir = PROJECT_DIR
        self.data_dir = DATA_DIR
        self.fp_in = os.path.join(DATA_DIR, paths.input)
        self.fp_out = os.path.join(DATA_DIR, paths.output)
        self.mongo_db = MongoDataBase(logger=logger)
        self.collection, self.doc_type = mongo_setup.collection, mongo_setup.doc_type
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
            pass
        else:
            if True not in self.output_flags:
                self.logger.warning('No data saved. Please set a data flag specifying which data format to save to.')
            return True

        return False

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
                self.logger.exception('Failed to write result_data to a json file')
                raise WritingJsonError
            else:
                self.logger.info(
                    f'Successfully wrote {len(self.result_data)} data points to {self.fp_out.split("/")[-1]}.json!')

        if self.output_flags.txt_flag:
            try:
                self.write_txt()
            except TypeError:
                self.logger.exception('Failed to write result_data to a txt file')
                raise WritingTxtError
            else:
                self.logger.info(
                    f'Successfully wrote {len(self.result_data)} data points to {self.fp_out.split("/")[-1]}.txt!')

        if self.output_flags.mongo_flag:
            try:
                self.save_mongo_db()
            except (DuplicateKeyError, ValueError, TypeError):
                self.logger.exception('Failed to save result_data to the Mongo database')
                raise SavingMongoError
            except BulkWriteError as err:
                self.logger.exception(f'Failed to save result_data to the Mongo database\n{err.details}')
                raise SavingMongoError
            else:
                self.logger.info(f'Successfully saved {len(self.result_data)} data points to the Mongo database!')

        if self.output_flags.sql_flag:
            try:
                self.save_sql_db()
            except Exception:
                self.logger.exception('Failed to save result_data to the SQL database')
                raise SavingSQLError
            else:
                self.logger.info(f'Successfully saved {len(self.result_data)} data points to the SQL database!')

    def write_json(self):
        """
        Writes data stored in self.result_data to a json file specified by self.fp_out.
        """

        with open(self.fp_out + '.json', 'w') as f:
            json.dump(self.result_data, f)

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

    def save_mongo_db(self):
        """
        Saves data stored in self.result_data to the Mongo database specified by self.mongo_db.

        Returns:
            bool: True if successful.
        """

        return self.mongo_db.insert(self.collection, self.result_data)

    def save_sql_db(self):
        """
        Saves data stored in self.result_data to the SQL database specified by self.sql_db.

        Returns:
            bool: True if successful.
        """

        return self.sql_db.insert(self.result_data)

    def load_data(self):
        """
        Method to be overloaded by derived classes. Loads data from format specified by self.input_flags and returns
        the results as a generator.

        Returns:
            generator: Contains the data.
        """

        if sum(self.input_flags) > 1:
            self.logger.warning(f'No data loaded. Multiple input flags were set, please specify one. self.input_flags: '
                                f'{self.input_flags}')
        elif True not in self.input_flags:
            self.logger.warning(f'No data loaded. No input flags were set, please specify one. self.input_flags: '
                                f'{self.input_flags}')

        data = None
        if self.input_flags.json_flag:
            try:
                data = self.load_json()
            except (OSError, json.JSONDecodeError):
                self.logger.exception(f'Failed to load the json file(s) \'{self.fp_in}\'')
            else:
                self.logger.info(f'Successfully loaded the json file(s) \'{self.fp_in}\'')

        if self.input_flags.txt_flag:
            try:
                data = self.load_txt()
            except OSError:
                self.logger.exception(f'Failed to load the txt file(s) \'{self.fp_in}\'')
            else:
                self.logger.info(f'Successfully loaded the txt file(s) \'{self.fp_in}\'')

        if self.input_flags.mongo_flag:
            try:
                data = self.load_mongo()
                count = data.count()
            except TypeError:
                self.logger.exception(f'Failed to load the data in collection \'{self.collection}\' on {self.mongo_db}')
            else:
                self.logger.info(f'Successfully loaded {count} data points from collection \'{self.collection}\' on '
                                 f'{self.mongo_db}')

        if self.input_flags.sql_flag:
            try:
                self.load_sql()
            except TypeError:
                self.logger.exception(f'Failed to load the data in SQL database')
            else:
                self.logger.info(f'Successfully loaded data from SQL database')

        return data

    def load_json(self):
        """
        Loads the data defined by the file in self.fp_in and returns a generator object containing the data.

        Yields:
            (dict): A single data point.
        """

        with open(self.fp_in, 'r') as f:
            for doc in json.load(f):
                yield doc

    def load_txt(self):
        """
        Loads the data defined by the file in self.fp_in and returns a generator object containing the data.

        Yields:
            (str): A single data point.
        """

        with open(self.fp_in, 'r') as f:
            while True:
                line = f.readline().rstrip('\n')
                if not line:
                    break

                yield line

    def load_mongo(self):
        """
        Queries the collection defined by self.collection in the database defined by self.mongo_db for all documents
        of the type = self.doc_type.

        Returns:
            pymongo cursor: The cursor containing the data.
        """

        return self.mongo_db[self.collection].find({'type': self.doc_type})

    def load_sql(self):
        """
        Loads the data in the SQL database.
        """

        raise OSError

########################################################################################################################
########################################################################################################################
########################################################################################################################


class DataInitializer(Base):
    """
    A class for converting data stored in chemdraw files (.sdf) into the proper formats to be operated on. The data is
    then stored in the specified file or database. Inherits from Base.
    """

    def __init__(self, f_in, collection, doc_type, logger=LOGGER, mongo_flag=True, sql_flag=False):

        super().__init__(IOPaths(f_in, '/'), Mongo(collection, doc_type), logger,
                         Flags(False, False, mongo_flag, sql_flag))

    def load_files(self, group=None):
        """
        Reads in all molecules in the sdf file defined by self.fp_in, formats the data, and stores the results in
        self.result_data.

        Args:
            group (str, optional): The group name assigned to each molecule in file. Defaults to the last word in the
                file name when split by '_' character.

        Returns:
            bool: True if successful
        """

        group = self.fp_in.split('_')[-1].split('.')[0] if group is None else group

        try:
            for mol in read_mols(self.fp_in):
                Chem.Kekulize(mol)
                doc = OrderedDict([('ID', None), ('type', self.doc_type), ('smiles', Chem.MolToSmiles(mol)),
                                   ('kekule', Chem.MolToSmiles(mol, kekuleSmiles=True)), ('group', group)])
                self.result_data.append(doc)
        except (OSError, Exception):
            self.result_data = []
            self.logger.exception(f'Variables: self.fp_in = {self.fp_in}, group = {group}')
            return False

        return True

    def save_data(self, create_id=False):
        """
        Saves data stored in self.result_data to the specified formats/locations defined by self.output_flags.

        Args:
            create_id (bool, optional): Determines whether to create new ID values for the data. Should only be used
                when inserting data into the database that is not derived. Defaults to False.

        Raises:
            SavingMongoError: Failure to save data to the Mongo database.
            SavingSQLError: Failure to save data to the SQL database.
        """

        if self.output_flags.mongo_flag:
            try:
                self.mongo_db.insert(self.collection, self.result_data, create_id=create_id)
            except (DuplicateKeyError, ValueError):
                self.logger.exception('Failed to save result_data to the Mongo database')
                raise SavingMongoError
            except BulkWriteError as err:
                self.logger.exception(f'Failed to save result_data to the Mongo database\n{err.details}')
                raise SavingMongoError
            else:
                self.logger.info(f'Successfully saved {len(self.result_data)} data points to the Mongo database!')

        if self.output_flags.sql_flag:
            try:
                self.sql_db.insert()
            except Exception:
                self.logger.exception('Failed to save result_data to the SQL database')
                raise SavingSQLError
            else:
                self.logger.info(f'Successfully saved {len(self.result_data)} data points to the SQL database!')

########################################################################################################################
########################################################################################################################
########################################################################################################################


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


class Database():
    """
    A class to establish a connection a MongoDB database

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

        # check if connected to correct database
        if self.db.name != 'rxn_templates':
            print('Not connected to "rxn_templates" database.')
            return False

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

        # check if connected to correct database
        if self.db.name != 'rxn_templates':
            print('Not connected to "rxn_templates" database.')
            return False

        return self.db[collection].insert_one({'reaction_smarts': reaction_smarts, 'template': template,
                                               'side_chain': side_chain, 'atom_idx': atom_idx})

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

        # check if connected to correct database
        if self.db.name != 'molecules':
            print('Not connected to "molecules" database.')
            return False

        return self.db[collection].insert_one({'reactant': reactant, 'products': products, 'num_products': num_products,
                                               'peptide': peptide, 'template': template, 'monomers': monomers,
                                               'reacting_side_chains': reacting_side_chains, 'atom_idx': atom_idx})

    def insert_filtered_candidates(self, reactant, filter_type, reacting_side_chains, collection='filtered_candidates'):

        if self.db.name != 'molecules':
            print('Not connected to "molecules" database.')
            return False

        return self.db[collection].insert_one({'reactant': reactant, 'filter': filter_type,
                                               'reacting_side_chains': reacting_side_chains})

    def insert_conformers(self, candidate, binary, convergences, energies, rmsd, collection='conformers'):

        if self.db.name != 'molecules':
            print('Not connected to "molecules" database.')
            return False

        return self.db[collection].insert_one({'candidate': candidate, 'binary': binary,
                                               'num_conformers': len(convergences), 'convergences': convergences,
                                               'energies': energies, 'avg_rmsd': rmsd})
