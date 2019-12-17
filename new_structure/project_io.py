import csv
import json
import os
from abc import ABC, abstractmethod

from bson import json_util
from pymongo import MongoClient, errors
from rdkit import Chem

import config
import utils


class IOInterface(ABC):
    """
    Interface for classes that handle IO of different molecule and reaction data.
    """

    @abstractmethod
    def load(self):
        """
        Abstract method for loading data handled by the specific derived class.
        """

    @abstractmethod
    def save(self, data):
        """
        Abstract method for saving the provided data in the location maintained by the derived class.

        Args:
            data (iterable[dict]): The data to be saved.
        """


class RawRegioSQMIO(IOInterface):

    _RESULT_FILEPATH = os.path.join(config.DATA_DIR, 'external', 'regiosqm_results_nm_3.csv')
    _SMILES_FILEPATH = os.path.join(config.DATA_DIR, 'external', 'regiosqm_sidechains.smiles')

    def load(self):

        with open(self._RESULT_FILEPATH, 'r') as file:
            return csv.reader(file, delimiter=',')

    def save(self, data):

        with open(self._SMILES_FILEPATH, 'w') as file:
            for sidechain in data:
                file.write(sidechain['_id'] + ' ' + sidechain['kekule'] + '\n')


class ChemDrawIO(IOInterface):

    def __init__(self, filepath):

        self.filepath = filepath

    def load(self):

        return Chem.SDMolSupplier(self.filepath)

    def save(self, data):

        writer = Chem.SDWriter(self.filepath)
        writer.write(data)
        writer.close()


class ChemDrawSideChainIO(ChemDrawIO):

    _FILEPATH = os.path.join(config.DATA_DIR, 'chemdraw', 'sidechains.sdf')

    def __init__(self):
        super().__init__(self._FILEPATH)


class ChemDrawMonomerIO(ChemDrawIO):

    _FILEPATH = os.path.join(config.DATA_DIR, 'chemdraw', 'monomers.sdf')

    def __init__(self):
        super().__init__(self._FILEPATH)


class AbstractJsonIO(IOInterface):
    """
    Abstract class for IO classes that read and write to json file format.
    """

    def __init__(self, file_num_range=(None, None)):
        """
        Intializer that sets instance variables 'low' and 'high' that determine the range of file numbers to load.

        Args:
            file_num_range (tuple[int]): A tuple containing two intergers specifying the range of file numbers.

        Raises:
            ValueError: Raised if file_num_range doesn't have at least two elements.
        """

        try:
            self.low, self.high = file_num_range
        except ValueError:
            raise ValueError('file_num_range must be a tuple containing two integers')

    def from_json(self, filepath):
        """
        Loads the data in the json file pointed to by filepath. The file name specified in filepath is treated as a
        base file name, meaning there could be a set of files within the directory with the same base file name
        distinquished only by a number before the extension. For example:

        base file name - test.json
        actual file names - test_0.json, test_1.json, ...

        The instance variables 'low' and 'high' dictate which files containing the base file name should be loaded. If
        'low' is not provided, this method will load data from all the numbered files in the directory with the provided
        base file name. If 'high' is not provided then it will load the data in the base file name marked with
        the specified low number. If both 'low' and 'high' are provided then it will load all data in the files with the
        base file name that are distinquished by the numbers within the range of low to high.

        Args:
            filepath (str): The path to the files containing the data to be loaded. The specific file name specified
                here should be the base file name.

        Yields:
            dict: A dictionary containing an entry of data in the json file.
        """

        # determine range based on the provided low and high arguments
        if self.low is None:
            self.low, self.high = utils.get_file_num_range(filepath)
        elif self.high is None:
            self.high = self.low + 1

        # yield all data from within the specified range
        for file_num in range(self.low, self.high):
            with open(utils.attach_file_num(filepath, file_num), 'r') as file:
                for doc in json_util.loads(json_util.dumps(json.load(file))):
                    yield doc

    def to_json(self, filepath, data):
        """
        Saves the provided data into a json file, where the file name is the base file name specified in the filepath,
        appended with a number before the extension so as to ensure uniques. See from_json() doc string for example.

        Args:
            filepath (str): The path to the file where the data is to be saved. The file name specified here is treated
                as a base file name rather than the absolute file name.
            data (iterable[dict]): The data to be saved.
        """

        with open(utils.file_rotator(filepath), 'w') as file:
            json.dump(json.loads(json_util.dumps(data)), file)


class JsonIDIO(AbstractJsonIO):
    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'ids.json')

    def load(self):
        with open(self._FILEPATH, 'r') as file:
            return json.load(file)

    def save(self, data):
        with open(self._FILEPATH, 'w', encoding='utf-8') as file:
            json.dump(data, file)


class JsonIndexIO(AbstractJsonIO):
    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'index.json')

    def load(self):
        with open(self._FILEPATH, 'r') as file:
            return json.load(file)

    def save(self, data):
        with open(self._FILEPATH, 'w', encoding='utf-8') as file:
            json.dump(data, file)


class JsonSideChainIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling sidechain data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'sidechains.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonMonomerIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling monomer data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'monomers.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonPeptideIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling peptide data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'peptides.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonTemplatePeptideIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling template_peptide data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'template_peptides.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonMacrocycleIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling macrocycle data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'macrocycles.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonReactionIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling reaction data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'reactions.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonRegioSQMIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling RegioSQM data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'regiosqm.json')

    def load(self):

        return super().from_json(self._FILEPATH)

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class MongoDataBase():
    """
    Class for setting up connection and accessing the MongoDataBase as well as intializing it.
    """

    def __init__(self):
        """
        Initializer that establishes a connection to the MongoDataBase based on the host, port, and database name
        specified in the config file.
        """

        self.client = MongoClient(config.HOST, config.PORT)
        self.database = self.client[config.DATABASE]

    def __del__(self):
        """
        Destructor that closes the connection to the database during garbage collection.
        """

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

        return self.database[collection]

    def setup(self, validation_level='strict', clear=True):
        """
        Set up the database scheme specified in the MongoDataBase section of config.py.

        Args:
            validation_level (str, optional): Set the validation level of each collection. Defaults to 'strict'.
            clear (bool, optional): If True, clears collections in the database. Defaults to True.
        """

        # clear existing collection in database
        if clear:
            self.clear()

        # create collections, add validators and indexes
        for collection, validator, index in zip(config.COLLECTIONS, config.VALIDATORS, config.INDICES):
            query = {'collMod': collection,
                     'validator': validator,
                     'validationLevel': validation_level}
            self.database.create_collection(collection)
            self.database.command(query)

            if index is not None:
                for ind in index:
                    self[collection].create_index(ind, unique=True)

    def clear(self):
        """
        Helper method that removes all collections from the database.
        """

        for collection in self.database.list_collection_names():
            self[collection].drop()


class AbstractMongoIO(IOInterface, MongoDataBase):
    """
    Abstract class for IO classes that read and write to a MongoDataBase.
    """

    def from_mongo(self, collection, query):
        """
        Loads the data specified by the query in the given collection.

        Args:
            collection (str): The collection in the MongoDataBase to load from.
            query (dict): The query that specifies what data to load.

        Returns:
            list: The requested data.
        """

        return list(self[collection].find(query))

    def to_mongo(self, collection, data):
        """
        Saves the data to the specified collection in the MongoDataBase.

        Args:
            collection (str): The name of the collection in the MongoDataBase to save to.
            data (iterable[dict]): The data to be saved.
        """

        # try:
        self[collection].insert_many(data, ordered=True)
        # except errors.BulkWriteError as err:
        #     print(err.details)
        #     raise


class MongoIDIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling id data.
    """

    _COLLECTION = config.COL4
    _QUERY = {'_id': 'id'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)[0]

    def save(self, data):

        try:
            super().to_mongo(self._COLLECTION, [data])
        except errors.BulkWriteError:
            self.update(data)

    def update(self, data):

        self[self._COLLECTION].find_one_and_replace(self._QUERY, data)


class MongoIndexIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling index data.
    """

    _COLLECTION = config.COL4
    _QUERY = {'_id': 'index'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)[0]

    def save(self, data):

        try:
            super().to_mongo(self._COLLECTION, [data])
        except errors.BulkWriteError:
            self.update(data)

    def update(self, data):
        self[self._COLLECTION].find_one_and_replace(self._QUERY, data)


class MongoSideChainIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling sidechain data.
    """

    _COLLECTION = config.COL1
    _QUERY = {'type': 'sidechain'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)


class MongoMonomerIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling monomer data.
    """

    _COLLECTION = config.COL1
    _QUERY = {'type': 'monomer'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)


class MongoPeptideIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling peptide data.
    """

    _COLLECTION = config.COL1
    _QUERY = {'type': 'peptide'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)


class MongoTemplatePeptideIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling template_peptide data.
    """

    _COLLECTION = config.COL1
    _QUERY = {'type': 'template_peptide'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)


class MongoMacrocycleIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling macrocycle data.
    """

    _COLLECTION = config.COL1
    _QUERY = {'type': 'macrocycle'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)


class MongoReactionIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling reaction data.
    """

    _COLLECTION = config.COL2
    _QUERY = {}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)


class MongoRegioSQMIO(AbstractMongoIO):

    _COLLECTION = config.COL3
    _QUERY = {'type': 'regiosqm'}

    def load(self):

        return super().from_mongo(self._COLLECTION, self._QUERY)

    def save(self, data):

        super().to_mongo(self._COLLECTION, data)
