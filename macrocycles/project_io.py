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


class PeptidePlannerIO(IOInterface):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'peptide_plan.txt')

    def __init__(self, peptide_length):
        self.peptide_length = peptide_length

    def load(self):
        with open(utils.attach_file_num(self.FILEPATH, self.peptide_length), 'r') as file:
            return file.readlines()

    def save(self, data):

        with open(utils.attach_file_num(self.FILEPATH, self.peptide_length), 'w') as file:
            for monomer_idxs in data:
                monomer_idxs = str(monomer_idxs).strip('(').strip(')').replace(' ', '')
                file.write(monomer_idxs + '\n')


class RawRegioSQMIO(IOInterface):

    RESULT_FILEPATH = os.path.join(config.DATA_DIR, 'external', 'regiosqm_results_nm_3.csv')
    SMILES_FILEPATH = os.path.join(config.DATA_DIR, 'external', 'regiosqm_mols.smiles')

    def load(self):

        with open(self.RESULT_FILEPATH, 'r') as file:
            return list(csv.reader(file, delimiter=','))

    def save(self, data):

        with open(self.SMILES_FILEPATH, 'w') as file:
            for mol in data:
                file.write(mol['_id'] + ' ' + mol['kekule'] + '\n')


class RawpKaIO(IOInterface):

    FILEPATH = os.path.join(config.DATA_DIR, 'external', 'pkas.txt')

    def load(self):
        with open(self.FILEPATH, 'r') as file:
            return list(file.readlines())

    def save(self, data):

        with open(self.FILEPATH, 'w') as file:
            for mol_id, kekule in data:
                file.write(mol_id + '; ' + kekule + '; ' + '\n')


class ChemDrawIO(IOInterface):

    def __init__(self, filepath):

        self.filepath = filepath

    def load(self):

        return Chem.SDMolSupplier(self.filepath)

    def save(self, data):

        writer = Chem.SDWriter(self.filepath)
        for mol in data:
            writer.write(mol)
        writer.close()


class ChemDrawSideChainIO(ChemDrawIO):

    FILEPATH = os.path.join(config.DATA_DIR, 'chemdraw', 'sidechains.sdf')

    def __init__(self):
        super().__init__(self.FILEPATH)


class ChemDrawMonomerIO(ChemDrawIO):

    FILEPATH = os.path.join(config.DATA_DIR, 'chemdraw', 'monomers.sdf')

    def __init__(self):
        super().__init__(self.FILEPATH)


class AbstractJsonIO(IOInterface):
    """
    Abstract class for IO classes that read and write to json file format.
    """

    def __init__(self, filepath, file_num_range=(None, None)):
        """
        Intializer that sets instance variables 'low' and 'high' that determine the range of file numbers to load.

        Args:
            file_num_range (tuple[int]): A tuple containing two intergers specifying the range of file numbers.

        Raises:
            ValueError: Raised if file_num_range doesn't have at least two elements.
        """

        try:
            self.filepath = filepath
            self.low, self.high = file_num_range
        except ValueError:
            raise ValueError('file_num_range must be a tuple containing two integers')

    def load(self):
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
            self.low, self.high = utils.get_file_num_range(self.filepath)
        elif self.high is None:
            self.high = self.low + 1

        # yield all data from within the specified range
        for file_num in range(self.low, self.high):
            with open(utils.attach_file_num(self.filepath, file_num), 'r') as file:
                for doc in json_util.loads(json_util.dumps(json.load(file))):
                    yield doc

    def save(self, data):
        """
        Saves the provided data into a json file, where the file name is the base file name specified in the filepath,
        appended with a number before the extension so as to ensure uniques. See from_json() doc string for example.

        Args:
            filepath (str): The path to the file where the data is to be saved. The file name specified here is treated
                as a base file name rather than the absolute file name.
            data (iterable[dict]): The data to be saved.
        """

        with open(utils.file_rotator(self.filepath), 'w') as file:
            json.dump(json.loads(json_util.dumps(data)), file)


class JsonIDIO(AbstractJsonIO):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'ids.json')

    def __init__(self, **kwargs):
        pass

    def load(self):
        with open(self.FILEPATH, 'r') as file:
            return json.load(file)

    def save(self, data):
        with open(self.FILEPATH, 'w', encoding='utf-8') as file:
            json.dump(data, file)


class JsonIndexIO(AbstractJsonIO):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'index.json')

    def __init__(self, **kwargs):
        pass

    def load(self):
        with open(self.FILEPATH, 'r') as file:
            return json.load(file)

    def save(self, data):
        with open(self.FILEPATH, 'w', encoding='utf-8') as file:
            json.dump(data, file)


class JsonSideChainIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling sidechain data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'sidechains.json')

    def __init__(self, file_num_range=(None, None)):
        super().__init__(self.FILEPATH, file_num_range)


class JsonMonomerIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling monomer data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'monomers.json')

    def __init__(self, file_num_range=(None, None)):
        super().__init__(self.FILEPATH, file_num_range)


class JsonPeptideIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling peptide data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'peptides.json')

    def __init__(self, file_num_range=(None, None), **kwargs):
        super().__init__(utils.attach_file_num(self.FILEPATH, kwargs['peptide_length']), file_num_range)


class JsonTemplatePeptideIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling template_peptide data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'template_peptides.json')

    def __init__(self, file_num_range=(None, None), **kwargs):
        super().__init__(utils.attach_file_num(self.FILEPATH, kwargs['peptide_length']), file_num_range)


class JsonMacrocycleIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling macrocycle data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'macrocycles.json')

    def __init__(self, file_num_range=(None, None), **kwargs):
        super().__init__(utils.attach_file_num(self.FILEPATH,
                                               kwargs['job_num'], kwargs['peptide_length']), file_num_range)


class JsonReactionIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling reaction data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'reactions.json')

    def __init__(self, file_num_range=(None, None)):
        super().__init__(self.FILEPATH, file_num_range)


class JsonRegioSQMIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling RegioSQM data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'regiosqm.json')

    def __init__(self, file_num_range=(None, None)):
        super().__init__(self.FILEPATH, file_num_range)


class JsonpKaIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling pKa data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'pka.json')

    def __init__(self, file_num_range=(None, None)):
        super().__init__(self.FILEPATH, file_num_range)


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

    def __init__(self, collection, query):

        self.collection = collection
        self.query = query

    def load(self):
        """
        Loads the data specified by the query in the given collection.

        Args:
            collection (str): The collection in the MongoDataBase to load from.
            query (dict): The query that specifies what data to load.

        Returns:
            list: The requested data.
        """

        return list(self[self.collection].find(self.query))

    def save(self, data):
        """
        Saves the data to the specified collection in the MongoDataBase.

        Args:
            collection (str): The name of the collection in the MongoDataBase to save to.
            data (iterable[dict]): The data to be saved.
        """

        try:
            self[self.collection].insert_many(data, ordered=True)
        except errors.BulkWriteError as err:
            print(err.details)
            raise


class MongoIDIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling id data.
    """

    COLLECTION = config.COL4
    QUERY = {'_id': 'id'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)

    def load(self):

        return super().load()[0]

    def save(self, data):

        try:
            super().save([data])
        except errors.BulkWriteError:
            self.update(data)

    def update(self, data):

        self[self.COLLECTION].find_one_and_replace(self.QUERY, data)


class MongoIndexIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling index data.
    """

    COLLECTION = config.COL4
    QUERY = {'_id': 'index'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)

    def load(self):

        return super().load()[0]

    def save(self, data):

        try:
            super().save([data])
        except errors.BulkWriteError:
            self.update(data)

    def update(self, data):
        self[self.COLLECTION].find_one_and_replace(self.QUERY, data)


class MongoSideChainIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling sidechain data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'sidechain'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


class MongoMonomerIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling monomer data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'monomer'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


class MongoPeptideIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling peptide data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'peptide'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


class MongoTemplatePeptideIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling template_peptide data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'template_peptide'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


class MongoMacrocycleIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling macrocycle data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'macrocycle'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


class MongoReactionIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling reaction data.
    """

    COLLECTION = config.COL2
    QUERY = {}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


class MongoRegioSQMIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling RegioSQM prediction data.
    """

    COLLECTION = config.COL3
    QUERY = {'type': 'regiosqm'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


class MongopKaIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling pKa prediction data.
    """

    COLLECTION = config.COL3
    QUERY = {'type': 'pka'}

    def __init__(self):
        super().__init__(self.COLLECTION, self.QUERY)


def get_io(json_io, mongo_io):
    def io_closure(**kwargs):
        if config.DATA_FORMAT == 'json':
            return json_io(**kwargs)
        if config.DATA_FORMAT == 'mongo':
            return mongo_io()

    return io_closure


get_id_io = get_io(JsonIDIO, MongoIDIO)
get_index_io = get_io(JsonIndexIO, MongoIndexIO)
get_sidechain_io = get_io(JsonSideChainIO, MongoSideChainIO)
get_monomer_io = get_io(JsonMonomerIO, MongoMonomerIO)
get_peptide_io = get_io(JsonPeptideIO, MongoPeptideIO)
get_template_peptide_io = get_io(JsonTemplatePeptideIO, MongoTemplatePeptideIO)
get_macrocycle_io = get_io(JsonMacrocycleIO, MongoMacrocycleIO)
get_reaction_io = get_io(JsonReactionIO, MongoReactionIO)
get_regiosqm_io = get_io(JsonRegioSQMIO, MongoRegioSQMIO)
get_pka_io = get_io(JsonpKaIO, MongopKaIO)


def get_hashed_predictions(predictions):

    def hasher():
        hashed_predictions = {}
        for prediction in predictions:
            hashed_predictions[prediction['reacting_mol']] = prediction['predictions']
        return hashed_predictions

    return hasher


get_hashed_regiosqm_predictions = get_hashed_predictions(get_regiosqm_io().load())
get_hashed_pka_predictions = get_hashed_predictions(get_pka_io().load())
