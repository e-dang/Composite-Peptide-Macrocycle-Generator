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
    Interface for classes that handle IO of different record data (molecules, reactions, ect...).
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
    """
    Class for handling the IO of the PeptidePlanner data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'peptide_plan.txt')

    def __init__(self, peptide_length):
        """
        Initializer.

        Args:
            peptide_length (int): The length of the peptides in the peptide plans.
        """

        self.peptide_length = peptide_length

    def load(self):
        """
        Method for loading the peptide plans from file.

        Returns:
            iterable: The lines of peptide plan file.
        """

        with open(utils.attach_file_num(self.FILEPATH, self.peptide_length), 'r') as file:
            return file.readlines()

    def save(self, data):
        """
        Method for saving the peptide plans to file.

        Args:
            data (iterable[tuple[int]]): The peptide plan data.
        """

        with open(utils.attach_file_num(self.FILEPATH, self.peptide_length), 'w') as file:
            for monomer_idxs in data:
                monomer_idxs = str(monomer_idxs).strip('(').strip(')').replace(' ', '')
                file.write(monomer_idxs + '\n')


class RawRegioSQMIO(IOInterface):
    """
    Class for saving the SMILES file needed as input for RegioSQM and loading of the output of RegioSQM.
    """

    RESULT_FILEPATH = os.path.join(config.DATA_DIR, 'external', 'regiosqm_results_nm_3.csv')
    SMILES_FILEPATH = os.path.join(config.DATA_DIR, 'external', 'regiosqm_mols.smiles')

    def load(self):
        """
        Method for loading the output of RegioSQM.

        Returns:
            list: The unformatted RegioSQM predicitions.
        """

        with open(self.RESULT_FILEPATH, 'r') as file:
            return list(csv.reader(file, delimiter=','))

    def save(self, data):
        """
        Method for saving the SMILES file needed as input for RegioSQM.

        Args:
            data (iterable[dict]): The SMILES data.
        """

        with open(self.SMILES_FILEPATH, 'w') as file:
            for mol in data:
                file.write(mol['_id'] + ' ' + mol['kekule'] + '\n')


class RawpKaIO(IOInterface):
    """
    Class for saving the SMILES file for tabulating pKas and loading the raw pKa predictions.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'external', 'pkas.txt')

    def load(self):
        """
        Method for loading the pKa predictions from txt file.

        Returns:
            list: The lines of the txt file containing the pKa predictions.
        """

        with open(self.FILEPATH, 'r') as file:
            return list(file.readlines())

    def save(self, data):
        """
        Method for saving the SMILES file for tabulating pKa predictions.

        Args:
            data (iterable[tuple[str]]): The molecules' id and kekule SMILES as tuples.
        """

        with open(self.FILEPATH, 'w') as file:
            for mol_id, kekule in data:
                file.write(mol_id + '; ' + kekule + '; ' + '\n')


class StructureDataFileIO(IOInterface):
    """
    Abstract class for IO classes that read and write to SDF files.
    """

    def __init__(self, filepath):
        """
        Initializer.

        Args:
            filepath (str): The filepath to an sdf file.
        """

        self.filepath = filepath

    def load(self):
        """
        Method for loading the molecules from the sdf file specified by self.filepath.

        Returns:
            RDKit MolSupplier: The molecules in the sdf file.
        """

        return Chem.SDMolSupplier(self.filepath)

    def save(self, data):
        """
        Method for saving molecules to an sdf file specified by self.filepath.

        Args:
            data (iterable[RDKit Mol]): The molecules to write to file.
        """

        writer = Chem.SDWriter(self.filepath)
        for mol in data:
            writer.write(mol)
        writer.close()


class SDFSideChainIO(StructureDataFileIO):
    """
    Class for handling the IO of sidechain data in sdf format.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'chemdraw', 'sidechains.sdf')

    def __init__(self):
        """
        Initializer that sets parent classes' member variable self.filepath to the class variable self.FILEPATH.
        """

        super().__init__(self.FILEPATH)


class SDFMonomerIO(StructureDataFileIO):
    """
    Class for handling the IO of monomer data in sdf format.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'chemdraw', 'monomers.sdf')

    def __init__(self):
        """
        Initializer that sets parent classes' member variable self.filepath to the class variable self.FILEPATH.
        """

        super().__init__(self.FILEPATH)


class AbstractJsonIO(IOInterface):
    """
    Abstract class for IO classes that read and write to json file format.
    """

    def __init__(self, filepath, file_num_range=(None, None)):
        """
        Intializer.

        Args:
            filepath (str): The filepath to a json file.
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
            data (iterable[dict]): The data to be saved.
        """

        with open(utils.file_rotator(self.filepath), 'w') as file:
            json.dump(json.loads(json_util.dumps(data)), file)


class JsonIDIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling id record data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'ids.json')

    def __init__(self, **kwargs):
        """
        Initializer.
        """

    def load(self):
        """
        Method for loading the id record data from file.

        Returns:
            iterable[dict]: The id record data.
        """

        with open(self.FILEPATH, 'r') as file:
            return json.load(file)

    def save(self, data):
        """
        Method for saving id record data to file.

        Args:
            data ([iterable[dict]]): The id record data.
        """

        with open(self.FILEPATH, 'w', encoding='utf-8') as file:
            json.dump(data, file)


class JsonIndexIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling index record data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'index.json')

    def __init__(self, **kwargs):
        """
        Initializer.
        """

    def load(self):
        """
        Method for loading the index record data from file.

        Returns:
            iterable[dict]: The index record data.
        """

        with open(self.FILEPATH, 'r') as file:
            return json.load(file)

    def save(self, data):
        """
        Method for saving the index record data to file.

        Args:
            data (iterable[dict]): The index record data.
        """

        with open(self.FILEPATH, 'w', encoding='utf-8') as file:
            json.dump(data, file)


class JsonSideChainIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling sidechain data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'sidechains.json')

    def __init__(self, file_num_range=(None, None)):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
        """

        super().__init__(self.FILEPATH, file_num_range)


class JsonMonomerIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling monomer data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'monomers.json')

    def __init__(self, file_num_range=(None, None)):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
        """

        super().__init__(self.FILEPATH, file_num_range)


class JsonPeptideIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling peptide data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'peptides.json')

    def __init__(self, file_num_range=(None, None), **kwargs):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
            peptide_length (int): The length of the peptides handled by this class.
        """

        super().__init__(utils.attach_file_num(self.FILEPATH, kwargs['peptide_length']), file_num_range)


class JsonTemplatePeptideIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling template_peptide data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'template_peptides.json')

    def __init__(self, file_num_range=(None, None), **kwargs):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
            peptide_length (int): The length of the peptides handled by this class.
        """

        super().__init__(utils.attach_file_num(self.FILEPATH, kwargs['peptide_length']), file_num_range)


class JsonMacrocycleIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling macrocycle data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'macrocycles.json')

    def __init__(self, file_num_range=(None, None), **kwargs):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
            peptide_length (int): The length of the peptides handled by this class.
            job_num (int): The job number assigned to the current process (if not ran in job array then set to 1).
        """

        super().__init__(utils.attach_file_num(self.FILEPATH,
                                               kwargs['job_num'], kwargs['peptide_length']), file_num_range)


class JsonReactionIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling reaction data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'reactions.json')

    def __init__(self, file_num_range=(None, None)):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
        """

        super().__init__(self.FILEPATH, file_num_range)


class JsonRegioSQMIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling formatted RegioSQM data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'regiosqm.json')

    def __init__(self, file_num_range=(None, None)):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
        """

        super().__init__(self.FILEPATH, file_num_range)


class JsonpKaIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling formatted pKa data.
    """

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'pka.json')

    def __init__(self, file_num_range=(None, None)):
        """
        Initializer.

        Args:
            file_num_range (tuple, optional): A tuple containing two intergers specifying the range of file numbers.
                See AbstractJsonIO's load() docstring for more information on what this argument does. Defaults to
                (None, None).
        """

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

    def setup(self, validation_level='strict'):
        """
        Set up the database scheme specified in the MongoDataBase section of config.py.

        Args:
            validation_level (str, optional): Set the validation level of each collection. Defaults to 'strict'.
        """

        # clear existing collection in database
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
        """
        Initializer.

        Args:
            collection (str): The collection within the database to access.
            query (dict): The query to use to load the data associated with the derived class.
        """

        self.collection = collection
        self.query = query

    def load(self):
        """
        Loads the data specified by the query in the given collection.

        Returns:
            list: The requested data.
        """

        return list(self[self.collection].find(self.query))

    def save(self, data):
        """
        Saves the data to the specified collection in the MongoDataBase.

        Args:
            data (iterable[dict]): The data to be saved.
        """

        try:
            self[self.collection].insert_many(data, ordered=True)
        except errors.BulkWriteError as err:
            print(err.details)
            raise


class MongoIDIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling id record data.
    """

    COLLECTION = config.COL4
    QUERY = {'_id': 'id'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)

    def load(self):
        """
        Method for loading the id record data from the MongoDataBase.

        Returns:
            dict: The id record data
        """

        return super().load()[0]

    def save(self, data):
        """
        Method for saving the id record data to the MongoDataBase.

        Args:
            data (dict): The id record data.
        """

        try:
            super().save([data])
        except errors.BulkWriteError:
            self.update(data)

    def update(self, data):
        """
        Helper method for updating the the id record data if it already exists in the MongoDataBase.

        Args:
            data (dict): The id record data.
        """

        self[self.COLLECTION].find_one_and_replace(self.QUERY, data)


class MongoIndexIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling index record data.
    """

    COLLECTION = config.COL4
    QUERY = {'_id': 'index'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)

    def load(self):
        """
        Method for loading the index record data from the MongoDataBase.

        Returns:
            dict: The index record data
        """

        return super().load()[0]

    def save(self, data):
        """
        Method for saving the index record data to the MongoDataBase.

        Args:
            data (dict): The index record data.
        """

        try:
            super().save([data])
        except errors.BulkWriteError:
            self.update(data)

    def update(self, data):
        """
        Helper method for updating the the index record data if it already exists in the MongoDataBase.

        Args:
            data (dict): The index record data.
        """

        self[self.COLLECTION].find_one_and_replace(self.QUERY, data)


class MongoSideChainIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling sidechain data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'sidechain'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


class MongoMonomerIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling monomer data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'monomer'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


class MongoPeptideIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling peptide data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'peptide'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


class MongoTemplatePeptideIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling template_peptide data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'template_peptide'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


class MongoMacrocycleIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling macrocycle data.
    """

    COLLECTION = config.COL1
    QUERY = {'type': 'macrocycle'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


class MongoReactionIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling reaction data.
    """

    COLLECTION = config.COL2
    QUERY = {}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


class MongoRegioSQMIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling RegioSQM prediction data.
    """

    COLLECTION = config.COL3
    QUERY = {'type': 'regiosqm'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


class MongopKaIO(AbstractMongoIO):
    """
    Implmentation of the AbstractMongoIO class for handling pKa prediction data.
    """

    COLLECTION = config.COL3
    QUERY = {'type': 'pka'}

    def __init__(self):
        """
        Initializer that calls base class initializer with the collection and query class constants.
        """

        super().__init__(self.COLLECTION, self.QUERY)


def get_io(json_io, mongo_io):
    """
    Closure for creating functions that returns the correct IO class depending on the data format specified in the
    config file.

    Args:
        json_io (AbstractJsonIO): An implementation of the AbstractJsonIO class.
        mongo_io (AbstractMongoIO): An implementation of the AbstractMongoIO class.

    Returns:
        func: A function that can be called in order to recieve the correct IO class.
    """

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
    """
    Closure for creating a function that will hash the provided predictions based on the reacting_mol attribute.

    Args:
        predictions (iterable[dict]): The predictions to hash.

    Returns:
        func: A function that can be called in order to recieve the hashed predictions.
    """

    def hasher():
        hashed_predictions = {}
        for prediction in predictions:
            hashed_predictions[prediction['reacting_mol']] = prediction['predictions']
        return hashed_predictions

    return hasher


get_hashed_regiosqm_predictions = get_hashed_predictions(get_regiosqm_io().load())
get_hashed_pka_predictions = get_hashed_predictions(get_pka_io().load())
