import csv
import glob
import json
import os
from abc import ABC, abstractmethod
from collections import defaultdict
from itertools import chain

from bson import json_util
from pymongo import MongoClient, errors
from rdkit import Chem

import macrocycles.config as config
import macrocycles.utils as utils


class FileReader:
    def setup(self, filepath, sort_order, *args, delimiter='_', glob='*'):
        self.globber = FileGlober()
        self.sorter = FileSorter()
        self.filepaths = []
        filepaths, unique_file_props = self.globber(self.add_glob(
            utils.attach_file_num(filepath, *args), glob), sort_order, delimiter)
        if len(filepaths):
            self.filepaths = self.sorter(filepaths, unique_file_props, sort_order)
        self.data = []

    def __iter__(self):
        if len(self.filepaths):
            self.data_count = 0
            self.read_data(self.filepaths[0])
            self.fp_count = 1
        return self

    def __next__(self):
        if self.data_count < len(self.data):
            self.data_count += 1
            return self.data[self.data_count - 1]
        elif (self.data_count == len(self.data) or len(self.data) == 0) and self.fp_count < len(self.filepaths):
            self.read_data(self.filepaths[self.fp_count])
            self.fp_count += 1
            self.data_count = 1
            if len(self.data) > 0:
                return self.data[self.data_count - 1]
            else:
                return next(self)
        else:
            raise StopIteration

    def add_glob(self, filepath, glob):
        path, filename = os.path.split(filepath)
        base_filename, ext = os.path.splitext(filename)
        return os.path.join(path, base_filename + glob + ext)

    def read_data(self, filepath):
        pass


class JsonFileReader(FileReader):
    def __init__(self, filepath, sort_order, delimiter='_', glob='*'):
        super().__init__(filepath, sort_order, delimiter, glob)

    def read_data(self, filepath):
        with open(filepath, 'r') as file:
            self.data = list(json_util.loads(json_util.dumps(json.load(file))))


class FileGlober:
    def __call__(self, glob_filepath, properties, delimiter):
        filepaths = glob.glob(glob_filepath)
        unique_file_props = defaultdict(set)
        file_props = []

        for filepath in filepaths:
            base_filename, *prop_vals = self.parse_filepath(filepath, delimiter)
            self.validate_properties(properties, prop_vals)

            file_prop = {}
            for prop, prop_val in zip(properties, prop_vals):
                file_prop[prop] = int(prop_val)
                unique_file_props[prop].add(int(prop_val))

            file_props.append(file_prop)

        return list(zip(filepaths, file_props)), unique_file_props

    def parse_filepath(self, filepath, delimiter):
        _, file = os.path.split(filepath)
        filename, _ = os.path.splitext(file)
        return filename.split(delimiter)

    def validate_properties(self, given_props, found_props):
        if len(given_props) != len(found_props):
            raise RuntimeError(
                'The number of specified properties on the base filename does not match the number of properties found')


class FileSorter:
    def __call__(self, filepaths, unique_file_props, sort_order):
        file_lists = self.split_list_on_property(filepaths, unique_file_props, sort_order[0])

        sorted_lists = []
        if len(sort_order[1:]) >= 1:
            for file_list, sort_key_val in file_lists:
                sorted_lists.append((self(file_list, unique_file_props, sort_order[1:]), sort_key_val))
        else:
            for file_list, sort_key_val in file_lists:
                ls = [filepath for filepath, file_prop in file_list]
                sorted_lists.append((ls, sort_key_val))

        sorted_lists.sort(key=lambda x: x[1])
        sorted_lists, _ = zip(*sorted_lists)
        return list(chain.from_iterable(sorted_lists))

    def split_list_on_property(self, filepaths, unique_file_props, sort_key):
        file_lists = []
        for unique_val in unique_file_props[sort_key]:
            file_list = list(filter(lambda x: x[1][sort_key] == unique_val, filepaths))
            if len(file_list):
                file_lists.append((file_list, unique_val))

        return file_lists


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


class ConformerPlannerIO(IOInterface):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'conformer_plan.txt')

    def __init__(self, peptide_length):
        self.peptide_length = peptide_length

    def load(self):
        with open(utils.attach_file_num(self.FILEPATH, self.peptide_length), 'r') as file:
            return file.readlines()

    def save(self, data):

        with open(utils.attach_file_num(self.FILEPATH, self.peptide_length), 'w') as file:
            for macrocycle_idx in data:
                file.write(str(macrocycle_idx) + '\n')


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


class JsonMacrocycleIO(AbstractJsonIO, JsonFileReader):
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

        super().__init__(utils.attach_file_num(self.FILEPATH, kwargs['peptide_length'], kwargs['job_num']),
                         file_num_range)
        super().setup(self.FILEPATH, ['peptide_length', 'job_num', 'file_num'], kwargs['peptide_length'])


class JsonConformerIO(AbstractJsonIO, JsonFileReader):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'conformers.json')

    def __init__(self, file_num_range=(None, None), **kwargs):

        super().__init__(utils.attach_file_num(self.FILEPATH, kwargs['peptide_length'], kwargs['job_num']),
                         file_num_range)
        super().setup(self.FILEPATH, ['peptide_length', 'job_num', 'file_num'], kwargs['peptide_length'])


class JsonEbejerConformerIO(AbstractJsonIO, JsonFileReader):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'ebejer.json')

    def __init__(self, file_num_range=(None, None), **kwargs):

        super().__init__(utils.attach_file_num(self.FILEPATH, kwargs['peptide_length'], kwargs['job_num']),
                         file_num_range)
        super().setup(self.FILEPATH, ['peptide_length', 'job_num', 'file_num'], kwargs['peptide_length'])


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


class JsonMWDescriptorIO(AbstractJsonIO, JsonFileReader):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'mw.json')

    def __init__(self, file_num_range=(None, None), **kwargs):

        super().__init__(utils.attach_file_num(self.FILEPATH,
                                               kwargs['peptide_length'], kwargs['job_num']), file_num_range)
        super().setup(self.FILEPATH, ['peptide_length', 'job_num', 'file_num'], kwargs['peptide_length'])


class JsonRBDescriptorIO(AbstractJsonIO, JsonFileReader):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'rb.json')

    def __init__(self, file_num_range=(None, None), **kwargs):

        super().__init__(utils.attach_file_num(self.FILEPATH,
                                               kwargs['peptide_length'], kwargs['job_num']), file_num_range)
        super().setup(self.FILEPATH, ['peptide_length', 'job_num', 'file_num'], kwargs['peptide_length'])


class JsonTPSADescriptorIO(AbstractJsonIO, JsonFileReader):

    FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'tpsa.json')

    def __init__(self, file_num_range=(None, None), **kwargs):

        super().__init__(utils.attach_file_num(self.FILEPATH,
                                               kwargs['peptide_length'], kwargs['job_num']), file_num_range)
        super().setup(self.FILEPATH, ['peptide_length', 'job_num', 'file_num'], kwargs['peptide_length'])


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
        super().__init__()

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

    def update(self, data):
        for doc in data:
            self[self.COLLECTION].find_one_and_replace({'_id': doc['_id']}, doc)


class MongoConformerIO(AbstractMongoIO):
    pass


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

        return None

    return io_closure


get_id_io = get_io(JsonIDIO, MongoIDIO)
get_index_io = get_io(JsonIndexIO, MongoIndexIO)
get_sidechain_io = get_io(JsonSideChainIO, MongoSideChainIO)
get_monomer_io = get_io(JsonMonomerIO, MongoMonomerIO)
get_peptide_io = get_io(JsonPeptideIO, MongoPeptideIO)
get_template_peptide_io = get_io(JsonTemplatePeptideIO, MongoTemplatePeptideIO)
get_macrocycle_io = get_io(JsonMacrocycleIO, MongoMacrocycleIO)
get_conformer_io = get_io(JsonConformerIO, MongoConformerIO)
get_reaction_io = get_io(JsonReactionIO, MongoReactionIO)
get_regiosqm_io = get_io(JsonRegioSQMIO, MongoRegioSQMIO)
get_pka_io = get_io(JsonpKaIO, MongopKaIO)


def get_hashed_predictions(prediction_io):
    """
    Closure for creating a function that will hash the provided predictions based on the reacting_mol attribute.

    Args:
        predictionIO (IOInterface): The uninstantiated IO class that will load the predictions.

    Returns:
        func: A function that can be called in order to recieve the hashed predictions.
    """

    def hasher():
        hashed_predictions = {}
        for prediction in prediction_io().load():
            hashed_predictions[prediction['reacting_mol']] = prediction['predictions']
        return hashed_predictions

    return hasher


get_hashed_regiosqm_predictions = get_hashed_predictions(get_regiosqm_io)
get_hashed_pka_predictions = get_hashed_predictions(get_pka_io)


def get_filtered_monomer_set():
    sdf_io = StructureDataFileIO(os.path.join(config.DATA_DIR, 'chemdraw', 'conservative_sidechains.sdf'))
    sidechains = get_sidechain_io().load()
    hashed_sidechains = hash_molecules(sidechains, 'shared_id', key='kekule')
    monomers = get_monomer_io().load()

    shared_ids = []
    for molecule in sdf_io.load():
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(molecule))
        Chem.Kekulize(mol)
        shared_ids.append(hashed_sidechains[Chem.MolToSmiles(mol, kekuleSmiles=True)]['shared_id'])

    filtered_monomers = []
    for monomer in monomers:
        if monomer['sidechain'] in shared_ids or monomer['imported']:
            filtered_monomers.append(monomer)

    return filtered_monomers


def hash_molecules(mols, *args, key='_id'):
    hashed_mols = defaultdict(dict)
    for mol in mols:
        for k, v in mol.items():
            if k in args:
                hashed_mols[mol[key]].update({k: v})

    return hashed_mols
