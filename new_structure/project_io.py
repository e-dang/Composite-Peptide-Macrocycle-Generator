import json
import os
from abc import ABC, abstractmethod

from bson import json_util

import config
import utils


class IOInterface(ABC):
    """
    Interface for classes that handle IO of different molecule and reaction data.
    """

    @abstractmethod
    def load(self, file_num_range):
        """
        Abstract method for loading data handled by the specific derived class.

        Args:
            file_num_range (tuple[int]): A tuple containing two intergers specifying the range of file numbers.
        """

    @abstractmethod
    def save(self, data):
        """
        Abstract method for saving the provided data in the location maintained by the derived class.

        Args:
            data (iterable[dict]): The data to be saved.
        """


class AbstractJsonIO(IOInterface):
    """
    Abstract class for IO classes that read and write to json file format.
    """

    def from_json(self, filepath, low, high):
        """
        Loads the data in the json file pointed to by filepath. The file name specified in filepath is treated as a
        base file name, meaning there could be a set of files within the directory with the same base file name
        distinquished only by a number before the extension. For example:

        base file name - test.json
        actual file names - test_0.json, test_1.json, ...

        The arguments low and high dictate which files containing the base file name should be loaded. If the argument
        low is not provided, this method will load data from all the numbered files in the directory with the provided
        base file name. If argument high is not provided then it will load the data in the base file name marked with
        the specified low number. If both arguments are provided then it will load all data in the files with the base
        file name that are distinquished by the numbers within the range of low to high.

        Args:
            filepath (str): The path to the files containing the data to be loaded. The specific file name specified
                here should be the base file name.
            low (int): The lower end of the range of file numbers attached to the base file name.
            high (int): The high end of the range of file numbers attached to the base file name.

        Yields:
            dict: A dictionary containing an entry of data in the json file.
        """

        # determine range based on the provided low and high arguments
        if low is None:
            low, high = utils.get_file_num_range(filepath)
        elif high is None:
            high = low + 1

        # yield all data from within the specified range
        for file_num in range(low, high):
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


class JsonHeterocycleIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling heterocycle data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'heterocycles.json')

    def load(self, file_num_range=(None, None)):

        try:
            return super().from_json(self._FILEPATH, file_num_range[0], file_num_range[1])
        except ValueError:
            raise ValueError('file_num_range must be a tuple containing two integers')

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonConnectionsIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling connection data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'connections.json')

    def load(self, file_num_range=(None, None)):

        try:
            return super().from_json(self._FILEPATH, file_num_range[0], file_num_range[1])
        except ValueError:
            raise ValueError('file_num_range must be a tuple containing two integers')

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonSideChainIO(AbstractJsonIO):
    """
    Implmentation of the AbstractJsonIO class for handling side_chain data.
    """

    _FILEPATH = os.path.join(config.DATA_DIR, 'generated', 'side_chains.json')

    def load(self, file_num_range=(None, None)):

        try:
            return super().from_json(self._FILEPATH, file_num_range[0], file_num_range[1])
        except ValueError:
            raise ValueError('file_num_range must be a tuple containing two integers')

    def save(self, data):

        super().to_json(self._FILEPATH, data)


class JsonMonomerIO(AbstractJsonIO):

    def load(self, file_num_range=(None, None)):
        pass

    def save(self, data):
        pass


class JsonPeptideIO(AbstractJsonIO):

    def load(self, file_num_range=(None, None)):
        pass

    def save(self, data):
        pass


class JsonTemplatePeptideIO(AbstractJsonIO):

    def load(self, file_num_range=(None, None)):
        pass

    def save(self, data):
        pass


class JsonMacrocycleIO(AbstractJsonIO):

    def load(self, file_num_range=(None, None)):
        pass

    def save(self, data):
        pass


class AbstractMongoIO(IOInterface):

    def to_mongo(self, collection):
        pass

    def from_mongo(self, collection):
        pass
