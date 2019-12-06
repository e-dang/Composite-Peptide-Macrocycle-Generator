from abc import ABC, abstractmethod
from collections import namedtuple

import mol_io


class IDataSupplier(ABC):
    """
    Interface for classes that couple together specific IO classes in order to supply data for a specific
    MolTransformer.
    """

    @abstractmethod
    def load(self):
        """
        Abstract method for calling load on each IO class of the specific derived DataSupplier and returning the data
        in a tuple.
        """


class SCGDataSupplier(IDataSupplier):
    """
    Implementation of a DataSupplier used to supply the data needed by the SideChainGenerator MolTransformer.
    """

    def __init__(self, source):
        """
        Initializer method used to instaniate the MolIO classes for loading heterocycles and connection molecules from
        either JSON files or MongoDB depending on the specified source.

        Args:
            source (str): The source to load data from. Can either be 'json' or 'mongo'
        """

        if source == 'json':
            self.heterocycles = mol_io.JsonHeterocycleIO()
            self.connections = mol_io.JsonConnectionsIO()
        elif source == 'mongo':
            pass

    def load(self):

        SideChainGeneratorData = namedtuple('SideChainGeneratorData', 'heterocycles connections')
        return SideChainGeneratorData(self.heterocycles.load(), self.connections.load())
