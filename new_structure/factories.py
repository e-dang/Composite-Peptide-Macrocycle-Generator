import multiprocessing
from abc import ABC, abstractmethod

import data_handlers
import generators


class IFactory(ABC):
    """
    Interface for classes that wrap the functionality of Generator and DataHandler classes together such that it
    passes the data loaded by the DataHandler class to the Generator class, and passes those results back to the
    DataHandler so it can be saved.
    """

    @abstractmethod
    def run(self, factory_arg):
        """
        Abstract method for wrapping the functionality of the Generator and DataHandler classes present in
        factory_arg such that it is executed in a parallel fashion using the multiprossesing module.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """

    @abstractmethod
    def run_serial(self, factory_arg):
        """
        Abstract method for wrapping the functionality of the Generator and DataHandler classes present in the
        factory_arg such that it is performed in a serial fashion.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """


class MolFactory(IFactory):
    """
    An implementation of IFactory that takes FactoryArguments that contains a Generator.
    """

    def run(self, factory_arg):

        generator, handler = self.split_args(factory_arg)

        result_data = []
        with multiprocessing.Pool() as pool:
            for result_mols in pool.imap_unordered(generator.generate, generator.get_args(handler.load())):
                for result_mol in result_mols:
                    result_data.append(result_mol)

        handler.save(result_data)

    def run_serial(self, factory_arg):

        generator, handler = self.split_args(factory_arg)

        result_data = []
        for arg in generator.get_args(handler.load()):
            for result_mol in generator.generate(arg):
                result_data.append(result_mol)

        handler.save(result_data)

    def split_args(self, factory_arg):
        """
        Helper method that extracts the Generator and DataHandler classes from the FactoryArgument class.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.

        Raises:
            AttributeError: Raised when the provided argument doesn't contain an attribute of 'generator' and/or
                'handler'.

        Returns:
            tuple: The Generator and DataHandler
        """

        try:
            return factory_arg.generator, factory_arg.handler
        except AttributeError:
            raise AttributeError(
                'The supplied factory_arg must have instance members \'generator\', \'loader\', and \'saver\'')


class FilterFactory(IFactory):
    def run(self):
        pass

    def run_serial(self):
        pass


class IFactoryArgument(ABC):
    """
    An interface for classes that wrap Generator and DataHandler classes together to be used by a Factory.
    """

    @abstractmethod
    def __init__(self):
        """
        An intializer method for creating the Generator and DataHandler classes.
        """


class SCConnectionModificationArgs(IFactoryArgument):
    """
    An implementation of a FactoryArgument containing the classes required to combine heterocycles and connection
    molecules into sidechain molecules.
    """

    def __init__(self):

        self.handler = data_handlers.SCCMDataHandler()
        self.generator = generators.SideChainConnectionModifier(self.handler.id_iterator)


class MonomerGenerationArgs(IFactoryArgument):
    """
    An implementation of a FactoryArgument containing the classes required to combine sidechains and backbone
    molecules into monomers.
    """

    def __init__(self):

        self.handler = data_handlers.MGDataHandler()
        self.generator = generators.MonomerGenerator(self.handler.index_iterator)


class PeptideGenerationArgs(IFactoryArgument):

    def __init__(self):

        self.generator = generators.PeptideGenerator()
        self.handler = data_handlers.PGDataHandler()


class TemplatePeptideGenerationArgs(IFactoryArgument):

    def __init__(self):

        self.generator = generators.TemplatePeptideGenerator()
        self.handler = data_handlers.TPGDataHandler()


class MacrocycleGenerationArgs(IFactoryArgument):

    def __init__(self):

        self.generator = generators.MacrocycleGenerator()
        self.handler = data_handlers.MCGDataHandler()


class UniMolecularReactionGenerationArgs(IFactoryArgument):

    def __init__(self):

        self.generator = generators.UniMolecularReactionGenerator()
        self.handler = data_handlers.UMRGDataHandler()


class BiMolecularReactionGenerationArgs(IFactoryArgument):

    def __init__(self):

        self.generator = generators.BiMolecularReactionGenerator()
        self.handler = data_handlers.BMRGDataHandler()
