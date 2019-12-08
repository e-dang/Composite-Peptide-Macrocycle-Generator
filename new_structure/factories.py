from abc import ABC, abstractmethod
import multiprocessing
import transformers
import data_handlers


class IFactory(ABC):
    """
    Interface for classes that wrap the functionality of Transformer and DataHandler classes together such that it
    passes the data loaded by the DataHandler class to the Transformer class, and passes those results back to the
    DataHandler so it can be saved.
    """

    @abstractmethod
    def run(self, factory_arg):
        """
        Abstract method for wrapping the functionality of the Transformer and DataHandler classes present in
        factory_arg such that it is executed in a parallel fashion using the multiprossesing module.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """

    @abstractmethod
    def run_serial(self, factory_arg):
        """
        Abstract method for wrapping the functionality of the Transformer and DataHandler classes present in the
        factory_arg such that it is performed in a serial fashion.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """


class MolFactory(IFactory):
    """
    An implementation of IFactory that takes FactoryArguments that contains a MolTransformer.
    """

    def run(self, factory_arg):

        mol_transformer, handler = self.split_args(factory_arg)

        result_data = []
        with multiprocessing.Pool() as pool:
            for result_mols in pool.imap_unordered(mol_transformer.transform, mol_transformer.get_args(handler.load())):
                for result_mol in result_mols:
                    result_data.append(result_mol)

        handler.save(result_data)

    def run_serial(self, factory_arg):

        mol_transformer, handler = self.split_args(factory_arg)

        result_data = []
        for arg in mol_transformer.get_args(handler.load()):
            for result_mol in mol_transformer.transform(arg):
                result_data.append(result_mol)

        handler.save(result_data)

    def split_args(self, factory_arg):
        """
        Helper method that extracts the MolTransformer and DataHandler classes from the FactoryArgument class.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.

        Raises:
            AttributeError: Raised when the provided argument doesn't contain an attribute of 'transformer' and/or
                'handler'.

        Returns:
            tuple: The MolTransformer and DataHandler
        """

        try:
            return factory_arg.transformer, factory_arg.handler
        except AttributeError:
            raise AttributeError(
                'The supplied factory_arg must have instance members \'transformer\', \'loader\', and \'saver\'')


class ReactionFactory(IFactory):
    def run(self):
        pass

    def run_serial(self):
        pass


class FilterFactory(IFactory):
    def run(self):
        pass

    def run_serial(self):
        pass


class IFactoryArgument(ABC):
    """
    An interface for classes that wrap Transformer and DataHandler classes together to be used by a Factory.
    """

    @abstractmethod
    def __init__(self, data_format):
        """
        An intializer method for creating the Transformer and DataHandler classes.

        Args:
            data_format (str): The format of the data to load from and save to.
        """


class SideChainGenerationArgs(IFactoryArgument):
    """
    An implementation of a FactoryArgument containing the classes required to combine heterocycles and connection
    molecules into sidechain molecules.
    """

    def __init__(self, data_format):

        self.transformer = transformers.SideChainGenerator()
        self.handler = data_handlers.SCGDataHandler(data_format)


class MonomerGenerationArgs(IFactoryArgument):
    """
    An implementation of a FactoryArgument containing the classes required to combine sidechains and backbone
    molecules into monomers.
    """

    def __init__(self, data_format):

        self.transformer = transformers.MonomerGenerator()
        self.handler = data_handlers.MGDataHandler(data_format)


class PeptideGenerationArgs(IFactoryArgument):

    def __init__(self, data_format):

        self.transformer = transformers.PeptideGenerator()
        self.handler = data_handlers.PGDataHandler(data_format)


class TemplatePeptideGenerationArgs(IFactoryArgument):

    def __init__(self, data_format):

        self.transformer = transformers.TemplatePeptideGenerator()
        self.handler = data_handlers.TPGDataHandler(data_format)


class ReactionGenerationArgs(IFactoryArgument):

    def __init__(self, data_format):

        self.transformer = transformers.JointReactionGenerator()
        self.handler = data_handlers.JRGDataHandler(data_format)
