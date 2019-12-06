from abc import ABC, abstractmethod
import multiprocessing
import transformers
import mol_io
import loaders


class IFactory(ABC):
    """
    Interface for classes that wrap the functionality of Transformer, DataSupplier, and IO classes together such that it
    passes the data loaded by the DataSupplier class to the Transformer class, and passes those results to the IO class
    so it can be saved.
    """

    @abstractmethod
    def run(self, factory_arg):
        """
        Abstract method for wrapping the functionality of the Transformer, DataSupplier, and IO class present in
        factory_arg such that it is executed in a parallel fashion using the multiprossesing module.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """

    @abstractmethod
    def run_serial(self, factory_arg):
        """
        Abstract method for wrapping the functionality of the Transformer, DataSupplier, and IO class present in the
        factory_arg such that it is performed in a serial fashion.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """


class MolFactory(IFactory):
    """
    An implementation of IFactory that takes FactoryArguments that contains a MolTransformer.
    """

    def run(self, factory_arg):

        mol_transformer, loader, saver = self.split_args(factory_arg)

        result_data = []
        with multiprocessing.Pool() as pool:
            for result_mols in pool.imap_unordered(mol_transformer.transform, mol_transformer.get_args(loader.load())):
                for result_mol in result_mols:
                    result_data.append(result_mol)

        saver.save(result_data)

    def run_serial(self, factory_arg):

        mol_transformer, loader, saver = self.split_args(factory_arg)

        result_data = []
        for arg in mol_transformer.get_args(loader.load()):
            for result_mol in mol_transformer.transform(arg):
                result_data.append(result_mol)

        saver.save(result_data)

    def split_args(self, factory_arg):
        """
        Helper method that extracts the MolTransformer, DataSupplier, and MolIO classes from the FactoryArgument class.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.

        Raises:
            AttributeError: Raised when the provided argument doesn't contain an attribute of transformer, loader,
            and/or saver.

        Returns:
            tuple: The MolTransformer, DataSupplier, and MolIO class
        """

        try:
            return factory_arg.transformer, factory_arg.loader, factory_arg.saver
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
    An interface for classes that wrap a Transformer, DataSupplier, and IO class to be used by a Factory.
    """

    @abstractmethod
    def __init__(self, source):
        """
        An intializer method for creating the Transformer, DataSupplier, and IO classes.

        Args:
            source (str): The format of the data to load/save to.
        """


class SideChainGenerationArgs(IFactoryArgument):
    """
    An implementation of a FactoryArgument containing the classes required to combine heterocycles and connection
    molecules into side_chain molecules.
    """

    def __init__(self, source):

        self.transformer = transformers.SideChainGenerator()
        self.loader = loaders.SCGDataSupplier(source)
        self.saver = mol_io.JsonSideChainIO()
