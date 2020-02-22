import multiprocessing
from abc import ABC, abstractmethod
from functools import partial

from rdkit import Chem
from rdkit.Chem import AllChem

import macrocycles.descriptors as descriptors
import macrocycles.argument_producers as argument_producers
import macrocycles.config as config
import macrocycles.data_handlers as data_handlers
import macrocycles.generators as generators


class IFactory(ABC):
    """
    Interface for classes that wrap the functionality of Generator, DataHandler, and ArgumentProducer classes together
    such that it passes the data loaded by the DataHandler class to the ArgumentProducer which then passes data to the
    Generator class, and passes those results back to the DataHandler so it can be saved.
    """

    @abstractmethod
    def run(self, factory_arg):
        """
        Abstract method for executing the functionality of the Generator, DataHandler, and ArgumentProducer classes
        present in the factory_arg such that it is executed in a parallel fashion using the multiprossesing module.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """

    @abstractmethod
    def run_serial(self, factory_arg):
        """
        Abstract method for executing the functionality of the Generator, DataHandler, and ArgumentProducer classes
        present in the factory_arg such that it is performed in a serial fashion.

        Args:
            factory_arg (FactoryArgument): A subclass of the FactoryArgument interface.
        """


class MolFactory(IFactory):
    """
    An implementation of IFactory that takes FactoryArguments that contains a Generator.
    """

    def run(self, factory_arg):

        generator, data_handler, arg_producer = self.split_args(factory_arg)

        self.count = 0
        self.result_data = []
        with multiprocessing.Pool(config.NUM_PROCS, maxtasksperchild=config.TASKS_PER_CHILD) as pool:
            for result_mols in pool.imap_unordered(generator.generate, arg_producer(data_handler.load())):
                for result_mol in result_mols:
                    self.result_data.append(result_mol)
                    self.checkpoint(data_handler)

        data_handler.save(self.result_data)
        self.count += len(self.result_data)

    def run_serial(self, factory_arg):

        generator, data_handler, arg_producer = self.split_args(factory_arg)

        self.count = 0
        self.result_data = []
        for arg in arg_producer(data_handler.load()):
            for result_mol in generator.generate(arg):
                self.result_data.append(result_mol)
                self.checkpoint(data_handler)

        data_handler.save(self.result_data)
        self.count += len(self.result_data)

    def split_args(self, factory_arg):
        """
        Helper method that extracts the Generator and DataHandler classes from the FactoryArgument class.

        Args:
            factory_arg (IFactoryArgument): A subclass of the IFactoryArgument interface.

        Raises:
            AttributeError: Raised when the provided argument doesn't contain an attribute of 'generator',
                'arg_producer', and/or 'data_handler'.

        Returns:
            tuple: The Generator and DataHandler
        """

        try:
            return factory_arg.generator, factory_arg.data_handler, factory_arg.arg_producer
        except AttributeError:
            raise AttributeError(
                'The supplied factory_arg must have instance members \'generator\', \'arg_producer\', and '
                '\'data_handler\'')

    def checkpoint(self, data_handler):
        """
        Helper method that saves the result_data once it reaches a certain capacity and resets result_data to an empty
        list.

        Args:
            data_handler (IDataHandler): An implementation of IDataHandler.
        """

        if len(self.result_data) >= config.CAPACITY:
            data_handler.save(self.result_data)
            self.count += len(self.result_data)
            self.result_data = []


class DescriptorFactory(IFactory):

    def run(self, factory_arg):

        descriptor, data_handler = self.split_args(factory_arg)

        self.count = 0
        self.result_data = []
        with multiprocessing.Pool(config.NUM_PROCS, maxtasksperchild=config.TASKS_PER_CHILD) as pool:
            for value in pool.imap_unordered(descriptor, data_handler.load()):
                self.result_data.append(value)

        data_handler.save([{descriptor.NAME: self.result_data}])
        self.count += len(self.result_data)

    def run_serial(self, factory_arg):

        descriptor, data_handler = self.split_args(factory_arg)

        self.count = 0
        self.result_data = []
        for mol in data_handler.load():
            self.result_data.append(descriptor(mol))

        data_handler.save([{descriptor.NAME: self.result_data}])
        self.count += len(self.result_data)

    def split_args(self, factory_arg):

        return factory_arg.descriptor, factory_arg.data_handler


class IFactoryArgument(ABC):
    """
    An interface for classes that wrap Generator, DataHandler, and ArgumentProducer classes together to be used by a
    Factory.
    """

    @abstractmethod
    def __init__(self, data_handler, arg_producer, generator):
        """
        Intializer method that assigns a Generator, DataHandler, and ArgumentProducer to instance variables.
        """

        self.data_handler = data_handler
        self.arg_producer = arg_producer
        self.generator = generator


class SCConnectionModificationArgs(IFactoryArgument):
    """
    An implementation of a IFactoryArgument containing the classes required to create new sidechain molecules from
    imported sidechains by changing the connection type.
    """

    def __init__(self, **kwargs):
        """
        Initializer method that creates instances of SCCMDataHandler, NullArgProducer, and SideChainConnectionModifier
        initialized with SCCMDataHandler's id_iterator.
        """

        data_handler = data_handlers.SCCMDataHandler(**kwargs)
        super().__init__(data_handler,
                         argument_producers.NullArgProducer(),
                         generators.SideChainConnectionModifier(data_handler.id_iterator))


class MonomerGenerationArgs(IFactoryArgument):
    """
    An implementation of IFactoryArgument containing the classes required to combine sidechains and backbone
    molecules into monomers.
    """

    def __init__(self, **kwargs):
        """
        Initializer method that creates instances of MGDataHandler, NullArgProducer, and MonomerGenerator
        initialized with MGDataHandler's id_iterator.
        """

        data_handler = data_handlers.MGDataHandler(**kwargs)
        super().__init__(data_handler,
                         argument_producers.NullArgProducer(),
                         generators.MonomerGenerator(data_handler.index_iterator))


class PeptideGenerationArgs(IFactoryArgument):
    """
    An implementation of IFactoryArgument containing the classes required to combine monomers into peptides.
    """

    def __init__(self, **kwargs):
        """
        Initializer method that creates instances of PGDataHandler, PeptideGeneratorArgProducer, and PeptideGenerator.

        Args:
            peptide_length (int): Keyword argument that defines the length of peptide to generate.
        """

        super().__init__(data_handlers.PGDataHandler(**kwargs),
                         argument_producers.PeptideGeneratorArgProducer(),
                         generators.PeptideGenerator(kwargs['peptide_length']))


class TemplatePeptideGenerationArgs(IFactoryArgument):
    """
    An implementation of IFactoryArgument containing the classes required to combine peptides and templates into
    template_peptide oligomers.
    """

    def __init__(self, **kwargs):
        """
        Initializer method that creates instances of TPGDataHandler, NullArgProducer, and TemplatePeptideGenerator.

        Args:
            peptide_length (int): The length of the peptides that should be loaded and combined into template_peptides.
                This keyword argument is only required if DATA_FORMAT is json.
        """

        super().__init__(data_handlers.TPGDataHandler(**kwargs),
                         argument_producers.NullArgProducer(),
                         generators.TemplatePeptideGenerator())


class MacrocycleGenerationArgs(IFactoryArgument):
    """
    An implementation of IFactoryArgument containing the classes required to create macrocycles from template_peptides
    and reactions.
    """

    def __init__(self, **kwargs):
        """
        Initializer method that creates instances of MCGDataHandler, MacrocycleGeneratorArgProducer, and
        MacrocycleGenerator.

        Args:
            peptide_length (int): The length of the peptides in the template_peptides that should be loaded and
                transformed into macrocycles. This keyword argument is only required if DATA_FORMAT is json.
        """

        super().__init__(data_handlers.MCGDataHandler(**kwargs),
                         argument_producers.MacrocycleGeneratorArgProducer(),
                         generators.MacrocycleGenerator())


class ConformerGenerationArgs(IFactoryArgument):

    def __init__(self, **kwargs):

        super().__init__(data_handlers.ConformerGeneratorDataHandler(**kwargs),
                         argument_producers.ConformerGeneratorArgProducer(**kwargs),
                         generators.MacrocycleConformerGenerator())


class EbejerConformerGenerationArgs(IFactoryArgument):

    def __init__(self, **kwargs):

        super().__init__(data_handlers.ConformerGeneratorDataHandler(**kwargs),
                         argument_producers.NullArgProducer(),
                         generators.EbejerConformerGenerator())


class UniMolecularReactionGenerationArgs(IFactoryArgument):
    """
    An implementation of IFactoryArgument containing the classes required to create UniMolecularReactions.
    """

    def __init__(self, **kwargs):
        """
        Initializer method that creates instances of UMRGDataHandler, CartesianProductArgProducer, and
        UniMolecularReactionGenerator.
        """

        super().__init__(data_handlers.UMRGDataHandler(**kwargs),
                         argument_producers.CartesianProductArgProducer(),
                         generators.UniMolecularReactionGenerator())


class BiMolecularReactionGenerationArgs(IFactoryArgument):
    """
    An implementation of IFactoryArgument containing the classes required to create BiMolecularReactions.
    """

    def __init__(self, **kwargs):
        """
        Initializer method that creates instances of BMRGDataHandler, CartesianProductArgProducer, and
        BiMolecularReactionGenerator.
        """
        super().__init__(data_handlers.BMRGDataHandler(**kwargs),
                         argument_producers.CartesianProductArgProducer(),
                         generators.BiMolecularReactionGenerator())


class MWDescriptorArgs:

    def __init__(self, **kwargs):

        self.descriptor = descriptors.MolWeightDescriptor()
        self.data_handler = data_handlers.MWDescriptorDataHandler(**kwargs)


class RBDescriptorArgs:

    def __init__(self, **kwargs):

        self.descriptor = descriptors.RotatableBondsDescriptor()
        self.data_handler = data_handlers.RBDescriptorDataHandler(**kwargs)


class TPSADescriptorArgs:

    def __init__(self, **kwargs):

        self.descriptor = descriptors.TPSADescriptor()
        self.data_handler = data_handlers.TPSADescriptorDataHandler(**kwargs)
