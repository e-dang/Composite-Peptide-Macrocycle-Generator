from abc import ABC, abstractmethod
from collections import namedtuple

import config
import iterators
import project_io
import utils


class IDataHandler(ABC):
    """
    Interface for classes that couple together specific IO classes in order to supply and write data for a specific
    Generator class.
    """

    @abstractmethod
    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):
        """
        Abstract initializer used to instantiate the IO classes for loading and saving the data required and generated
        by the derived DataHandler's coupled Generator class.

        Args:
            data_format (str): The format of the data to load from and save to. Can either be 'json' or 'mongo'
        """

    @abstractmethod
    def load(self):
        """
        Abstract method for calling load on each IO class of the specific derived DataHandler and returning the data
        in a tuple.
        """

    @abstractmethod
    def save(self, data):
        """
        Abstract method for calling save() on the IO class used for saving data in the specific derived DatSupplier.
        """


class SCCMDataHandler(IDataHandler):
    """
    Implementation of a DataHandler used to supply and write the data needed and generated by the SideChainGenerator
    Generator.
    """

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        if data_format == 'json':
            id_io = project_io.JsonIDIO()
            self._sidechain_io = project_io.JsonSideChainIO()
        elif data_format == 'mongo':
            id_io = project_io.MongoIDIO()
            self._sidechain_io = project_io.MongoSideChainIO()

        self.id_iterator = iterators.IDIterator(id_io)

    def load(self):

        return self._sidechain_io.load()

    def save(self, data):

        self._sidechain_io.save(data)
        self.id_iterator.save()


class MGDataHandler(IDataHandler):
    """
    Implementation of a DataHandler used to supply and write the data needed and generated by the MonomerGenerator
    Generator.
    """

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        if data_format == 'json':
            index_io = project_io.JsonIndexIO()
            self._sidechain_loader = project_io.JsonSideChainIO()
            self._monomer_saver = project_io.JsonMonomerIO()
        elif data_format == 'mongo':
            index_io = project_io.MongoIndexIO()
            self._sidechain_loader = project_io.MongoSideChainIO()
            self._monomer_saver = project_io.MongoMonomerIO()

        self.index_iterator = iterators.IndexIterator(index_io)

    def load(self):

        return self._sidechain_loader.load()

    def save(self, data):

        self._monomer_saver.save(data)
        self.index_iterator.save()


class TestPGDataHandler(IDataHandler):
    """
    Implementation of a DataHandler used to supply and write the data needed and generated by the PeptideGenerator
    Generator.
    """

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        if data_format == 'json':
            self._monomer_loader = project_io.JsonMonomerIO()
            self._peptide_saver = project_io.JsonPeptideIO()
        elif data_format == 'mongo':
            self._monomer_loader = project_io.MongoMonomerIO()
            self._peptide_saver = project_io.MongoPeptideIO()

        self.peptide_length = kwargs['peptide_length']

    def load(self):
        return self._monomer_loader.load(), self.peptide_length

    def save(self, data):

        self._peptide_saver.save(data)


class PublicationPGDataHandler(IDataHandler):
    """
    Implementation of a DataHandler used to supply and write the data needed and generated by the PeptideGenerator
    Generator.
    """

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        if data_format == 'json':
            self._peptide_saver = project_io.JsonPeptideIO()
        elif data_format == 'mongo':
            self._peptide_saver = project_io.MongoPeptideIO()

        self.peptide_length = kwargs['peptide_length']
        self.plan_loader = project_io.PeptidePlannerIO(self.peptide_length)

    def load(self):
        return self.plan_loader.load()

    def save(self, data):

        self._peptide_saver.save(data, peptide_length=self.peptide_length)


class TPGDataHandler(IDataHandler):

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        if data_format == 'json':
            self._peptide_loader = project_io.JsonPeptideIO()
            self._template_peptide_saver = project_io.JsonTemplatePeptideIO()
        elif data_format == 'mongo':
            self._peptide_loader = project_io.MongoPeptideIO()
            self._template_peptide_saver = project_io.MongoTemplatePeptideIO()

        self.peptide_length = kwargs['peptide_length']

    def load(self):

        return self._peptide_loader.load(peptide_length=self.peptide_length)

    def save(self, data):

        self._template_peptide_saver.save(data, peptide_length=self.peptide_length)


class MCGDataHandler(IDataHandler):

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        if data_format == 'json':
            self._template_peptide_loader = project_io.JsonTemplatePeptideIO()
            self._reaction_loader = project_io.JsonReactionIO()
            self._macrocycle_saver = project_io.JsonMacrocycleIO()
        elif data_format == 'mongo':
            self._template_peptide_loader = project_io.MongoTemplatePeptideIO()
            self._reaction_loader = project_io.MongoReactionIO()
            self._macrocycle_saver = project_io.MongoMacrocycleIO()

        self.peptide_length = kwargs['peptide_length']
        self.start = kwargs['start']
        self.end = kwargs['end']

    def load(self):
        MacrocycleGeneratorData = namedtuple('MacrocycleGeneratorData', 'template_peptides reactions')

        return MacrocycleGeneratorData(self.load_template_peptides(), self._reaction_loader.load())

    def save(self, data):

        self._macrocycle_saver.save(data)

    def load_template_peptides(self):

        for i, template_peptide in enumerate(self._template_peptide_loader.load(peptide_length=self.peptide_length)):
            if i < self.start:
                continue
            elif i >= self.end:
                break
            else:
                yield template_peptide


class UMRGDataHandler(IDataHandler):

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        self._templates = utils.get_templates()
        self._reactions = utils.get_unimolecular_reactions()
        if data_format == 'json':
            self._reaction_saver = project_io.JsonReactionIO()
        elif data_format == 'mongo':
            self._reaction_saver = project_io.MongoReactionIO()

    def load(self):

        return self._templates, self._reactions

    def save(self, data):

        self._reaction_saver.save(data)


class BMRGDataHandler(IDataHandler):

    def __init__(self, data_format=config.DATA_FORMAT, **kwargs):

        self._reactions = utils.get_bimolecular_reactions()
        if data_format == 'json':
            self._sidechain_loader = project_io.JsonSideChainIO()
            self._monomer_loader = project_io.JsonMonomerIO()
            self._reaction_saver = project_io.JsonReactionIO()
        elif data_format == 'mongo':
            self._sidechain_loader = project_io.MongoSideChainIO()
            self._monomer_loader = project_io.MongoMonomerIO()
            self._reaction_saver = project_io.MongoReactionIO()

    def load(self):

        reacting_mols = list(filter(lambda x: x['connection'] == 'methyl', self._sidechain_loader.load()))
        reacting_mols.extend(list(filter(lambda x: x['required'] and x['imported'], self._monomer_loader.load())))
        return reacting_mols, self._reactions

    def save(self, data):

        self._reaction_saver.save(data)
