from abc import ABC, abstractmethod
from collections import namedtuple

import project_io
import reactions
import utils


class IDataHandler(ABC):
    """
    Interface for classes that couple together specific IO classes in order to supply and write data for a specific
    Generator class.
    """

    @abstractmethod
    def __init__(self, data_format):
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


class SCGDataHandler(IDataHandler):
    """
    Implementation of a DataHandler used to supply and write the data needed and generated by the SideChainGenerator
    Generator.
    """

    def __init__(self, data_format):

        self._connections = utils.get_connections()
        if data_format == 'json':
            self._heterocycle_loader = project_io.JsonHeterocycleIO()
            self._sidechain_saver = project_io.JsonSideChainIO()
        elif data_format == 'mongo':
            self._heterocycle_loader = project_io.MongoHeterocycleIO()
            self._sidechain_saver = project_io.MongoSideChainIO()

    def load(self):

        SideChainGeneratorData = namedtuple('SideChainGeneratorData', 'heterocycles connections')
        return SideChainGeneratorData(self._heterocycle_loader.load(), self._connections)

    def save(self, data):

        self._sidechain_saver.save(data)


class MGDataHandler(IDataHandler):
    """
    Implementation of a DataHandler used to supply and write the data needed and generated by the MonomerGenerator
    Generator.
    """

    def __init__(self, data_format):

        self._backbones = utils.get_backbones()
        if data_format == 'json':
            self._sidechain_loader = project_io.JsonSideChainIO()
            self._monomer_saver = project_io.JsonMonomerIO()
        elif data_format == 'mongo':
            self._sidechain_loader = project_io.MongoSideChainIO()
            self._monomer_saver = project_io.MongoMonomerIO()

    def load(self):

        MonomerGeneratorData = namedtuple('MonomerGeneratorData', 'sidechains backbones')
        return MonomerGeneratorData(self._sidechain_loader.load(), self._backbones)

    def save(self, data):

        self._monomer_saver.save(data)


class PGDataHandler(IDataHandler):
    """
    Implementation of a DataHandler used to supply and write the data needed and generated by the PeptideGenerator
    Generator.
    """

    def __init__(self, data_format):

        if data_format == 'json':
            self._monomer_loader = project_io.JsonMonomerIO()
            self._peptide_saver = project_io.JsonPeptideIO()
        elif data_format == 'mongo':
            self._monomer_loader = project_io.MongoMonomerIO()
            self._peptide_saver = project_io.MongoPeptideIO()

    def load(self):
        PeptideGeneratorData = namedtuple('PeptideGeneratorData', 'monomers peptide_length')
        return PeptideGeneratorData(self._monomer_loader.load(), 3)

    def save(self, data):

        self._peptide_saver.save(data)


class TPGDataHandler(IDataHandler):

    def __init__(self, data_format):

        self._templates = utils.get_templates()
        if data_format == 'json':
            self._peptide_loader = project_io.JsonPeptideIO()
            self._template_peptide_saver = project_io.JsonTemplatePeptideIO()
        elif data_format == 'mongo':
            self._peptide_loader = project_io.MongoPeptideIO()
            self._template_peptide_saver = project_io.MongoTemplatePeptideIO()

    def load(self):

        TemplatePeptideGeneratorData = namedtuple('TemplatePeptideGeneratorData', 'peptides templates')
        return TemplatePeptideGeneratorData(self._peptide_loader.load(), self._templates)

    def save(self, data):

        self._template_peptide_saver.save(data)


class BMRGDataHandler(IDataHandler):

    def __init__(self, data_format):

        self._templates = utils.get_templates()
        self._reactions = utils.get_bimolecular_reactions()
        if data_format == 'json':
            self._sidechain_loader = project_io.JsonSideChainIO()
            self._reaction_saver = project_io.JsonReactionIO()
        elif data_format == 'mongo':
            self._sidechain_loader = project_io.MongoSideChainIO()
            self._reaction_saver = project_io.MongoReactionIO()

    def load(self):

        BiMolecularReactionGeneratorData = namedtuple(
            'BiMolecularReactionGeneratorData', 'sidechains templates reactions')
        return BiMolecularReactionGeneratorData(self._sidechain_loader.load(), self._templates, self._reactions)

    def save(self, data):

        self._reaction_saver.save(data)