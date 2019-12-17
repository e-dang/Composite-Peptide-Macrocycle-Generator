from abc import ABC, abstractmethod

from rdkit import Chem

import iterators
import project_io


class IDataInitializer(ABC):

    @abstractmethod
    def initialize(self):
        pass


class DataInitializer(IDataInitializer):

    def __init__(self, data_format):

        self.data_format = data_format

    def initialize(self):

        self.phase_1()
        self.phase_2()

    def phase_1(self):

        self.record_initializers = [IDInitializer(self.data_format), IndexInitializer(self.data_format)]
        self.select_initialization(self.record_initializers)

    def phase_2(self):

        self.id_iterator = iterators.IDIterator(self.data_format)
        self.index_iterator = iterators.IndexIterator(self.data_format)
        self.mol_initializers = [SideChainInitializer(self.data_format, self.id_iterator),
                                 MonomerInitializer(self.data_format, self.id_iterator, self.index_iterator)]
        self.select_initialization(self.mol_initializers)

    def select_initialization(self, initializers):

        for initializer in initializers:
            initializer.initialize()


class SideChainInitializer(IDataInitializer):

    def __init__(self, data_format, id_iterator):
        if data_format == 'json':
            self.saver = project_io.JsonSideChainIO()
        elif data_format == 'mongo':
            self.saver = project_io.MongoSideChainIO()

        self.loader = project_io.ChemDrawSideChainIO()
        self.id_generator = id_iterator

    def initialize(self):

        data = []
        for sidechain in self.loader.load():
            _id = self.id_generator.get_next()
            binary = sidechain.ToBinary()
            Chem.Kekulize(sidechain)
            data.append({'_id': _id,
                         'type': 'sidechain',
                         'binary': binary,
                         'kekule': Chem.MolToSmiles(sidechain, kekuleSmiles=True),
                         'connection': 'methyl',
                         'shared_id': _id})

        self.saver.save(data)
        self.id_generator.save()


class MonomerInitializer(IDataInitializer):

    def __init__(self, data_format, id_iterator, index_iterator):

        if data_format == 'json':
            self.saver = project_io.JsonMonomerIO()
        elif data_format == 'mongo':
            self.saver = project_io.MongoMonomerIO()

        self.loader = project_io.ChemDrawMonomerIO()
        self.id_generator = id_iterator
        self.index_generator = index_iterator

    def initialize(self):

        data = []
        for monomer in self.loader.load():
            binary = monomer.ToBinary()
            Chem.Kekulize(monomer)
            data.append({'_id': self.id_generator.get_next().upper(),
                         'type': 'monomer',
                         'binary': binary,
                         'kekule': Chem.MolToSmiles(monomer, kekuleSmiles=True),
                         'index': self.index_generator.get_next(),
                         'required': False,
                         'backbone': 'alpha',
                         'sidechain': None})

        self.saver.save(data)
        self.index_generator.save()


class IDInitializer(IDataInitializer):

    def __init__(self, data_format):

        if data_format == 'json':
            self.saver = project_io.JsonIDIO()
        elif data_format == 'mongo':
            self.saver = project_io.MongoIDIO()

    def initialize(self):

        self.saver.save({'_id': 'id',
                         'count': 0,
                         'prefix': ''})


class IndexInitializer(IDataInitializer):

    def __init__(self, data_format):
        if data_format == 'json':
            self.saver = project_io.JsonIndexIO()
        elif data_format == 'mongo':
            self.saver = project_io.MongoIndexIO()

    def initialize(self):
        self.saver.save({'_id': 'index',
                         'index': 1})
