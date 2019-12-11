from abc import ABC, abstractmethod

from rdkit import Chem

import iterators
import project_io


class IDataInitializer(ABC):

    @abstractmethod
    def initialize(self):
        pass


class DataInitializer(IDataInitializer):

    def initialize(self):

        self.phase_1()
        self.phase_2()

    def phase_1(self):

        self.record_initializers = [IDInitializer(), IndexInitializer()]
        self.select_initialization(self.record_initializers)

    def phase_2(self):

        self.mol_initializers = [SideChainInitializer(), MonomerInitializer()]
        self.select_initialization(self.mol_initializers)

    def select_initialization(self, initializers):

        for initializer in initializers:
            initializer.initialize()


class SideChainInitializer(IDataInitializer):

    def __init__(self):
        self.loader = project_io.ChemDrawSideChainIO()
        self.saver = project_io.JsonSideChainIO()
        self.id_generator = iterators.IDIterator()

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


class MonomerInitializer(IDataInitializer):

    def __init__(self):

        self.loader = project_io.ChemDrawMonomerIO()
        self.saver = project_io.JsonMonomerIO()
        self.id_generator = iterators.IDIterator()
        self.index_generator = iterators.IndexIterator()

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


class IDInitializer(IDataInitializer):

    def __init__(self):

        self.saver = project_io.JsonIDIO()

    def initialize(self):

        self.saver.save([{'_id': 'id',
                          'count': 0,
                          'prefix': ''}])


class IndexInitializer(IDataInitializer):

    def __init__(self):

        self.saver = project_io.JsonIndexIO()

    def initialize(self):
        self.saver.save([{'_id': 'index',
                          'index': 1}])
