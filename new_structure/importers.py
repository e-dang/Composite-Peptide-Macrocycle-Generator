from abc import ABC, abstractmethod

from rdkit import Chem

import iterators
import project_io


class IImporter(ABC):

    @abstractmethod
    def import_data(self):
        pass


class DataImporter(IImporter):

    def __init__(self, data_format):
        if data_format == 'json':
            sidechain_io = project_io.JsonSideChainIO()
            monomer_io = project_io.JsonMonomerIO()
            id_io = project_io.JsonIDIO()
            index_io = project_io.JsonIndexIO()
            regiosqm_io = project_io.JsonRegioSQMIO()
        elif data_format == 'mongo':
            sidechain_io = project_io.MongoSideChainIO()
            monomer_io = project_io.MongoMonomerIO()
            id_io = project_io.MongoIDIO()
            index_io = project_io.MongoIndexIO()
            regiosqm_io = project_io.MongoRegioSQMIO()

        id_iterator = iterators.IDIterator(id_io)
        index_iterator = iterators.IndexIterator(index_io)
        self.molecule_importers = [SideChainImporter(sidechain_io, id_iterator),
                                   MonomerImporter(monomer_io, id_iterator, index_iterator)]
        self.prediction_importers = [RegioSQMImporter(regiosqm_io)]

    def import_data(self):
        self.import_molecules()
        self.import_predictions()

    def import_molecules(self):
        for importer in self.molecule_importers:
            importer.import_data()

    def import_predictions(self):
        for importer in self.prediction_importers:
            importer.import_data()


class SideChainImporter(IImporter):

    def __init__(self, sidechain_io, id_iterator):

        self.loader = project_io.ChemDrawSideChainIO()
        self.saver = sidechain_io
        self.id_iterator = id_iterator

    def import_data(self):
        data = []
        for sidechain in self.loader.load():
            _id = self.id_iterator.get_next()
            binary = sidechain.ToBinary()
            Chem.Kekulize(sidechain)
            data.append({'_id': _id,
                         'type': 'sidechain',
                         'binary': binary,
                         'kekule': Chem.MolToSmiles(sidechain, kekuleSmiles=True),
                         'connection': 'methyl',
                         'shared_id': _id})

        self.saver.save(data)
        self.id_iterator.save()


class MonomerImporter(IImporter):

    def __init__(self, monomer_io, id_iterator, index_iterator):

        self.loader = project_io.ChemDrawMonomerIO()
        self.saver = monomer_io
        self.id_iterator = id_iterator
        self.index_iterator = index_iterator

    def import_data(self):

        data = []
        for monomer in self.loader.load():
            binary = monomer.ToBinary()
            Chem.Kekulize(monomer)
            data.append({'_id': self.id_iterator.get_next().upper(),
                         'type': 'monomer',
                         'binary': binary,
                         'kekule': Chem.MolToSmiles(monomer, kekuleSmiles=True),
                         'index': self.index_iterator.get_next(),
                         'required': False,
                         'backbone': 'alpha',
                         'sidechain': None})

        self.saver.save(data)
        self.index_iterator.save()


class RegioSQMImporter(IImporter):

    def __init__(self, regiosqm_io, solvent='nitromethane', cutoff=3):

        self.raw_regiosqm_io = project_io.RawRegioSQMIO()
        self.regiosqm_io = regiosqm_io
        self.solvent = solvent
        self.cutoff = cutoff

    def import_data(self):

        data = []
        for i, row in enumerate(self.raw_regiosqm_io.load()):
            if i % 3 == 0:
                _id = row[0]
            elif i % 3 == 1:
                predictions = [int(row[j]) for j in range(0, len(row), 2)]
            else:
                data.append({'_id': _id,
                             'type': 'regiosqm',
                             'predictions': predictions,
                             'solvent': self.solvent,
                             'cutoff': self.cutoff})

        self.regiosqm_io.save(data)


class pKaImporter(IImporter):

    def __init__(self):
        pass

    def import_data(self):
        pass
