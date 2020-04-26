import glob
import json
import os
import uuid
from collections import namedtuple

from rdkit import Chem

import macrocycles.config as config
import new_architecture.models as models
import new_architecture.repository.repository as repo


class JsonImporter:
    def __init__(self):
        self.search_dir = os.path.join(config.DATA_DIR, 'imports')

    def load(self, data_type):

        for filepath in self._assemble_filepaths(data_type):
            with open(filepath, 'r') as file:
                for doc in json.load(file):
                    yield doc

    def _assemble_filepaths(self, data_type):
        return glob.glob(os.path.join(self.search_dir, data_type + '*.json'))


class ConnectionImporter:
    def __init__(self, loader):
        self.loader = loader
        self.saver = repo.create_connection_repository()

    def import_data(self):
        data = [models.Connection.from_mol(Chem.MolFromSmiles(connection['kekule']))
                for connection in self.loader.load(self.saver.CATEGORY)]
        return self.saver.save(data)


class BackboneImporter:
    def __init__(self, loader):
        self.loader = loader
        self.saver = repo.create_backbone_repository()

    def import_data(self):
        data = []
        for backbone in self.loader.load(self.saver.CATEGORY):
            backbone['binary'] = Chem.MolFromSmiles(backbone['mapped_kekule']).ToBinary()
            data.append(models.Backbone.from_dict(backbone))

        return self.saver.save(data)


class TemplateImporter:
    def __init__(self, loader):
        self.loader = loader
        self.saver = repo.create_template_repository()

    def import_data(self):
        data = [models.Template.from_mol(Chem.MolFromSmiles(template['kekule']))
                for template in self.loader.load(self.saver.CATEGORY)]
        return self.saver.save(data)


class SidechainImporter:
    def __init__(self, loader):
        self.loader = loader
        self.saver = repo.create_sidechain_repository()

    def import_data(self):
        self._load_connections()
        self._check_connections()

        data = []
        for sidechain in self.loader.load(self.saver.CATEGORY):
            sidechain = self._match_connection(sidechain)
            sidechain['binary'] = Chem.MolFromSmiles(sidechain['kekule']).ToBinary()
            sidechain['shared_id'] = str(uuid.uuid4())
            data.append(models.Sidechain.from_dict(sidechain))

        return self.saver.save(data)

    def _load_connections(self):
        self.connections = {}
        for connection in repo.create_connection_repository().load():
            self.connections[connection.kekule] = connection._id

    def _check_connections(self):
        if len(self.connections) == 0:
            raise RuntimeError(
                'No connection molecules were found in the repository! Connections must be imported before sidechains '
                'can be imported.')

    def _match_connection(self, sidechain):
        try:
            sidechain['connection'] = self.connections[sidechain['connection']]
        except KeyError:
            raise KeyError(f'Unrecognized connection specified in imported sidechain: {sidechain}')
        else:
            return sidechain


class MonomerImporter:

    def __init__(self, loader):
        self.loader = loader
        self.saver = repo.create_monomer_repository()

    def import_data(self):
        self._load_backbones()
        self._check_backbones()
        mock_sidechain = namedtuple('sidechain', 'shared_id connection')

        data = []
        for monomer in self.loader.load(self.saver.CATEGORY):
            data.append(models.Monomer.from_mol(Chem.MolFromSmiles(
                monomer['kekule']), self._match_backbone(monomer['backbone']), mock_sidechain(None, None), True))

        return self.saver.save(data)

    def _load_backbones(self):
        self.backbones = {}
        for backbone in repo.create_backbone_repository().load():
            self.backbones[backbone.kekule] = backbone

    def _check_backbones(self):
        if len(self.backbones) == 0:
            raise RuntimeError(
                'No backbone molecules were found in the repository! Backbones must be imported before '
                'monomers can be imported.')

    def _match_backbone(self, backbone):
        try:
            return self.backbones[backbone]
        except KeyError:
            raise KeyError(f'Unrecognized backbone specified in imported monomers: {backbone}')


def create_importers(data_format_importer=JsonImporter()):
    return [ConnectionImporter(data_format_importer),
            BackboneImporter(data_format_importer),
            TemplateImporter(data_format_importer),
            SidechainImporter(data_format_importer),
            MonomerImporter(data_format_importer)]
