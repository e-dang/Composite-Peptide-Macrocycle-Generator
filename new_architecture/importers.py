import new_architecture.repository.repository as repo
import new_architecture.models as models
import macrocycles.config as config
import os
import json
import glob
from rdkit import Chem
import uuid
from collections import namedtuple


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
        # data = [models.Backbone.from_mol(Chem.MolFromSmiles(backbone['kekule']))
        #         for backbone in self.loader.load(self.saver.CATEGORY)]

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


'''
sidechains should be imported from json format with following information : {'smiles': XXX, 'connection': XXX} where both smiles and connection are smiles strings.
monomers should be imported from json format with following information : {'smiles': XXX, 'connection': XXX, 'backbone': XXX, 'is_proline': XXX}
'''


class SidechainImporter:
    def __init__(self, loader):
        self.loader = loader
        self.saver = repo.create_sidechain_repository()
        self.connections = self._hash_connections(repo.create_connection_repository().load())
        self._validate_connections()

    def import_data(self):
        data = []
        for sidechain in self.loader.load(self.saver.CATEGORY):
            sidechain = self._match_connection(sidechain)
            sidechain['binary'] = Chem.MolFromSmiles(sidechain['kekule']).ToBinary()
            sidechain['shared_id'] = str(uuid.uuid4())
            data.append(models.Sidechain.from_dict(sidechain))

        return self.saver.save(data)

    def _hash_connections(self, connections):
        connection_dict = {}
        for connection in connections:
            connection_dict[connection.kekule] = connection._id

        return connection_dict

    def _validate_connections(self):
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
        self.backbones = self._hash_backbones(repo.create_backbone_repository().load())
        self._validate_backbones()

    def import_data(self):
        mock_sidechain = namedtuple('sidechain', 'shared_id connection')

        data = []
        for monomer in self.loader.load(self.saver.CATEGORY):
            data.append(models.Monomer.from_mol(Chem.MolFromSmiles(
                monomer['kekule']), self._match_backbone(monomer['backbone']), mock_sidechain(None, None), True))

        return self.saver.save(data)

    def _hash_backbones(self, backbones):
        backbone_dict = {}
        for backbone in backbones:
            backbone_dict[backbone.kekule] = backbone

        return backbone_dict

    def _validate_backbones(self):
        if len(self.backbones) == 0:
            raise RuntimeError(
                'No backbone molecules were found in the repository! Backbones must be imported before '
                'monomers can be imported.')

    def _match_backbone(self, backbone):
        try:
            return self.backbones[backbone]
        except KeyError:
            raise KeyError(f'Unrecognized backbone specified in imported monomers: {backbone}')
