import os

import pytest

import macrocycles.config as config
import new_architecture.importers as importers
import new_architecture.repository.repository as repo
from tests.new_architecture.test_hdf5 import filepath, initialize_repo
from tests.new_architecture.test_repository import repository_patch


@pytest.fixture
def import_path_patch(monkeypatch):
    filepath = os.path.join(config.PROJECT_DIR, 'tests', 'new_architecture', 'data')
    monkeypatch.setattr(importers.config, 'DATA_DIR', filepath)
    yield filepath


@pytest.fixture
def json_importer(import_path_patch):
    yield importers.JsonImporter()


def test_json_importer_assemble_filepaths(json_importer):
    filepaths = json_importer._assemble_filepaths('test_type')

    assert(len(filepaths) == 2)
    assert(filepaths[0] == os.path.join(json_importer.search_dir, 'test_type1.json'))
    assert(filepaths[1] == os.path.join(json_importer.search_dir, 'test_type2.json'))


def test_json_importer_load(json_importer):
    docs = list(json_importer.load('test_type'))

    assert(len(docs) == 2)
    assert(docs[0] == {'A': 1})
    assert(docs[1] == {'B': 2})


def test_connection_importer(json_importer, repository_patch):
    connection_importer = importers.ConnectionImporter(json_importer)
    ids = connection_importer.import_data()

    connection_repo = repo.create_connection_repository()
    data = list(connection_repo.load(ids))
    docs = json_importer.load(connection_importer.saver.CATEGORY)
    kekules = [doc['kekule'] for doc in docs]

    assert(len(data) == 2)
    for mol in data:
        assert(mol._id != None)
        assert(mol.kekule in kekules)
        kekules.remove(mol.kekule)


def test_bacbone_importer(json_importer, repository_patch):
    backbone_importer = importers.BackboneImporter(json_importer)
    ids = backbone_importer.import_data()

    backbone_repo = repo.create_backbone_repository()
    data = list(backbone_repo.load(ids))
    docs = json_importer.load(backbone_importer.saver.CATEGORY)
    kekules = [doc['kekule'] for doc in docs]

    assert(len(data) == 3)
    for mol in data:
        assert(mol._id != None)
        assert(mol.kekule in kekules)
        kekules.remove(mol.kekule)


def test_template_importer(json_importer, repository_patch):
    template_importer = importers.TemplateImporter(json_importer)
    ids = template_importer.import_data()

    template_repo = repo.create_template_repository()
    data = list(template_repo.load(ids))
    docs = json_importer.load(template_importer.saver.CATEGORY)
    kekules = [doc['kekule'] for doc in docs]

    assert(len(data) == 3)
    for mol in data:
        assert(mol._id != None)
        assert(mol.kekule in kekules)
        kekules.remove(mol.kekule)
