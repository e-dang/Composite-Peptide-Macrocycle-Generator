import os

import pytest

import new_architecture.importers as importers
import new_architecture.repository.repository as repo
from new_architecture.initializer import CPMGInitializer
from tests.new_architecture.test_hdf5 import filepath
from tests.new_architecture.test_importers import import_path_patch


@pytest.fixture
def repository_patch(monkeypatch, filepath, import_path_patch):
    monkeypatch.setattr(repo.create_backbone_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_connection_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_template_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_sidechain_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_monomer_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_peptide_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_template_peptide_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_repository_initializer, '__defaults__', (repo.HDF5,))
    initializer = CPMGInitializer()
    for importer in initializer.importers:
        importer.loader.search_dir = os.path.join(import_path_patch, 'imports')
    yield initializer


def test_initializer(repository_patch):
    initializer = repository_patch
    initializer.initialize()

    assert(len(list(repo.create_backbone_repository().load())) == 3)
    assert(len(list(repo.create_connection_repository().load())) == 2)
    assert(len(list(repo.create_template_repository().load())) == 3)
    assert(len(list(repo.create_sidechain_repository().load())) == 2)
    assert(len(list(repo.create_monomer_repository().load())) == 2)
