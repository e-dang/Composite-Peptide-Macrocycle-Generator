import os

import pytest
from rdkit.Chem import AllChem

import macrocycles.config as config
import new_architecture.importers as importers
import new_architecture.repository.repository as repo
from new_architecture.models import PROLINE_N_TERM
from tests.new_architecture.test_hdf5 import filepath, initialize_repo
from tests.new_architecture.test_repository import repository_patch
import macrocycles.exceptions as exceptions


@pytest.fixture
def import_path_patch(monkeypatch):
    filepath = os.path.join(config.PROJECT_DIR, 'tests', 'new_architecture', 'data')
    monkeypatch.setattr(importers.config, 'DATA_DIR', filepath)
    yield filepath


@pytest.fixture
def json_importer(import_path_patch):
    yield importers.JsonImporter()


@pytest.fixture
def independent_importers(json_importer, repository_patch):
    connection_importer = importers.ConnectionImporter(json_importer)
    backbone_importer = importers.BackboneImporter(json_importer)
    template_importer = importers.TemplateImporter(json_importer)
    connection_importer.import_data()
    backbone_importer.import_data()
    template_importer.import_data()


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


def test_backbone_importer(json_importer, repository_patch):
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


def test_backbone_importer_fail(monkeypatch, json_importer, repository_patch):
    monkeypatch.setattr(importers.repo.BackboneRepository, 'CATEGORY', 'invalid_backbones')
    backbone_importer = importers.BackboneImporter(json_importer)
    with pytest.raises(exceptions.InvalidMolecule):
        ids = backbone_importer.import_data()


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


def test_sidechain_importer(json_importer, independent_importers):
    sidechain_importer = importers.SidechainImporter(json_importer)
    ids = sidechain_importer.import_data()

    sidechain_repo = repo.create_sidechain_repository()
    connection_repo = repo.create_connection_repository()
    sidechain_data = list(sidechain_repo.load(ids))
    sidechain_docs = json_importer.load(sidechain_importer.saver.CATEGORY)
    connection_data = list(connection_repo.load())
    kekules = [doc['kekule'] for doc in sidechain_docs]
    connection_ids = [mol._id for mol in connection_data]

    assert(len(sidechain_data) == 2)
    for mol in sidechain_data:
        assert(mol._id != None)
        assert(mol.shared_id != None)
        assert(mol.kekule in kekules)
        assert(mol.connection in connection_ids)
        kekules.remove(mol.kekule)


def test_monomer_importer(json_importer, independent_importers):
    monomer_importer = importers.MonomerImporter(json_importer)
    ids = monomer_importer.import_data()

    monomer_repo = repo.create_monomer_repository()
    backbone_repo = repo.create_backbone_repository()
    monomer_data = list(monomer_repo.load(ids))
    monomer_docs = json_importer.load(monomer_importer.saver.CATEGORY)
    backbone_data = list(backbone_repo.load())
    kekules = [doc['kekule'] for doc in monomer_docs]
    backbone_ids = [mol._id for mol in backbone_data]

    assert(len(monomer_data) == 2)
    for mol in monomer_data:
        rdkit_mol = mol.mol
        assert(mol._id != None)
        assert(mol.required == bool(AllChem.CalcNumAromaticRings(rdkit_mol)))
        assert(mol.backbone in backbone_ids)
        assert(mol.sidechain is None)
        assert(mol.connection is None)
        assert(mol.proline == bool(AllChem.CalcNumAliphaticRings(
            rdkit_mol) and rdkit_mol.HasSubstructMatch(PROLINE_N_TERM)))
        assert(mol.imported == True)


def test_create_importers(json_importer, repository_patch):
    instances = importers.create_importers()
    types = list(map(type, instances))
    assert(len(instances) == 5)
    assert(importers.ConnectionImporter in types)
    assert(importers.BackboneImporter in types)
    assert(importers.TemplateImporter in types)
    assert(importers.SidechainImporter in types)
    assert(importers.MonomerImporter in types)
