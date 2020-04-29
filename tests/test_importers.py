import os

import pytest
from rdkit.Chem import AllChem

import conftest
import cpmg.config as config
import cpmg.exceptions as exceptions
import cpmg.importers as importers
import cpmg.repository as repo
import cpmg.exporters as exporters
from cpmg.models import PROLINE_N_TERM
from data.mols import *


@pytest.fixture()
def json_importer(hdf5_repository):
    yield importers.JsonImporter()


@pytest.fixture()
def independent_importers(json_importer):
    connection_importer = importers.ConnectionImporter(json_importer)
    backbone_importer = importers.BackboneImporter(json_importer)
    template_importer = importers.TemplateImporter(json_importer)
    connection_importer.import_data()
    backbone_importer.import_data()
    template_importer.import_data()


@pytest.fixture()
def mol_importers(json_importer, independent_importers):
    sidechain_importer = importers.SidechainImporter(json_importer)
    monomer_importer = importers.MonomerImporter(json_importer)
    sidechain_importer.import_data()
    monomer_importer.import_data()


def test_json_importer_assemble_filepaths(json_importer):
    filepaths = json_importer._assemble_filepaths('test_type')

    assert(len(filepaths) == 2)
    assert(filepaths[0] == os.path.join(config.IMPORT_DIR, 'test_type1.json'))
    assert(filepaths[1] == os.path.join(config.IMPORT_DIR, 'test_type2.json'))


def test_json_importer_load(json_importer):
    docs = list(json_importer.load('test_type'))

    assert(len(docs) == 2)
    assert(docs[0] == {'A': 1})
    assert(docs[1] == {'B': 2})


def test_connection_importer(json_importer):
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


def test_backbone_importer(json_importer):
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


def test_backbone_importer_fail(monkeypatch, json_importer):
    monkeypatch.setattr(importers.repo.BackboneRepository, 'CATEGORY', 'invalid_backbones')
    backbone_importer = importers.BackboneImporter(json_importer)
    with pytest.raises(exceptions.InvalidMolecule):
        ids = backbone_importer.import_data()


def test_template_importer(json_importer):
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
    connections = [mol.kekule for mol in connection_data]

    assert(len(sidechain_data) == 3)
    for mol in sidechain_data:
        assert(mol._id != None)
        assert(mol.shared_id != None)
        assert(mol.kekule in kekules)
        assert(mol.connection in connections)
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
    backbones = [mol.to_reduced_dict() for mol in backbone_data]

    assert(len(monomer_data) == 4)
    for mol in monomer_data:
        rdkit_mol = mol.mol
        assert(mol._id != None)
        assert(mol.required == bool(AllChem.CalcNumAromaticRings(rdkit_mol)))
        assert(mol.backbone in backbones)
        assert(mol.sidechain is None)
        assert(mol.connection is None)
        assert(mol.proline == bool(AllChem.CalcNumAliphaticRings(
            rdkit_mol) and rdkit_mol.HasSubstructMatch(PROLINE_N_TERM)))
        assert(mol.imported == True)
        assert(mol.kekule in kekules)
        kekules.remove(mol.kekule)


def test_regiosqm_prediction_importer(mol_importers):
    exporter = exporters.RegioSQMExporter()
    exporter.export_regiosqm_smiles_file()

    regiosqm_importer = importers.RegioSQMPredictionImporter()

    ids = regiosqm_importer.import_data()

    regiosqm_repo = repo.create_regiosqm_repository()
    regiosqm_data = list(regiosqm_repo.load(ids))

    assert(len(regiosqm_data) == 3)
    for prediction in regiosqm_data:
        assert(prediction.solvent == 'nitromethane')
        assert(prediction.cutoff == 3.0)
        assert(prediction.reacting_mol in ('CC1=CC=C[NH]1', 'CC1=CC=C(O)C=C1',
                                           'O=C(O)[C@@H]1C[C@H](OC2=CC=NC3=C2SC=C3)CN1'))
        assert(prediction.predictions in ([3, 6], [2, 3, 4], [15]))


def test_pka_prediction_importer(mol_importers):
    pka_importer = importers.pKaPredictionImporter()

    ids = pka_importer.import_data()

    pka_repo = repo.create_pka_repository()
    pka_data = list(pka_repo.load(ids))

    test_data = [TEST_PKA_PREDICTION_1, TEST_PKA_PREDICTION_2, TEST_PKA_PREDICTION_3]
    assert(len(pka_data) == 3)
    for prediction in pka_data:
        print(prediction.to_dict())

        assert(prediction.to_dict() in test_data)


def test_create_importers(hdf5_repository):
    instances = importers.create_importers()
    types = list(map(type, instances))
    assert(len(instances) == 5)
    assert(importers.ConnectionImporter in types)
    assert(importers.BackboneImporter in types)
    assert(importers.TemplateImporter in types)
    assert(importers.SidechainImporter in types)
    assert(importers.MonomerImporter in types)
