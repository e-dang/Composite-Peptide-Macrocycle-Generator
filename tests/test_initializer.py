from cpmg.initializer import CPMGInitializer
from cpmg.exporters import RegioSQMExporter
import cpmg.repository as repo


def test_initialize_mol_only(partial_initialized_repository):
    assert len(list(repo.create_backbone_repository().load())) == 3
    assert len(list(repo.create_connection_repository().load())) == 2
    assert len(list(repo.create_template_repository().load())) == 3
    assert len(list(repo.create_sidechain_repository().load())) == 3
    assert len(list(repo.create_monomer_repository().load())) == 4


def test_initialize_predictions_only(partial_initialized_repository):
    exporter = RegioSQMExporter()
    exporter.export_regiosqm_smiles_file()

    initializer = CPMGInitializer()
    initializer.initialize_predictions_only()

    assert len(list(repo.create_backbone_repository().load())) == 3
    assert len(list(repo.create_connection_repository().load())) == 2
    assert len(list(repo.create_template_repository().load())) == 3
    assert len(list(repo.create_sidechain_repository().load())) == 3
    assert len(list(repo.create_monomer_repository().load())) == 4
    assert len(list(repo.create_regiosqm_repository().load())) == 3
    assert len(list(repo.create_pka_repository().load())) == 3


def test_initialize(initialized_repository):
    assert len(list(repo.create_backbone_repository().load())) == 3
    assert len(list(repo.create_connection_repository().load())) == 2
    assert len(list(repo.create_template_repository().load())) == 3
    assert len(list(repo.create_sidechain_repository().load())) == 3
    assert len(list(repo.create_monomer_repository().load())) == 4
    assert len(list(repo.create_regiosqm_repository().load())) == 3
    assert len(list(repo.create_pka_repository().load())) == 3
