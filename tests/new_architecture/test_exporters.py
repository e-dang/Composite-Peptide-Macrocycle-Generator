import os
import pytest
import new_architecture.exporters as exporters
import macrocycles.config as config
from rdkit import Chem
from rdkit.Chem import AllChem
from new_architecture.io_formats import json_load, text_load
from tests.new_architecture.test_repository import filepath, initialize_repo
from tests.new_architecture.test_initializer import repository_patch, import_path_patch


@pytest.fixture
def export_path_patch(monkeypatch, repository_patch):
    initializer = repository_patch
    initializer.initialize()
    export_filepath = os.path.join(config.PROJECT_DIR, 'tests', 'new_architecture', 'data', 'regiosqm_smiles.smiles')
    sidechain_import_path = os.path.join(config.PROJECT_DIR, 'tests',
                                         'new_architecture', 'data', 'imports', 'sidechains.json')
    monomer_import_path = os.path.join(config.PROJECT_DIR, 'tests', 'new_architecture',
                                       'data', 'imports', 'monomers.json')
    monkeypatch.setattr(exporters.config, 'REGIOSQM_SMILES_FILEPATH', export_filepath)
    yield export_filepath, sidechain_import_path, monomer_import_path
    os.remove(export_filepath)


def is_required(monomer):
    return bool(AllChem.CalcNumAromaticRings(Chem.MolFromSmiles(monomer['kekule'])))


def test_reqiosqm_exporter(export_path_patch):
    export_filepath, sidechain_filepath, monomer_filepath = export_path_patch

    exporter = exporters.RegioSQMExporter()

    exporter.export_regiosqm_smiles_file()

    sidechains = list(filter(lambda x: x['connection'] != 'C', json_load(sidechain_filepath)))
    monomers = list(filter(is_required, json_load(monomer_filepath)))
    kekules = [doc['kekule'] for doc in sidechains + monomers]

    with open(export_filepath, 'r') as file:
        lines = list(file.readlines())
        assert(len(lines) == len(sidechains) + len(monomers))
        for line in lines:
            kekule = line.split(' ')[-1].strip('\n')
            assert(kekule in kekules)
            kekules.remove(kekule)
