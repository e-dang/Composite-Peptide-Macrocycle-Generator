import os

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import conftest
import cpmg.config as config
import cpmg.exporters as exporters
from cpmg.io_formats import load_json


@pytest.fixture()
def import_paths(monkeypatch, hdf5_repository_with_imports):
    sidechain_import_path = os.path.join(config.IMPORT_DIR, 'sidechains.json')
    monomer_import_path = os.path.join(config.IMPORT_DIR, 'monomers.json')
    yield sidechain_import_path, monomer_import_path


def is_required(monomer):
    return bool(AllChem.CalcNumAromaticRings(Chem.MolFromSmiles(monomer['kekule'])))


def test_reqiosqm_exporter(import_paths):
    sidechain_filepath, monomer_filepath = import_paths

    exporter = exporters.RegioSQMExporter()

    exporter.export_regiosqm_smiles_file()

    sidechains = list(filter(lambda x: x['connection'] == 'C', load_json(sidechain_filepath)))
    monomers = list(filter(is_required, load_json(monomer_filepath)))
    kekules = [doc['kekule'] for doc in sidechains + monomers]

    with open(config.REGIOSQM_SMILES_FILEPATH, 'r') as file:
        lines = list(file.readlines())
        assert(len(lines) == len(sidechains) + len(monomers))
        for line in lines:
            kekule = line.split(' ')[-1].strip('\n')
            assert(kekule in kekules)
            kekules.remove(kekule)
