import os
import pytest
import cpmg.hdf5 as hdf5
import cpmg.config as config
from cpmg.initializer import CPMGInitializer


@pytest.fixture(autouse=True)
def patched_configs(tmpdir, monkeypatch):
    monkeypatch.setattr('cpmg.config.DATA_DIR', os.path.join(config.TEST_DIR, 'data'))
    monkeypatch.setattr('cpmg.config.IMPORT_DIR', os.path.join(config.TEST_DIR, 'data', 'imports'))
    monkeypatch.setattr('cpmg.config.DATA_FORMAT', 'hdf5')
    monkeypatch.setattr('cpmg.config.HDF5_FILEPATH', os.path.join(str(tmpdir), 'test_hdf5_repo.hdf5'))
    monkeypatch.setattr('cpmg.config.REGIOSQM_SMILES_FILEPATH', os.path.join(str(tmpdir), 'regiosqm_smiles.txt'))


@pytest.fixture()
def hdf5_repository():
    initializer = hdf5.HDF5Initializer()
    initializer.initialize()
    yield hdf5.HDF5Repository()


@pytest.fixture()
def hdf5_repository_with_imports(hdf5_repository):
    initializer = CPMGInitializer()
    initializer.initialize()
