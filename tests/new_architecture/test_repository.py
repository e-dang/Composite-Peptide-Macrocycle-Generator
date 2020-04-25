import pytest

import new_architecture.models as models
import new_architecture.repository.repository as repo
from tests.new_architecture.data.mols import *
from tests.new_architecture.test_hdf5 import filepath, initialize_repo


@pytest.fixture()
def repository_patch(monkeypatch, initialize_repo):
    yield repo.HDF5


@pytest.mark.parametrize('impl,expected_result', [(repo.HDF5, repo.HDF5Repository())])
def test_repository_impl_from_string(impl, expected_result):
    assert(type(repo.repository_impl_from_string(impl)) == type(expected_result))


def test_repository_impl_from_string_fail():
    with pytest.raises(ValueError):
        repo.repository_impl_from_string('dne')


def test_sidechain_repository(repository_patch):
    sc_repo = repo.create_sidechain_repository(repository_patch)
    _ids = sc_repo.save(list(map(models.Sidechain.from_dict, [TEST_SIDECHAIN_1, TEST_SIDECHAIN_2, TEST_SIDECHAIN_3])))
    data = sc_repo.load(_ids)
    data = list(map(lambda x: x.to_dict(), data))
    assert(len(data) == 3)
    assert(TEST_SIDECHAIN_1 in data)
    assert(TEST_SIDECHAIN_2 in data)
    assert(TEST_SIDECHAIN_3 in data)


def test_sidechain_repository_fail(repository_patch):
    sc_repo = repo.create_sidechain_repository(repository_patch)
    with pytest.raises(TypeError):
        _ids = sc_repo.save(['dne'])


def test_monomer_repository(repository_patch):
    monomer_repo = repo.create_monomer_repository(repository_patch)
    _ids = monomer_repo.save(list(map(models.Monomer.from_dict, [TEST_MONOMER_1, TEST_MONOMER_2, TEST_MONOMER_3])))
    data = monomer_repo.load(_ids)
    data = list(map(lambda x: x.to_dict(), data))
    assert(len(data) == 3)
    assert(TEST_MONOMER_1 in data)
    assert(TEST_MONOMER_2 in data)
    assert(TEST_MONOMER_3 in data)
