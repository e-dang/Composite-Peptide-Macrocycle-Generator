import pytest

import new_architecture.models as models
import new_architecture.repository.repository as repo
import tests.new_architecture.utils as utils
from tests.new_architecture.data.mols import *
from tests.new_architecture.test_hdf5 import filepath, initialize_repo


@pytest.fixture()
def repository_patch(monkeypatch, initialize_repo):
    monkeypatch.setattr(repo.create_backbone_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_connection_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_template_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_sidechain_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_monomer_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_peptide_repository, '__defaults__', (repo.HDF5,))
    monkeypatch.setattr(repo.create_template_peptide_repository, '__defaults__', (repo.HDF5,))


@pytest.mark.parametrize('impl,expected_result', [(repo.HDF5, repo.HDF5Repository())])
def test_repository_impl_from_string(impl, expected_result):

    assert(type(repo.repository_impl_from_string(impl)) == type(expected_result))


def test_repository_impl_from_string_fail():
    with pytest.raises(ValueError):
        repo.repository_impl_from_string('dne')


def test_backbone_repository(repository_patch):
    backbone_repo = repo.create_backbone_repository()
    _ids = backbone_repo.save(list(map(models.Backbone.from_dict, [TEST_BACKBONE_1, TEST_BACKBONE_2, TEST_BACKBONE_3])))
    data = backbone_repo.load(_ids)
    data = utils.mols_to_dict(data)

    assert(len(data) == 3)
    assert(TEST_BACKBONE_1 in data)
    assert(TEST_BACKBONE_2 in data)
    assert(TEST_BACKBONE_3 in data)


def test_backbone_repository_fail(repository_patch):
    backbone_repo = repo.create_backbone_repository()

    _ids = backbone_repo.save(['dne'])

    assert(len(backbone_repo.failed_instances) == 1)


def test_connection_repository(repository_patch):
    connection_repo = repo.create_connection_repository()
    _ids = connection_repo.save(list(map(models.Connection.from_dict, [TEST_CONNECTION_1])))
    data = connection_repo.load(_ids)
    data = utils.mols_to_dict(data)

    assert(len(data) == 1)
    assert(TEST_CONNECTION_1 in data)


def test_connection_repository_fail(repository_patch):
    connection_repo = repo.create_connection_repository()

    _ids = connection_repo.save(['dne'])

    assert(len(connection_repo.failed_instances) == 1)


def test_template_repository(repository_patch):
    template_repo = repo.create_template_repository()
    _ids = template_repo.save(list(map(models.Template.from_dict, [TEST_TEMPLATE_1, TEST_TEMPLATE_2, TEST_TEMPLATE_3])))
    data = template_repo.load(_ids)
    data = utils.mols_to_dict(data)

    assert(len(data) == 3)
    assert(TEST_TEMPLATE_1 in data)
    assert(TEST_TEMPLATE_2 in data)
    assert(TEST_TEMPLATE_3 in data)


def test_template_repository_fail(repository_patch):
    template_repo = repo.create_template_repository()

    _ids = template_repo.save(['dne'])

    assert(len(template_repo.failed_instances) == 1)


def test_sidechain_repository(repository_patch):
    sc_repo = repo.create_sidechain_repository()
    _ids = sc_repo.save(list(map(models.Sidechain.from_dict, [TEST_SIDECHAIN_1, TEST_SIDECHAIN_2, TEST_SIDECHAIN_3])))
    data = sc_repo.load(_ids)
    data = utils.mols_to_dict(data)

    assert(len(data) == 3)
    assert(TEST_SIDECHAIN_1 in data)
    assert(TEST_SIDECHAIN_2 in data)
    assert(TEST_SIDECHAIN_3 in data)


def test_sidechain_repository_fail(repository_patch):
    sc_repo = repo.create_sidechain_repository()

    _ids = sc_repo.save(['dne'])

    assert(len(sc_repo.failed_instances) == 1)


def test_monomer_repository(repository_patch):
    monomers = [TEST_MONOMER_1, TEST_MONOMER_2, TEST_MONOMER_3]
    monomer_repo = repo.create_monomer_repository()
    _ids = monomer_repo.save(list(map(models.Monomer.from_dict, monomers)))
    data = monomer_repo.load(_ids)
    data = utils.mols_to_dict(data)

    for i, monomer in enumerate(monomers):
        monomer['index'] = i

    monomers.sort(key=lambda x: x['index'])
    data.sort(key=lambda x: x['index'])

    assert(len(data) == len(monomers))
    for doc, monomer in zip(data, monomers):
        assert(doc == monomer)


def test_monomer_repository_fail(repository_patch):
    monomer_repo = repo.create_monomer_repository()

    _ids = monomer_repo.save(['dne'])

    assert(len(monomer_repo.failed_instances) == 1)


def test_peptide_repository(repository_patch):
    peptide_repo = repo.create_peptide_repository()
    _ids = peptide_repo.save(list(map(models.Peptide.from_dict, [TEST_PEPTIDE_1])))
    data = peptide_repo.load(_ids)
    data = utils.mols_to_dict(data)

    assert(len(data) == 1)
    assert(TEST_PEPTIDE_1 == data[0])


def test_peptide_repository_fail(repository_patch):
    peptide_repo = repo.create_peptide_repository()

    _ids = peptide_repo.save(['dne'])

    assert(len(peptide_repo.failed_instances) == 1)


def test_template_peptide_repository(repository_patch):
    template_peptide_repo = repo.create_template_peptide_repository()
    _ids = template_peptide_repo.save(list(map(models.TemplatePeptide.from_dict, [TEST_TEMPLATE_PEPTIDE_1])))
    data = template_peptide_repo.load(_ids)
    data = utils.mols_to_dict(data)

    assert(len(data) == 1)
    assert(TEST_TEMPLATE_PEPTIDE_1 == data[0])


def test_template_peptide_repository_fail(repository_patch):
    template_peptide_repo = repo.create_template_peptide_repository()

    _ids = template_peptide_repo.save(['dne'])

    assert(len(template_peptide_repo.failed_instances) == 1)
