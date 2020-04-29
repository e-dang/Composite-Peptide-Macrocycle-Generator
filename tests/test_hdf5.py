import os
from copy import copy

import h5py
import pytest
from rdkit import Chem

import conftest
import cpmg.config as config
import cpmg.hdf5 as hdf5
import cpmg.ranges as ranges
from data.mols import *

TEST_DICT = {'A': 1, 'B': True, 'C': 1.0, 'D': b'test_bin_string',
             'E': 'test_string', 'F': Chem.MolFromSmiles('CCC').ToBinary()}
TEST_LIST = [TEST_DICT, TEST_DICT, TEST_DICT]


@pytest.fixture
def sidechain_repo(hdf5_repository):
    group = 'sidechains'
    dataset_1 = [TEST_SIDECHAIN_1, TEST_SIDECHAIN_2]
    dataset_2 = [TEST_SIDECHAIN_3, TEST_SIDECHAIN_4]
    ids = hdf5_repository.save(group, dataset_1)
    ids.extend(hdf5_repository.save(group, dataset_2))
    yield hdf5_repository, ids, group, dataset_1, dataset_2


@pytest.fixture
def monomer_repo(hdf5_repository):
    group = 'monomers'
    dataset_1 = [TEST_MONOMER_1, TEST_MONOMER_2, TEST_MONOMER_3]
    dataset_2 = [TEST_MONOMER_4, TEST_MONOMER_5, TEST_MONOMER_6]
    ids = hdf5_repository.save(group, dataset_1)
    ids.extend(hdf5_repository.save(group, dataset_2))
    yield hdf5_repository, ids, group


def test_serialize_deserialize():
    serialized_data = hdf5.serialize(TEST_DICT)

    assert(isinstance(serialized_data, str))

    deserialized_data = hdf5.deserialize(serialized_data)

    assert(TEST_DICT == deserialized_data)


def test_serialize_deserialize_chunk():
    serialized_data = hdf5.serialize_chunk(TEST_LIST)

    for serialized_str in serialized_data:
        assert(isinstance(serialized_str, str))

    deserialized_data = hdf5.deserialize_chunk(serialized_data)

    for doc, deserialized_doc in zip(TEST_LIST, deserialized_data):
        assert(doc == deserialized_doc)


def test_hdf5_file_regular_open():
    file = hdf5.HDF5File()

    assert(os.path.exists(config.HDF5_FILEPATH))
    assert(file)
    assert(len(list(file.keys())) == 0)

    file.close()

    assert(not file)


def test_hdf5_file_context_manager():
    with hdf5.HDF5File() as file:
        assert(os.path.exists(config.HDF5_FILEPATH))
        assert(file)
        assert(len(list(file.keys())) == 0)

    assert(not file)


def test_hdf5_file_create_group():
    group_name = 'test'
    with hdf5.HDF5File() as file:
        group1 = file.create_group(group_name)
        group2 = file.create_group(group_name)

        assert(isinstance(group1, h5py.Group))
        assert(isinstance(group2, h5py.Group))
        assert(group1.name == group2.name)


def test_hdf5_initializer():
    initializer = hdf5.HDF5Initializer()
    initializer.initialize()
    with hdf5.HDF5File() as file:
        assert(sorted(list(file.keys())) == sorted(initializer.data_types))


def test_hdf5_repository_save_load_range_single(hdf5_repository):
    group = 'sidechains'
    hdf5_repository.save(group, TEST_SIDECHAIN_1)

    with hdf5.HDF5File() as file:
        assert(list(file[group].keys()) == ['0'])
        assert(len(file[group]['0']) == 1)

    _, data = zip(*list(hdf5_repository.load(group, ranges.WholeRange())))
    assert(len(data) == 1)
    assert(data[0] == TEST_SIDECHAIN_1)


def test_hdf5_repository_save_load_range_multi(sidechain_repo):
    repo, ids, group, dataset_1, dataset_2 = sidechain_repo
    test_dataset = sorted(dataset_1 + dataset_2, key=lambda x: x['kekule'])

    with hdf5.HDF5File() as file:
        assert(list(file[group].keys()) == ['0', '1'])
        assert(len(file[group]['0']) == 2)
        assert(len(file[group]['1']) == 2)

    _, data = zip(*list(repo.load(group, ranges.WholeRange())))
    data = sorted(list(data), key=lambda x: x['kekule'])
    assert(len(data) == 4)
    for doc, test_doc in zip(data, test_dataset):
        assert(doc == test_doc)


def test_hdf5_repository_save_load_discrete_chunk_multi(sidechain_repo):
    repo, ids, group, dataset_1, dataset_2 = sidechain_repo
    test_dataset = sorted([dataset_1[0], dataset_2[1]], key=lambda x: x['kekule'])

    with hdf5.HDF5File() as file:
        assert(list(file[group].keys()) == ['0', '1'])
        assert(len(file[group]['0']) == 2)
        assert(len(file[group]['1']) == 2)

    _, data = zip(*list(repo.load(group, ranges.DiscreteDataChunk([0, 3]))))
    data = sorted(list(data), key=lambda x: x['kekule'])
    assert(len(data) == 2)
    for doc, test_doc in zip(data, test_dataset):
        assert(doc == test_doc)


def test_hdf5_repository_save_load_ids_single(hdf5_repository):
    group = 'sidechains'
    ids = hdf5_repository.save(group, TEST_SIDECHAIN_1)

    with hdf5.HDF5File() as file:
        assert(list(file[group].keys()) == ['0'])
        assert(len(file[group]['0']) == 1)

    _, data = zip(*list(hdf5_repository.load(group, ids)))
    assert(len(data) == 1)
    assert(data[0] == TEST_SIDECHAIN_1)


def test_hdf5_repository_save_load_ids_multi(sidechain_repo):
    repo, ids, group, dataset_1, dataset_2 = sidechain_repo
    test_dataset = sorted(dataset_1 + dataset_2, key=lambda x: x['kekule'])

    with hdf5.HDF5File() as file:
        assert(list(file[group].keys()) == ['0', '1'])
        assert(len(file[group]['0']) == 2)
        assert(len(file[group]['1']) == 2)

    _, data = zip(*list(repo.load(group, ids)))
    data = sorted(list(data), key=lambda x: x['kekule'])
    assert(len(data) == 4)
    for doc, test_doc in zip(data, test_dataset):
        assert(doc == test_doc)


@pytest.mark.parametrize('dataset1,dataset2,expected_result,hdf5_repository', [
    ([TEST_MONOMER_1, TEST_MONOMER_2, TEST_MONOMER_3], [TEST_MONOMER_4, TEST_MONOMER_5, TEST_MONOMER_6], 6, ''),
    ([], [], 0, '')], indirect=['hdf5_repository'])
def test_hdf5_get_num_records(dataset1, dataset2, expected_result, hdf5_repository):
    ids = hdf5_repository.save('monomers', dataset1)
    ids.extend(hdf5_repository.save('monomers', dataset2))

    assert(hdf5_repository.get_num_records('monomers') == expected_result)


def test_hdf5_find(monomer_repo):
    repo, ids, group = monomer_repo

    key = [ids[0], ids[3], ids[4]]
    locations = repo.find(group, key)
    values = sorted(list(locations.values()))
    for value in values:
        value.sort()

    assert(len(key) == 0)
    assert(len(locations) == 2)
    assert(sorted(list(locations.keys())) == ['/' + group + '/0', '/' + group + '/1'])
    assert(values == [[0], [0, 1]])


def test_hdf5_find_fail(monomer_repo):
    repo, ids, group = monomer_repo

    key = ['dne']
    locations = repo.find(group, key)

    assert(len(key) == 1)
    assert(len(locations) == 0)


def test_hdf5_remove(monomer_repo):
    repo, ids, group = monomer_repo
    key = [ids[0], ids[3], ids[4]]

    assert(repo.remove(group, copy(key)))
    with hdf5.HDF5File() as file:
        assert(len(file[group]['0']) == 2)
        assert(len(file[group]['1']) == 1)

    locations = repo.find(group, key)
    assert(len(key) == 3)
    assert(len(locations) == 0)


def test_hdf5_move(monomer_repo):
    repo, ids, group = monomer_repo
    key = [ids[0], ids[3], ids[4]]
    dest_group = 'misc/monomers'

    # assert remove operation worked
    assert(repo.move(group, copy(key), dest_group))
    with hdf5.HDF5File() as file:
        assert(len(file[group]['0']) == 2)
        assert(len(file[group]['1']) == 1)
    locations = repo.find(group, key)
    assert(len(key) == 3)
    assert(len(locations) == 0)

    # assert copy operation worked
    locations = repo.find(dest_group, key)
    assert(len(key) == 0)
    assert(len(locations) == 1)


def test_hdf5_deactivate_records(monomer_repo):
    repo, ids, group = monomer_repo
    key = [ids[0], ids[3], ids[4]]
    dest_group = 'inactives/monomers'

    # assert remove operation worked
    assert(repo.deactivate_records(group, copy(key)))
    with hdf5.HDF5File() as file:
        assert(len(file[group]['0']) == 2)
        assert(len(file[group]['1']) == 1)
    locations = repo.find(group, key)
    assert(len(key) == 3)
    assert(len(locations) == 0)

    # assert copy operation worked
    locations = repo.find(dest_group, key)
    assert(len(key) == 0)
    assert(len(locations) == 1)


def test_hdf5_activate_records(monomer_repo):
    repo, ids, group = monomer_repo
    key = [ids[0], ids[3], ids[4]]
    dest_group = 'inactives/monomers'
    repo.deactivate_records(group, copy(key))

    # assert remove operation worked
    assert(repo.activate_records(group, copy(key)))
    with hdf5.HDF5File() as file:
        with pytest.raises(KeyError):
            file[dest_group]['0']
            file[dest_group]['1']
        assert(len(file[group]) == 3)
    locations = repo.find(dest_group, key)
    assert(len(key) == 3)
    assert(len(locations) == 0)

    # assert copy operation worked
    locations = repo.find(group, key)
    assert(len(key) == 0)
    assert(len(locations) == 1)
