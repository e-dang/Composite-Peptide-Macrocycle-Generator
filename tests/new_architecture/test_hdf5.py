import pytest
import os
import new_architecture.repository.hdf5 as hdf5
from new_architecture.repository.repository import WholeRange
import macrocycles.config as config
from tests.new_architecture.data.mols import *
from rdkit import Chem

GROUP = 'sidechains'
TEST_DICT = {'A': 1, 'B': True, 'C': 1.0, 'D': b'test_bin_string',
             'E': 'test_string', 'F': Chem.MolFromSmiles('CCC').ToBinary()}
TEST_LIST = [TEST_DICT, TEST_DICT, TEST_DICT]


@pytest.fixture()
def filepath(monkeypatch):
    fp = os.path.join(config.PROJECT_DIR, 'tests', 'new_architecture', 'data', 'test_file.hdf5')
    monkeypatch.setattr(hdf5.config, "HDF5_FILEPATH", fp)
    yield fp
    os.remove(fp)


@pytest.fixture()
def initialize_repo(filepath):
    initializer = hdf5.HDF5Initializer()
    initializer.initialize()
    yield (initializer, filepath)


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


def test_hdf5_file_regular_open(filepath):
    file = hdf5.HDF5File()

    assert(os.path.exists(filepath))
    assert(file)
    assert(len(list(file.keys())) == 0)
    file.create_group('test_group')
    assert(len(list(file.keys())) == 1)

    file.close()

    assert(not file)


def test_hdf5_context_manager(filepath):
    with hdf5.HDF5File() as file:
        assert(os.path.exists(filepath))
        assert(file)
        assert(len(list(file.keys())) == 0)
        file.create_group('test_group')
        assert(len(list(file.keys())) == 1)

    assert(not file)


def test_hdf5_initializer(initialize_repo):
    initializer, _ = initialize_repo
    with hdf5.HDF5File() as file:
        assert(sorted(list(file.keys())) == sorted(initializer.data_types))


def test_hdf5_repository_save_load_range_single(initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    repo.save('sidechains', TEST_SIDECHAIN_1)

    with hdf5.HDF5File(filepath) as file:
        assert(list(file[GROUP].keys()) == ['0'])
        assert(file[GROUP]['0'].size == (1,))

    _, data = zip(*list(repo.load(GROUP, WholeRange())))
    assert(len(data) == 1)
    assert(data[0] == TEST_SIDECHAIN_1)


def test_hdf5_repository_save_load_ids_single(initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    _ids = repo.save('sidechains', TEST_SIDECHAIN_1)

    with hdf5.HDF5File(filepath) as file:
        assert(list(file[GROUP].keys()) == ['0'])
        assert(file[GROUP]['0'].size == (1,))

    _, data = zip(*list(repo.load(GROUP, _ids)))
    assert(len(data) == 1)
    assert(data[0] == TEST_SIDECHAIN_1)


def test_hdf5_repository_save_load_range_multi(initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    repo.save('sidechains', [TEST_SIDECHAIN_1, TEST_SIDECHAIN_2])

    with hdf5.HDF5File(filepath) as file:
        assert(list(file[GROUP].keys()) == ['0'])
        assert(file[GROUP]['0'].size == (2,))

    _, data = zip(*list(repo.load(GROUP, WholeRange())))
    data = sorted(list(data), key=lambda x: x['kekule'])
    assert(len(data) == 2)
    assert(data[0] == TEST_SIDECHAIN_1)
    assert(data[1] == TEST_SIDECHAIN_2)


def test_hdf5_repository_save_load_ids_multi(initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    _ids = repo.save('sidechains', [TEST_SIDECHAIN_1, TEST_SIDECHAIN_2])

    with hdf5.HDF5File(filepath) as file:
        assert(list(file[GROUP].keys()) == ['0'])
        assert(file[GROUP]['0'].size == (2,))

    _, data = zip(*list(repo.load(GROUP, _ids)))
    data = sorted(list(data), key=lambda x: x['kekule'])
    assert(len(data) == 2)
    assert(data[0] == TEST_SIDECHAIN_1)
    assert(data[1] == TEST_SIDECHAIN_2)


def test_hdf5_repository_save_load_range_multi_separate(initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    repo.save('sidechains', TEST_SIDECHAIN_1)
    repo.save('sidechains', TEST_SIDECHAIN_2)

    with hdf5.HDF5File(filepath) as file:
        assert(list(file[GROUP].keys()) == ['0', '1'])
        assert(file[GROUP]['0'].size == (1,))
        assert(file[GROUP]['1'].size == (1,))

    _, data = zip(*list(repo.load(GROUP, WholeRange())))
    assert(len(data) == 2)
    assert(data[0] == TEST_SIDECHAIN_1)
    assert(data[1] == TEST_SIDECHAIN_2)


def test_hdf5_repository_save_load_ids_multi_separate(initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    _ids = repo.save('sidechains', TEST_SIDECHAIN_1)
    _ids.extend(repo.save('sidechains', TEST_SIDECHAIN_2))

    with hdf5.HDF5File(filepath) as file:
        assert(list(file[GROUP].keys()) == ['0', '1'])
        assert(file[GROUP]['0'].size == (1,))
        assert(file[GROUP]['1'].size == (1,))

    _, data = zip(*list(repo.load(GROUP, _ids)))
    assert(len(data) == 2)
    assert(data[0] == TEST_SIDECHAIN_1)
    assert(data[1] == TEST_SIDECHAIN_2)


@pytest.mark.parametrize('dataset1,dataset2,expected_result,initialize_repo', [([TEST_MONOMER_1, TEST_MONOMER_2, TEST_MONOMER_3], [TEST_MONOMER_4, TEST_MONOMER_5, TEST_MONOMER_6], 6, ''), ([], [], 0, '')], indirect=['initialize_repo'])
def test_hdf5_get_num_records(dataset1, dataset2, expected_result, initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    _ids = repo.save('monomers', dataset1)
    _ids.extend(repo.save('monomers', dataset2))

    assert(repo.get_num_records('monomers') == expected_result)


def test_hdf5_find(initialize_repo):
    _, filepath = initialize_repo
    repo = hdf5.HDF5Repository()
    _ids = repo.save('monomers', [TEST_MONOMER_1, TEST_MONOMER_2, TEST_MONOMER_3])
    _ids.extend(repo.save('monomers', [TEST_MONOMER_4, TEST_MONOMER_5, TEST_MONOMER_6]))

    locations = repo.find('monomers', [_ids[0], _ids[3], _ids[4]])
    values = sorted(list(locations.values()))
    for value in values:
        value.sort()

    assert(len(locations) == 2)
    assert(sorted(list(locations.keys())) == ['0', '1'])
    assert(values == [[0], [0, 1]])
