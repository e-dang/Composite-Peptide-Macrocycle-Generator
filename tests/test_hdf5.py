import os

import h5py
import pytest
from rdkit import Chem

# import conftest
import cpmg.config as config
import cpmg.hdf5 as hdf5
import cpmg.ranges as ranges


@pytest.fixture()
def dict1():
    return {'A': 15, 'B': False, 'C': -1.0, 'D': b'test_bin_string',
            'E': 'test_string', 'F': Chem.MolFromSmiles('C[NH]C').ToBinary()}


@pytest.fixture()
def dict2():
    return {'A': 1, 'B': True, 'C': 11.0, 'D': b'test_bin_string',
            'E': 'another_test_string', 'F': Chem.MolFromSmiles('C([OH])C').ToBinary()}


@pytest.fixture()
def dict3():
    return {'A': -1, 'B': True, 'C': 1.0, 'D': b'test_bin_string',
            'E': 'test_string2', 'F': Chem.MolFromSmiles('CCC').ToBinary()}


@pytest.fixture()
def list1(dict1, dict2, dict3):
    key = lambda x: x['A']
    return sorted([dict1, dict2, dict3], key=key), key


@pytest.fixture(params=[
    ('connection_dicts', hdf5.ConnectionHDF5Repository, ['0', 'index', 'completed', 'kekule_index']),
    ('backbone_dicts', hdf5.BackboneHDF5Repository, ['0', 'index', 'completed', 'mapped_kekule_index']),
    ('template_dicts', hdf5.TemplateHDF5Repository, ['0', 'index', 'completed', 'kekule_index']),
    ('sidechain_dicts', hdf5.SidechainHDF5Repository, ['0', 'index', 'completed', 'kekule_index']),
    ('monomer_dicts', hdf5.MonomerHDF5Repository, ['0', 'index', 'completed', 'kekule_index']),
    ('peptide_len_3_dicts', hdf5.PeptideHDF5Repository, ['3/0', 'index', 'completed', 'kekule_index']),
    ('peptide_len_4_dicts', hdf5.PeptideHDF5Repository, ['4/0', 'index', 'completed', 'kekule_index']),
    ('peptide_len_5_dicts', hdf5.PeptideHDF5Repository, ['5/0', 'index', 'completed', 'kekule_index']),
    ('template_peptide_len_3_dicts', hdf5.TemplatePeptideHDF5Repository, ['3/0', 'index', 'completed', 'kekule_index']),
    ('template_peptide_len_4_dicts', hdf5.TemplatePeptideHDF5Repository, ['4/0', 'index', 'completed', 'kekule_index']),
    ('template_peptide_len_5_dicts', hdf5.TemplatePeptideHDF5Repository, ['5/0', 'index', 'completed', 'kekule_index']),
    ('reaction_dicts', hdf5.ReactionHDF5Repository, ['0', 'index', 'completed']),
    ('regiosqm_dicts', hdf5.RegioSQMHDF5Repository, ['0', 'index', 'completed', 'reacting_mol_index']),
    ('pka_dicts', hdf5.pKaHDF5Repository, ['0', 'index', 'completed', 'reacting_mol_index'])
],
    ids=['connections', 'backbones', 'templates', 'sidechains', 'monomers', 'peptides_3', 'peptides_4', 'peptides_5',
         'template_peptides_3', 'template_peptides_4', 'template_peptides_5', 'reactions', 'regiosqm', 'pka'
         ])
def obj_data(request):
    fixture, repo, datasets = request.param
    fixture = request.getfixturevalue(fixture)
    return fixture, repo(), sorted(datasets)


@pytest.fixture(params=[
    ('peptide_plan_tuple_len_3', hdf5.PeptidePlanHDF5Repository, ['3/cap/0', '3/no_cap/0', 'completed', 'index']),
    ('peptide_plan_tuple_len_4', hdf5.PeptidePlanHDF5Repository, ['4/cap/0', '4/no_cap/0', 'completed', 'index']),
    ('peptide_plan_tuple_len_5', hdf5.PeptidePlanHDF5Repository, ['5/cap/0', '5/no_cap/0', 'completed', 'index'])
],
    ids=['peptide_plan_tuple_3', 'peptide_plan_tuple_4', 'peptide_plan_tuple_5'])
def array_data(request):
    fixture, repo, datasets = request.param
    fixture = request.getfixturevalue(fixture)
    return fixture, repo(), datasets


@pytest.fixture()
def obj_repository_w_saved_data(obj_data):
    dicts, repo, datasets = obj_data

    ids = repo.save(dicts)

    return dicts, repo, datasets, ids


@pytest.fixture()
def array_repository_w_saved_data(array_data):
    arrays, repo, datasets = array_data

    ids = repo.save(arrays)

    return arrays, repo, datasets, ids


@pytest.fixture()
def obj_repository_w_inactive_records(obj_repository_w_saved_data):
    dicts, repo, datasets, ids = obj_repository_w_saved_data

    repo.deactivate_records(ids)

    return dicts, repo, datasets, ids


@pytest.fixture()
def array_repository_w_inactive_records(array_repository_w_saved_data):
    arrays, repo, datasets, ids = array_repository_w_saved_data

    repo.deactivate_records(ids)

    return arrays, repo, datasets, ids


def recursive_search(group, group_name):
    found_datasets = []
    for obj in group:
        obj = group[obj]
        if isinstance(obj, h5py.Dataset):
            name = obj.name
            name = name.split('/')
            name.remove(group_name)
            name = '/'.join(name)
            found_datasets.append(name[1:])
        elif isinstance(obj, h5py.Group):
            found_datasets.extend(recursive_search(obj, group_name))

    return found_datasets


def remove_meta_data_datasets(repo, datasets):
    for dataset in list(datasets):
        if repo.INDEX_DATASET in dataset:
            datasets.remove(dataset)
        elif repo.COMPLETE_DATASET in dataset:
            datasets.remove(dataset)


def test_serialize_deserialize(list1):
    ls, _ = list1
    for doc in ls:
        serialized_data = hdf5.serialize(doc)

        assert isinstance(serialized_data, str)

        deserialized_data = hdf5.deserialize(serialized_data)

        assert doc == deserialized_data


def test_serialize_deserialize_chunk(list1):
    ls, key = list1
    serialized_data = hdf5.serialize_chunk(ls)

    for serialized_str in serialized_data:
        assert isinstance(serialized_str, str)

    deserialized_data = hdf5.deserialize_chunk(serialized_data)
    deserialized_data.sort(key=key)

    assert deserialized_data == ls


def test_hdf5_file_regular_open():
    file = hdf5.HDF5File()

    assert os.path.exists(config.HDF5_FILEPATH)
    assert file
    assert len(list(file.keys())) == 0

    file.close()

    assert not file


def test_hdf5_file_context_manager():
    with hdf5.HDF5File() as file:
        assert os.path.exists(config.HDF5_FILEPATH)
        assert file
        assert len(list(file.keys())) == 0

    assert not file


def test_hdf5_file_create_group():
    group_name = 'test'
    with hdf5.HDF5File() as file:
        group1 = file.create_group(group_name)
        group2 = file.create_group(group_name)

        assert isinstance(group1, h5py.Group)
        assert isinstance(group2, h5py.Group)
        assert group1.name == group2.name


def test_hdf5_obj_repository_save(obj_repository_w_saved_data):
    dicts, repo, datasets, _ = obj_repository_w_saved_data

    # verify datasets have been made properly
    with hdf5.HDF5File() as file:
        found_datasets = sorted(recursive_search(file[repo.GROUP], repo.GROUP))
        assert found_datasets == datasets
        assert len(file[repo.GROUP][datasets[0]]) == len(dicts)

    # verify indices have been made properly
    for index in repo.indices:
        stored_index = repo.get_index(index)
        for doc in dicts:
            attr_val = doc[index]
            assert attr_val in stored_index
            stored_index.pop(attr_val)


def test_hdf5_obj_repository_load_from_whole_range(obj_repository_w_saved_data, dict_sort_key):
    dicts, repo, _, ids = obj_repository_w_saved_data

    key = ranges.Key(ranges.WholeRange())
    loaded_ids, loaded_data = zip(*list(repo.load(key)))
    loaded_ids = list(loaded_ids)
    loaded_data = list(loaded_data)

    loaded_ids.sort()
    loaded_data.sort(key=dict_sort_key)
    dicts.sort(key=dict_sort_key)

    assert loaded_data == dicts
    assert ids == loaded_ids


def test_hdf5_obj_repository_load_from_ids(obj_repository_w_saved_data, dict_sort_key):
    dicts, repo, _, ids = obj_repository_w_saved_data

    loaded_ids, loaded_data = zip(*list(repo.load(ids)))
    loaded_ids = list(loaded_ids)
    loaded_data = list(loaded_data)

    loaded_ids.sort()
    loaded_data.sort(key=dict_sort_key)
    dicts.sort(key=dict_sort_key)

    assert loaded_data == dicts
    assert ids == loaded_ids


def test_hdf5_obj_repository_load_from_index_vals(obj_repository_w_saved_data, extract_attr_from_dicts, dict_sort_key):
    dicts, repo, _, ids = obj_repository_w_saved_data

    for index in repo.DEFAULT_INDICES:
        key = ranges.Key(extract_attr_from_dicts(dicts, index), index=index)
        loaded_ids, loaded_data = zip(*list(repo.load(key)))
        loaded_ids = list(loaded_ids)
        loaded_data = list(loaded_data)

        loaded_ids.sort()
        loaded_data.sort(key=dict_sort_key)
        dicts.sort(key=dict_sort_key)

        assert loaded_data == dicts
        assert ids == loaded_ids


def test_hd5f_array_repository_save(array_repository_w_saved_data):
    arrays, repo, datasets, _ = array_repository_w_saved_data

    num_records = len(arrays.reg_combos) + len(arrays.cap_combos)

    # verify datasets have been made properly
    with hdf5.HDF5File() as file:
        found_datasets = recursive_search(file[repo.GROUP], repo.GROUP)
        assert found_datasets == datasets
        for dataset in datasets:
            if repo.INDEX_DATASET not in dataset:
                num_records -= len(file[repo.GROUP][dataset])
        assert num_records == 0


def test_hdf5_array_repository_load_from_whole_range(array_repository_w_saved_data):
    arrays, repo, _, ids = array_repository_w_saved_data

    arrays = [list(record) for record in arrays.reg_combos + arrays.cap_combos]

    key = ranges.Key(ranges.WholeRange())
    loaded_ids, loaded_data = zip(*list(repo.load(key)))
    loaded_ids = list(loaded_ids)
    loaded_data = [list(record) for record in loaded_data]

    loaded_ids.sort()
    loaded_data.sort()
    arrays.sort()

    assert loaded_data == arrays
    assert ids == loaded_ids


def test_hdf5_array_repository_load_from_ids(array_repository_w_saved_data):
    arrays, repo, _, ids = array_repository_w_saved_data

    arrays = [list(record) for record in arrays.reg_combos + arrays.cap_combos]

    loaded_ids, loaded_data = zip(*list(repo.load(ids)))
    loaded_ids = list(loaded_ids)
    loaded_data = [list(record) for record in loaded_data]

    loaded_ids.sort()
    loaded_data.sort()
    arrays.sort()

    assert loaded_data == arrays
    assert ids == loaded_ids


def test_hdf5_obj_repository_get_num_records(obj_repository_w_saved_data):
    _, repo, _, ids = obj_repository_w_saved_data

    assert repo.get_num_records() == len(ids)


def test_hdf5_array_repository_get_num_records(array_repository_w_saved_data):
    _, repo, _, ids = array_repository_w_saved_data

    assert repo.get_num_records() == len(ids)


def test_hdf5_obj_repository_get_datasets(obj_repository_w_saved_data):
    _, repo, datasets, ids = obj_repository_w_saved_data

    with hdf5.HDF5File() as file:
        found_datasets = recursive_search(file[repo.GROUP], repo.GROUP)
        remove_meta_data_datasets(repo, found_datasets)

    assert repo.get_datasets() == found_datasets


def test_hdf5_array_repository_get_datasets(array_repository_w_saved_data):
    _, repo, datasets, ids = array_repository_w_saved_data

    with hdf5.HDF5File() as file:
        found_datasets = recursive_search(file[repo.GROUP], repo.GROUP)
        remove_meta_data_datasets(repo, found_datasets)

    assert repo.get_datasets() == found_datasets


def test_hdf5_obj_repository_find(obj_repository_w_saved_data):
    _, repo, _, ids = obj_repository_w_saved_data

    locations = repo.find(ids)

    found_ids = []
    for path, indices in locations.items():
        with hdf5.HDF5File() as file:
            dataset = file[path]
            for _id, idx in dataset.attrs.items():
                if idx in indices:
                    found_ids.append(_id)

    assert sorted(found_ids) == ids


def test_hdf5_array_repository_find(array_repository_w_saved_data):
    _, repo, _, ids = array_repository_w_saved_data

    locations = repo.find(ids)

    found_ids = []
    for path, indices in locations.items():
        with hdf5.HDF5File() as file:
            dataset = file[path]
            for _id, idx in dataset.attrs.items():
                if idx in indices:
                    found_ids.append(_id)

    assert sorted(found_ids) == ids


def test_hdf5_obj_repository_remove_records(obj_repository_w_saved_data):
    _, repo, _, ids = obj_repository_w_saved_data

    prev_indices = {}
    for index in repo.indices:
        prev_indices[index] = repo.get_index(index)

    removed_ids = ranges.Key([ids[0]])
    repo.remove_records(removed_ids)

    removed_locations = repo.find(removed_ids)
    locations = repo.find(ids)

    found_ids = []
    for path, indices in locations.items():
        with hdf5.HDF5File() as file:
            dataset = file[path]
            for _id, idx in dataset.attrs.items():
                if idx in indices:
                    found_ids.append(_id)

    assert repo.get_num_records() == len(ids) - 1
    assert len(removed_locations) == 0
    assert len(found_ids) == len(ids) - 1
    assert sorted(found_ids) != ids

    for index in repo.indices:
        assert prev_indices[index] != repo.get_index(index)


def test_hdf5_array_repository_remove_records(array_repository_w_saved_data):
    _, repo, _, ids = array_repository_w_saved_data

    prev_indices = {}
    for index in repo.indices:
        prev_indices[index] = repo.get_index(index)

    removed_ids = ranges.Key([ids[0]])
    repo.remove_records(removed_ids)

    removed_locations = repo.find(removed_ids)
    locations = repo.find(ids)

    found_ids = []
    for path, indices in locations.items():
        with hdf5.HDF5File() as file:
            dataset = file[path]
            for _id, idx in dataset.attrs.items():
                if idx in indices:
                    found_ids.append(_id)

    assert repo.get_num_records() == len(ids) - 1
    assert len(removed_locations) == 0
    assert len(found_ids) == len(ids) - 1
    assert sorted(found_ids) != ids

    for index in repo.indices:
        assert prev_indices[index] != repo.get_index(index)


def test_hdf5_obj_repository_remove_dataset(obj_repository_w_saved_data):
    _, repo, datasets, ids = obj_repository_w_saved_data

    prev_indices = {}
    for index in repo.indices:
        prev_indices[index] = repo.get_index(index)

    repo.remove_dataset(datasets[0])

    with hdf5.HDF5File() as file:
        group = file[repo.GROUP]
        found_datasets = sorted(list(group.keys()))

        assert datasets != found_datasets
        for index in repo.indices:
            assert prev_indices[index] != repo.get_index(index)


def test_hdf5_array_repository_remove_dataset(array_repository_w_saved_data):
    _, repo, datasets, ids = array_repository_w_saved_data

    prev_indices = {}
    for index in repo.indices:
        prev_indices[index] = repo.get_index(index)

    repo.remove_dataset(datasets[0])

    with hdf5.HDF5File() as file:
        group = file[repo.GROUP]
        found_datasets = sorted(list(group.keys()))

        assert datasets != found_datasets
        for index in repo.indices:
            assert prev_indices[index] != repo.get_index(index)


def test_hdf5_obj_repository_remove_index(obj_repository_w_saved_data):
    _, repo, datasets, ids = obj_repository_w_saved_data

    prev_indices = list(repo.indices)

    if len(prev_indices) > 0:
        repo.remove_index(prev_indices[0])

        assert list(repo.indices) != prev_indices


def test_hdf5_array_repository_remove_index(array_repository_w_saved_data):
    _, repo, datasets, ids = array_repository_w_saved_data

    prev_indices = list(repo.indices)

    if len(prev_indices) > 0:
        repo.remove_index(prev_indices[0])

        assert list(repo.indices) != prev_indices


def test_hdf5_obj_repository_deactivate_records(obj_repository_w_inactive_records):
    _, repo, datasets, ids = obj_repository_w_inactive_records

    found_datasets = repo.get_datasets()
    num_records = repo.get_num_records()
    remove_meta_data_datasets(repo, datasets)

    assert found_datasets != datasets
    assert num_records != len(ids)
    for index in repo.indices:
        assert len(repo.get_index(index)) == 0


def test_hdf5_array_repository_deactivate_records(array_repository_w_inactive_records):
    _, repo, datasets, ids = array_repository_w_inactive_records

    found_datasets = repo.get_datasets()
    num_records = repo.get_num_records()
    remove_meta_data_datasets(repo, datasets)

    assert found_datasets != datasets
    assert num_records != len(ids)
    for index in repo.indices:
        assert len(repo.get_index(index)) == 0


def test_hdf5_obj_repository_activate_records(obj_repository_w_inactive_records):
    _, repo, datasets, ids = obj_repository_w_inactive_records

    repo.activate_records(ids)

    found_datasets = repo.get_datasets()
    num_records = repo.get_num_records()
    remove_meta_data_datasets(repo, datasets)

    assert found_datasets == datasets
    assert num_records == len(ids)
    for index in repo.indices:
        assert len(repo.get_index(index)) == len(ids)


def test_hdf5_array_repository_activate_records(array_repository_w_inactive_records):
    _, repo, datasets, ids = array_repository_w_inactive_records

    repo.activate_records(ids)

    found_datasets = repo.get_datasets()
    num_records = repo.get_num_records()
    remove_meta_data_datasets(repo, datasets)

    assert found_datasets == datasets
    assert num_records == len(ids)
    for index in repo.indices:
        assert len(repo.get_index(index)) == len(ids)
