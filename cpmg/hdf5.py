import pickle
import uuid
from collections import defaultdict

import h5py
import numpy as np

import cpmg.config as config
import cpmg.utils as utils
from cpmg.ranges import IndexKey, Key, Range, WholeRange
import cpmg.models as models


def serialize(data):
    return pickle.dumps(data).hex()


def deserialize(data):
    return pickle.loads(bytes.fromhex(data))


def serialize_chunk(data):
    return [serialize(record) for record in data]


def deserialize_chunk(data):
    return [deserialize(record) for record in data]


class HDF5File(h5py.File):
    def __init__(self, filepath=None, mode='a'):
        self.filepath = filepath or config.HDF5_FILEPATH
        super().__init__(self.filepath, mode, libver='latest')

    def __enter__(self):
        return self

    def __del__(self):
        self.close()

    def create_group(self, group, track_order=None):
        if group in self:
            return self[group]

        print(f'Created group {group}')
        return super().create_group(group, track_order=track_order)


class AbstractHDF5RepositoryImpl:
    GROUP = None
    DEFAULT_INDICES = None
    INDEX_DATASET = 'index'
    INDENT = '   '

    def __repr__(self):
        def recursive_print(name, obj):
            path = obj.name.split('/')
            path.remove('')
            num = len(path)
            print((self.INDENT * num) + name + ' - ' + str(len(obj)))

        with HDF5File() as file:
            print(f'{self.INDENT}{self.GROUP}')
            group = file[self.GROUP]
            group.visititems(recursive_print)

        return ''

    def __init__(self):
        with HDF5File() as file:
            group = file.create_group(self.GROUP)
            self._init_index(group)

    def load(self, key):
        group = self._refine_group_from_key(self.GROUP, key)
        return self._load(group, key)

    def save(self, data):
        data = utils.to_list(data)
        if len(data) == 0:
            return []

        group = self._refine_group_from_data(self.GROUP, data)
        return Key(self._save(group, data))

    def remove_records(self, key):
        return self._remove_records(self.GROUP, key)

    def remove_dataset(self, datasets):
        datasets = utils.to_list(datasets)
        with HDF5File() as file:
            group = file[self.GROUP]
            for dataset in self._get_datasets(group):
                if dataset in datasets:
                    self._remove_index_data(self.GROUP, Key(list(group[dataset].attrs.keys())))
                    del group[dataset]

        return True

    def remove_index(self, index):
        with HDF5File() as file:
            index_dataset = file[self.GROUP][self.INDEX_DATASET]
            index_name = index_dataset.attrs[index]
            del index_dataset.attrs[index]
            del file[self.GROUP][index_name]
            del self.indices[index]

    def remove_inactive_dataset(self):
        with HDF5File() as file:
            del file['inactives'][self.GROUP]

    def find(self, key):
        return self._find(self.GROUP, key)

    def deactivate_records(self, key):
        src_group = self._refine_group_from_key(self.GROUP, key)
        dest_group = self._refine_group('inactives', src_group)
        return self._move(src_group, key, dest_group)

    def activate_records(self, key):
        src_group = self._refine_group_from_key(self._refine_group('inactives', self.GROUP), key)
        return self._move(src_group, key, self.GROUP)

    def get_num_records(self, dataset=None):
        return self._get_num_records(self.GROUP, dataset)

    def load_inactivate_records(self, key):
        return self._load(self._refine_group_from_key(self._refine_group('inactives', self.GROUP), key), key)

    def get_index(self, index):
        try:
            with HDF5File() as file:
                dataset = file[self.GROUP][self.indices[index]]
                return dict(dataset.attrs.items())
        except KeyError:
            raise KeyError(f'Index {index} does not exist!')

    def get_datasets(self):
        with HDF5File() as file:
            return sorted(self._get_datasets(file[self.GROUP]))

    def create_index(self, index):
        with HDF5File() as file:
            group = file[self.GROUP]
            index_dataset = self._create_index_dataset(file[self.GROUP], group[self.INDEX_DATASET], index)
            zipped_data = zip(*self.load(Key(WholeRange())))
            ids, data = zipped_data
            for record_index, _id in zip(self._extract_index_data(index, data), ids):
                index_dataset.attrs[record_index] = _id

    def _save(self, group, data, ids=None):
        """
        To be overriden by subclasses
        """

        pass

    def _get_record(self, dataset, idx, _id):
        """
        To be overriden by subclasses
        """

        pass

    def _init_index(self, group):
        indices_dataset = group.require_dataset(self.INDEX_DATASET, shape=(0,), dtype='f', exact=True)
        self.indices = dict(indices_dataset.attrs.items())
        for index in self.DEFAULT_INDICES:
            if index not in self.indices:
                self._create_index_dataset(group, indices_dataset, index)

    def _create_index_dataset(self, group, indices_dataset, index):
        index_name = '_'.join([index, self.INDEX_DATASET])
        indices_dataset.attrs[index] = index_name
        self.indices[index] = index_name
        return group.create_dataset(index_name, shape=(0,), dtype='f')

    def _get_num_records(self, group, dataset=None):
        num_records = 0
        with HDF5File() as file:
            if dataset is None:
                for dataset in self._get_datasets(file[group]):
                    num_records += len(file[group][dataset])
            else:
                num_records += len(file[group][dataset])

        return num_records

    def _load(self, group, key):
        if isinstance(key.key, IndexKey):
            return self._load_from_ids(group, self._convert_to_id_key(key))
        else:
            return self._load_from_range(group, key)

    def _convert_to_id_key(self, key):
        if key.index == 'id':
            return key

        return Key(list(self._get_ids_from_index(key)), peptide_length=key.peptide_length)

    def _get_ids_from_index(self, key):
        try:
            with HDF5File() as file:
                dataset = file[self.GROUP][self.indices[key.index]]
                for index_val in key:
                    try:
                        yield dataset.attrs[index_val]
                    except KeyError:
                        print(f'The index {key.index} does not contain the value {index_val}! Skipping...')
        except KeyError:
            raise KeyError(f'{key.index} is not a recognized index!')

    def _load_from_range(self, group, key):
        with HDF5File() as file:
            range_index = 0
            try:
                for dataset in self._get_datasets(file[group]):
                    dataset = file[group][dataset]
                    if key in Range(range_index, range_index + len(dataset)):
                        # items must be sorted manually until track_order argument is fixed in h5py
                        for _id, idx in sorted(dataset.attrs.items(), key=lambda x: x[1]):
                            if range_index in key:
                                yield self._get_record(dataset, idx, _id)
                            range_index += 1
                    else:
                        range_index += len(dataset)
            except KeyError:
                return []

    def _load_from_ids(self, group, key):
        with HDF5File() as file:
            try:
                for dataset in self._get_datasets(file[group]):
                    dataset = file[group][dataset]
                    for _id, idx in dataset.attrs.items():
                        if _id in key:
                            yield self._get_record(dataset, idx, _id)
            except KeyError:
                return []

    def _get_datasets(self, group):
        paths = []

        def find_dataset_path(name, obj):
            if isinstance(obj, h5py.Dataset) and self.INDEX_DATASET not in name:
                paths.append(name)

        group.visititems(find_dataset_path)

        return np.sort(paths)

    def _find(self, group, key):
        id_key = self._convert_to_id_key(key)
        locations, func = self._get_finder(id_key)

        with HDF5File() as file:
            file[group].visititems(func)

        return locations

    def _get_finder(self, key):
        locations = defaultdict(list)

        def find_ids(name, obj):
            if isinstance(obj, h5py.Dataset):
                for _id, idx in obj.attrs.items():
                    if _id in key:
                        locations[obj.name].append(idx)
                if obj.name in locations:
                    locations[obj.name].sort()

        return locations, find_ids

    def _remove_records(self, group, key):
        id_key = self._convert_to_id_key(key)
        locations = self._find(group, id_key)
        with HDF5File() as file:
            for path, indices in locations.items():
                dataset = file[path]
                indices = np.delete(np.arange(len(dataset)), indices)

                if len(indices) != 0:
                    self._remove_values(dataset, indices)
                    self._remove_ids(dataset, indices)
                else:
                    del file[path]

        self._remove_index_data(group, id_key)

        return True

    def _remove_index_data(self, group, key):
        with HDF5File() as file:
            for _, index_dataset in self.indices.items():
                try:
                    index_dataset = file[group][index_dataset]
                    for index_value, _id in index_dataset.attrs.items():
                        if _id in key:
                            del index_dataset.attrs[index_value]
                except KeyError:
                    continue

    def _move(self, src_group, key, dest_group):
        ids, data = zip(*self._load(src_group, key))
        data = list(data)

        dest_group = self._refine_group_from_data(dest_group, data)
        self._save(dest_group, data, ids=ids)
        self._remove_records(src_group, key)

        return True

    def _remove_values(self, dataset, indices):
        data = dataset[indices]
        dataset[:len(data)] = data
        dataset.resize(data.shape)

    def _remove_ids(self, dataset, indices):
        items = list(filter(lambda x: x[1] in indices, dataset.attrs.items()))
        items.sort(key=lambda x: x[1])
        ids, _ = zip(*items)
        for new_idx, _id in enumerate(ids):
            dataset.attrs[_id] = new_idx

        for _id in set(dataset.attrs.keys()).difference(ids):
            del dataset.attrs[_id]

    def _get_unique_ids(self, dataset):
        while True:
            _id = str(uuid.uuid4())
            if _id not in dataset.attrs:
                yield _id

    def _get_next_dataset_name(self, group):
        max_idx = utils.get_maximum(group, int)
        max_idx = -1 if max_idx is None else max_idx
        return str(max_idx + 1)

    def _refine_group(self, group, *args):
        return '/'.join([group, *args])

    def _refine_group_from_key(self, group, key):
        return group

    def _refine_group_from_data(self, group, data):
        return group

    def _extract_index_data(self, index, data):
        for record in data:
            yield record[index]


class HDF5ObjRepositoryImpl(AbstractHDF5RepositoryImpl):
    ENCODING = 'ascii'

    def _save(self, group, data, ids=None):
        serialized_data = serialize_chunk(data)
        max_bin_len = utils.get_maximum(serialized_data, len)

        with HDF5File() as file:
            dataset = self._create_dataset(file.create_group(group), len(serialized_data), str(max_bin_len))
            ids = self._write(dataset, serialized_data, ids or self._get_unique_ids(dataset))

        self._update_indices(data, ids)

        return ids

    def _create_dataset(self, group, rows, dlen):
        dataset = self._get_next_dataset_name(group)
        return group.create_dataset(dataset, (rows,), maxshape=(None,), chunks=True, dtype='S' + dlen,
                                    compression=config.COMPRESSION, compression_opts=config.COMPRESSION_OPTS)

    def _update_indices(self, data, ids):
        with HDF5File() as file:
            for index, index_dataset in self.indices.items():
                index_dataset = file[self.GROUP][index_dataset]
                for record_index, _id in zip(self._extract_index_data(index, data), ids):
                    index_dataset.attrs[record_index] = _id

    def _write(self, dataset, data, ids):
        used_ids = []
        for i, (record, _id) in enumerate(zip(data, ids)):
            dataset[i] = record.encode(self.ENCODING)
            dataset.attrs[_id] = i
            used_ids.append(_id)

        return used_ids

    def _get_record(self, dataset, idx, _id):
        return (_id, deserialize(dataset[idx].decode(self.ENCODING)))


class HDF5ArrayRepositoryImpl(AbstractHDF5RepositoryImpl):

    def _save(self, group, data, ids=None):
        with HDF5File() as file:
            dataset = self._create_dataset(file.create_group(group), data)
            return self._write(dataset, ids or self._get_unique_ids(dataset))

    def _create_dataset(self, group, data):
        dataset = self._get_next_dataset_name(group)
        return group.create_dataset(dataset, data=data, chunks=True, compression=config.COMPRESSION,
                                    compression_opts=config.COMPRESSION_OPTS)

    def _write(self, dataset, ids):
        used_ids = []
        for i, (_, _id) in enumerate(zip(dataset, ids)):
            dataset.attrs[_id] = i
            used_ids.append(_id)

        return used_ids

    def _get_record(self, dataset, idx, _id):
        return (_id, dataset[idx])


class ConnectionHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.Connection.STRING
    DEFAULT_INDICES = ['kekule']


class BackboneHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.Backbone.STRING
    DEFAULT_INDICES = ['mapped_kekule']


class TemplateHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.Template.STRING
    DEFAULT_INDICES = ['kekule']


class SidechainHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.Sidechain.STRING
    DEFAULT_INDICES = ['kekule']


class MonomerHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.Monomer.STRING
    DEFAULT_INDICES = ['kekule']


class PeptideHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.Peptide.STRING
    DEFAULT_INDICES = ['kekule']

    def _refine_group_from_key(self, group, key):
        try:
            peptide_length = key.peptide_length
            if peptide_length is not None:
                group = self._refine_group(group, str(peptide_length))
        except AttributeError:
            pass

        return group

    def _refine_group_from_data(self, group, data):
        peptide_length = str(data[0]['length'])
        return self._refine_group(group, peptide_length)


class TemplatePeptideHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.TemplatePeptide.STRING
    DEFAULT_INDICES = ['kekule']

    def _refine_group_from_key(self, group, key):
        try:
            peptide_length = key.peptide_length
            if peptide_length is not None:
                group = self._refine_group(group, str(peptide_length))
        except AttributeError:
            pass

        return group

    def _refine_group_from_data(self, group, data):
        peptide_length = str(data[0]['peptide']['length'])
        return self._refine_group(group, peptide_length)


class ReactionHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.Reaction.STRING
    DEFAULT_INDICES = []


class RegioSQMHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.RegioSQMPrediction.STRING
    DEFAULT_INDICES = ['reacting_mol']


class pKaHDF5Repository(HDF5ObjRepositoryImpl):
    GROUP = models.pKaPrediction.STRING
    DEFAULT_INDICES = ['reacting_mol']


class PeptidePlanHDF5Repository(HDF5ArrayRepositoryImpl):
    GROUP = models.PeptidePlan.STRING
    DEFAULT_INDICES = []

    def _save(self, group, data, ids=None):
        try:
            temp_group = self._refine_group(group, 'no_cap')
            _ids = super()._save(temp_group, data.reg_combos, ids)

            temp_group = self._refine_group(group, 'cap')
            _ids.extend(super()._save(temp_group, data.cap_combos, ids))
        except AttributeError:  # used when moving data
            peptide_length = [int(s) for s in group.split('/') if s.isdigit()][0]
            no_cap_data, cap_data = utils.split(data, lambda x: len(x) == peptide_length)

            if len(no_cap_data) != 0:
                temp_group = self._refine_group(group, 'no_cap')
                _ids = super()._save(temp_group, np.array(no_cap_data), ids[:len(no_cap_data)])

            if len(cap_data) != 0:
                temp_group = self._refine_group(group, 'cap')
                _ids.extend(super()._save(temp_group, np.array(cap_data), ids[len(no_cap_data):]))

        return _ids

    def _refine_group_from_key(self, group, key):
        try:
            peptide_length = key.peptide_length
            if peptide_length is not None:
                group = self._refine_group(group, str(peptide_length))
        except AttributeError:
            pass

        return group

    def _refine_group_from_data(self, group, data):
        try:
            peptide_length = str(data.length)
        except AttributeError:
            lengths = np.unique(list(map(len, data)))
            if len(lengths) == 1:
                peptide_length = str(lengths[0])
            else:
                peptide_length = str(min(lengths))

        return self._refine_group(group, peptide_length)


class HDF5Repository:
    def __init__(self):
        self.connection_repo = ConnectionHDF5Repository()
        self.backbone_repo = BackboneHDF5Repository()
        self.template_repo = TemplateHDF5Repository()
        self.sidechain_repo = SidechainHDF5Repository()
        self.monomer_repo = MonomerHDF5Repository()
        self.peptide_repo = PeptideHDF5Repository()
        self.template_peptide_repo = TemplatePeptideHDF5Repository()
        self.reaction_repo = ReactionHDF5Repository()
        self.regiosqm_repo = RegioSQMHDF5Repository()
        self.pka_repo = pKaHDF5Repository()
        self.peptide_plan_repo = PeptidePlanHDF5Repository()

    def __repr__(self):

        print('/')
        for _, instance in self.__dict__.items():
            instance.__repr__()

        return ''

    @classmethod
    def instance(cls):
        if not hasattr(cls, '_instance'):
            cls._instance = cls()
        return cls._instance
