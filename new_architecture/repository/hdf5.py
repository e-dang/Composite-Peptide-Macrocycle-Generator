import pickle
import uuid
from collections import defaultdict

import h5py
import numpy as np

import macrocycles.config as config
import new_architecture.utils as utils


def serialize(data):
    return pickle.dumps(data).hex()


def deserialize(data):
    return pickle.loads(bytes.fromhex(data))


def serialize_chunk(data):
    return [serialize(doc) for doc in data]


def deserialize_chunk(data):
    return [deserialize(doc) for doc in data]


class HDF5File(h5py.File):
    def __init__(self, filepath=None, mode='a'):
        self.filepath = config.HDF5_FILEPATH if filepath is None else filepath
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


class HDF5Initializer:
    def __init__(self):
        self.data_types = ['backbone', 'connections', 'templates',
                           'sidechains', 'monomers', 'peptides', 'template_peptides']

    def initialize(self):
        with HDF5File() as file:
            for data_type in self.data_types:
                file.create_group(data_type)


class HDF5Repository:
    ENCODING = 'ascii'

    def load(self, group, key):
        if isinstance(key, list) and isinstance(key[0], str):
            return self._load_ids(group, key)
        else:
            return self._load_range(group, key)

    def save(self, group, data):
        data = utils.to_list(data)
        if len(data) == 0:
            return []

        serialized_data = serialize_chunk(data)
        max_bin_len = utils.get_maximum(serialized_data, len)

        with HDF5File() as file:
            dataset = self._create_dataset(file[group], len(serialized_data), str(max_bin_len))
            ids = self._write(dataset, serialized_data, self._get_unique_ids(dataset))

        return ids

    def find(self, group, key):
        locations = defaultdict(list)
        with HDF5File() as file:
            for dataset in file[group]:
                dataset = file[group][dataset]
                for _id, idx in dataset.attrs.items():
                    if _id in key:
                        locations[dataset.name].append(idx)
                        key.remove(_id)
                if dataset.name in locations:
                    locations[dataset.name].sort()

        return locations

    def remove(self, group, key):
        locations = self.find(group, key)
        with HDF5File() as file:
            for path, indices in locations.items():
                dataset = file[path]
                indices = np.delete(np.arange(len(dataset)), indices)

                if len(indices) != 0:
                    self._remove_values(dataset, indices)
                    self._remove_ids(dataset, indices)
                    # # remove values
                    # data = dataset[indices]
                    # dataset[:len(data)] = data
                    # dataset.resize((len(data),))

                    # # remove ids
                    # items = list(filter(lambda x: x[1] in indices, dataset.attrs.items()))
                    # items.sort(key=lambda x: x[1])
                    # ids, _ = zip(*items)
                    # for new_idx, _id in enumerate(ids):
                    #     dataset.attrs[_id] = new_idx

                    # for _id in set(dataset.attrs.keys()).difference(ids):
                    #     del dataset.attrs[_id]
                else:
                    del file[path]

        return True

    def move(self, src_group, key, dest_group):
        ids, data = zip(*self.load(src_group, key))
        data = serialize_chunk(data)
        max_bin_len = utils.get_maximum(data, len)

        with HDF5File() as file:
            dataset = self._create_dataset(file.create_group(dest_group), len(ids), str(max_bin_len))
            self._write(dataset, data, ids)
            self.remove(src_group, key)

        return True

    def deactivate_records(self, group, key):
        return self.move(group, key, '/'.join(['inactives', group]))

    def activate_records(self, group, key):
        return self.move('/'.join(['inactives', group]), key, group)

    def get_num_records(self, group):
        num_records = 0
        with HDF5File() as file:
            for dataset in file[group]:
                num_records += len(file[group][dataset])

        return num_records

    def _write(self, dataset, data, ids):
        used_ids = []
        for i, (doc, _id) in enumerate(zip(data, ids)):
            dataset[i] = doc.encode(self.ENCODING)
            dataset.attrs[_id] = i
            used_ids.append(_id)

        return used_ids

    def _remove_values(self, dataset, indices):
        data = dataset[indices]
        dataset[:len(data)] = data
        dataset.resize((len(data),))

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

    def _create_dataset(self, group, rows, dlen):
        max_idx = utils.get_maximum(group, int)
        max_idx = -1 if max_idx is None else max_idx
        return group.create_dataset(str(max_idx + 1), (rows,), maxshape=(None,), chunks=True, dtype='S' + dlen, compression=config.COMPRESSION, compression_opts=config.COMPRESSION_OPTS)

    def _load_range(self, group, key):
        with HDF5File() as file:
            range_index = 1
            for dataset in file[group]:
                dataset = file[group][dataset]
                for _id, idx in dataset.attrs.items():
                    if range_index in key:
                        yield (_id, deserialize(dataset[idx].decode(self.ENCODING)))
                    range_index += 1

    def _load_ids(self, group, key):
        with HDF5File() as file:
            for dataset in file[group]:
                dataset = file[group][dataset]
                for _id, idx in dataset.attrs.items():
                    if _id in key:
                        yield (_id, deserialize(dataset[idx].decode(self.ENCODING)))
