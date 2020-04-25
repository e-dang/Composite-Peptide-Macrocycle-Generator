import pickle
import h5py
import numpy as np
import uuid
import macrocycles.config as config


def get_maximum(data, func):
    try:
        return np.max(list(map(func, data)))
    except ValueError:
        return None


def serialize(data):
    return pickle.dumps(data).hex()


def deserialize(data):
    return pickle.loads(bytes.fromhex(data))


def serialize_chunk(data):
    return [serialize(doc) for doc in data]


def deserialize_chunk(data):
    return [deserialize(doc) for doc in data]


def to_list(data):
    if isinstance(data, dict):
        return [data]

    return data


class HDF5File(h5py.File):
    def __init__(self, filepath=None, mode='a'):
        self.filepath = config.HDF5_FILEPATH if filepath is None else filepath
        super().__init__(self.filepath, mode, libver='latest')

    def __enter__(self):
        return self

    def __del__(self):
        self.close()


class HDF5Initializer:
    def __init__(self):
        self.data_types = ['connections', 'templates', 'sidechains', 'monomers', 'peptides']

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
        serialized_data = serialize_chunk(to_list(data))
        max_bin_len = get_maximum(serialized_data, len)

        _ids = []
        with HDF5File() as file:
            dataset = self._create_dataset(file[group], len(serialized_data), str(max_bin_len))
            for i, s_doc in enumerate(serialized_data):
                dataset[i] = s_doc.encode(self.ENCODING)
                _ids.append(str(uuid.uuid4()))
                dataset.attrs[_ids[-1]] = i

        return _ids

    def _create_dataset(self, group, rows, dlen):
        max_idx = get_maximum(group, int)
        max_idx = -1 if max_idx is None else max_idx
        return group.create_dataset(str(max_idx + 1), (rows,), dtype='S' + dlen, compression=config.COMPRESSION, compression_opts=config.COMPRESSION_OPTS)

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
