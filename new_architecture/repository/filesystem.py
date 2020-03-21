import macrocycles.config as config
import macrocycles.utils as utils
from bson import json_util
import json
import glob
import os


class GenericDAO:
    def __init__(self, io_context):
        self.io_context = io_context

    def save(self, save_request):
        self.io_context.save(save_request)

    def load(self, load_request):
        return self.io_context.load(load_request)


class PeptidePlannerDAO(GenericDAO):
    TYPE = 'macrocycles'

    def __init__(self, io_context):
        super().__init__(io_context)

    def load(self, load_request):
        data = super().load(load_request)
        return self.format_peptide_plan(data)

    def format_peptide_plan(self, raw_peptide_plan):
        pass


class Filesystem:

    def __init__(self, loader, saver):
        self.loader = loader
        self.saver = saver

    def load(self, load_request):
        return self.loader.load(load_request)

    def save(self, save_request):
        return self.saver.save(save_request)


class FilesystemComponent:
    def __init__(self, typeio):
        self.typeio = typeio
        self.filepaths = config.FILEPATHS


class FilesystemLoader(FilesystemComponent):

    def __init__(self, typeio, file_iterator):
        super().__init__(typeio)
        self.file_iterator = file_iterator
        self.meta_data = self.load(LoadRequest('metadata'))

    def load(self, load_request):
        filepath = self._get_filepath(load_request.data_type)
        self._initialize_file_iterator(filepath, self.meta_data, load_request.data_chunk)
        return self._load_data(load_request.data_chunk)

    def _load_data(self, data_chunk):
        for count, filepath in self.file_iterator.iterate():
            for data_unit in self.typeio.load(filepath):
                if count == data_chunk:
                    yield data_unit

    def _get_filepath(self, load_request):
        filepath = self.filepaths[load_request.data_type]
        return utils.attach_file_num(filepath, load_request.peptide_length)

    def _initialize_file_iterator(self, filepath, meta_data, data_range):
        self.file_iterator.load_files(filepath)
        self.file_iterator.skip_to(meta_data, data_range)


class FilesystemSaver(FilesystemComponent):

    def __init__(self, typeio):
        super().__init__(typeio)

    def save(self, save_request):
        data_filename, meta_data_filepath = self._get_filepaths(save_request)
        self.typeio.save(data_filename, save_request.data)
        self.typeio.save(meta_data_filepath, save_request.meta_data.to_json())

    def _get_filepaths(self, save_request):
        data_filepath = self.filepaths[save_request.data_type]
        data_filepath = utils.attach_file_num(data_filepath, save_request.peptide_length, save_request.job_num)
        _, data_filename = os.path.split(data_filepath)
        meta_data_filepath = os.path.join(config.STDOUT_DIR, data_filename)
        return data_filename, meta_data_filepath


class JsonIO:
    def load(self, filepath):
        with open(filepath, 'r') as file:
            return json_util.loads(json_util.dumps(json.load(file)))

    def save(self, filepath, data):
        with open(filepath, 'w') as file:
            json.dump(json.loads(json_util.dumps(data)), file)


class FileIterator:

    def __init__(self, glob_char='*'):
        self.filepaths = []
        self.filenames = []
        self.glob_char = glob_char
        self.count = 0

    def load_files(self, base_filepath):
        globbed_filepath = utils.attach_file_num(base_filepath, self.glob_char)
        self.filepaths = sorted(glob.glob(globbed_filepath))
        self._extract_filenames()

    def skip_to(self, meta_data, data_range):

        for i, filename in enumerate(self.filenames()):
            curr_range = Range(self.count, meta_data[filename]['count'])
            if curr_range in data_range:
                self.filepaths = self.filepaths[i:]
                self.filenames = self.filenames[i:]
                break
            self.count += meta_data[filename]['count']

    def iterate(self, meta_data):
        for filepath, filename in zip(self.filepaths, self.filenames):
            yield filepath, self.count
            self.count += meta_data[filename]['count']

    def _extract_filenames(self):
        for filepath in self.filepaths:
            _, filename = os.path.split(filepath)
            self.filenames.append(filename)


class MetaDataAggregater:
    def __init__(self, typeio):
        self.typeio = typeio
        self.input_filepaths = glob.glob(utils.attach_file_num(config.STDOUT_DIR, '*'))
        self.output_filepath = config.FILEPATHS['metadata']

    def aggregate(self):
        data = self._parse_raw_meta_data()
        self._output(data)

    def _parse_raw_meta_data(self):
        data = {}
        for filepath in self.input_filepaths:
            data.update(self.typeio.load(filepath))

        return data

    def _output(self, aggregated_data):
        stored_data = self.typeio.load(self.output_filepath)
        stored_data.update(aggregated_data)
        self.typeio.save(stored_data)


# request = LoadRequest('')
# repository.load()
