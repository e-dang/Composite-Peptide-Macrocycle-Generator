import math
from collections import deque


class Range:
    """
    Defines a range of numbers in interval [start, end).
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __contains__(self, other):
        return self.start < other.end or other.start < self.end


class WholeRange(Range):

    def __init__(self):
        super().__init__(None, None)

    def __contains__(self, other):
        return True

    def __eq__(self, other):
        return True


class DiscreteDataChunk(Range):
    def __init__(self, indices):
        self.indices = deque(sorted(indices))
        super().__init__(min(self.indices), max(self.indices))

    def __eq__(self, number):
        try:
            if self.indices[0] == number:
                self.pop()
                return True
        except IndexError:
            raise StopIteration

        return False

    def pop(self):
        self.indices.popleft()


class ContinuousDataChunk(Range):

    def __init__(self, total_num_data, total_num_jobs, job_num):
        self.calc_chunk_size(total_num_data, total_num_jobs)
        super().__init__(self.calc_start(job_num), self.calc_end(total_num_data, job_num))

    def __eq__(self, number):
        if self.start <= number < self.end:
            return True
        elif self.end <= number:
            raise StopIteration

        return False

    def calc_chunk_size(self, total_num_data, total_num_jobs):
        self.chunk_size = math.ceil(total_num_data / total_num_jobs)

    def calc_start(self, job_num):
        return self.chunk_size * (job_num - 1)

    def calc_end(self, total_num_data, job_num):
        end = self.chunk_size * job_num

        if end > total_num_data:
            end = total_num_data

        return end


class MetaData:
    def __init__(self, data_type, peptide_length, job_num):
        self.data_type = data_type
        self.peptide_length = peptide_length
        self.job_num = job_num

    def to_json(self, filename):
        json_data = {filename: self.__dict__}
        return json_data


class SaveRequest:
    def __init__(self, data, data_type, peptide_length=None, job_num=None):
        self.data = data
        self.data_type = data_type
        self.peptide_length = peptide_length
        self.job_num = job_num

    @classmethod
    def create_request(cls, data, meta_data):
        return cls(data, meta_data.data_type, meta_data.peptide_length, meta_data.job_num)


class LoadRequest:

    def __init__(self, data_type, peptide_length=None, data_chunk=WholeRange()):
        self.data_type = data_type
        self.peptide_length = peptide_length
        self.data_chunk = data_chunk
        self.next_request = None

    def add_request(self, request):
        if self.next_request is None:
            self.next_request = request
        else:
            self.next_request.add_request(request)

    def get_next_request(self, result):
        try:
            self.next_request.update_data_chunk(result)
        except AttributeError:
            pass

        return self.next_request

    def update_data_chunk(self, data_chunk):
        self.data_chunk = data_chunk

    def format_result(self, func):
        self.result = func(self.result)


class Repository:

    def __init__(self, daos):
        self.daos = daos
        self.data = []
        self.meta_data = None

    def add(self, data):
        self.data.append(data)

    def add_meta_data(self, meta_data):
        self.meta_data = meta_data

    def load(self, load_request):
        while True:
            data = self.daos[load_request.data_type].load(load_request)
            load_request = load_request.get_next_request(data)
            if load_request is None:
                return data

    def save(self):
        self.daos[self.meta_data.data_type].save(SaveRequest.create_request(self.data, self.meta_data))


def create_repository():
    import sys
    import inspect
    from collections import defaultdict

    clsmembers = inspect.getmembers(sys.modules[__name__], lambda x: inspect.isclass(x) and 'DAO' in x.__name__)

    daos = defaultdict(GenericDAO)
    for _, instance in clsmembers:
        try:
            daos[instance.TYPE] = instance
        except AttributeError:
            continue

    return Repository(daos)


repository = create_repository()
