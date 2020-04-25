import math
from collections import deque
import new_architecture.models as models
from new_architecture.repository.hdf5 import HDF5Repository, to_list
import macrocycles.config as config
import random


HDF5 = 'hdf5'


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


class AbstractRepository:
    TYPE = None
    CATEGORY = None

    def __init__(self, impl):
        self.impl = impl

    def load(self, key):
        for _id, data in self.impl.load(self.CATEGORY, key):
            yield self.TYPE.from_dict(data, _id=_id)

    def save(self, data):
        if not isinstance(random.choice(data), self.TYPE):
            print(self.TYPE, random.choice(data))
            raise TypeError('Error! Repository - type mismatch')

        return self.impl.save(self.CATEGORY, map(lambda x: x.to_dict(), data))


class TemplateRepository(AbstractRepository):
    TYPE = models.Template
    CATEGORY = 'templates'

    def __init__(self, impl):
        super().__init__(impl)


class SidechainRepository(AbstractRepository):
    TYPE = models.Sidechain
    CATEGORY = 'sidechains'

    def __init__(self, impl):
        super().__init__(impl)


class MonomerRepository(AbstractRepository):
    TYPE = models.Monomer
    CATEGORY = 'monomers'

    def __init__(self, impl):
        super().__init__(impl)


class PeptideRepository(AbstractRepository):
    TYPE = models.Peptide
    CATEGORY = 'peptides'

    def __init__(self, impl):
        super().__init__(impl)


def repository_impl_from_string(impl):
    if impl == HDF5:
        return HDF5Repository()
    else:
        raise ValueError('Unrecognized repository implementation')


def get_repository(repository):
    def repository_closure(impl=config.DATA_FORMAT):
        return repository(repository_impl_from_string(impl))

    return repository_closure


create_template_repository = get_repository(TemplateRepository)
create_sidechain_repository = get_repository(SidechainRepository)
create_monomer_repository = get_repository(MonomerRepository)
create_peptide_repository = get_repository(PeptideRepository)

# def create_sidechain_repository(impl=config.DATA_FORMAT):
#     return SidechainRepository(repository_impl_from_string(impl))


# class SaveRequest:
#     def __init__(self, data, data_type, peptide_length=None, job_num=None):
#         self.data = data
#         self.data_type = data_type
#         self.peptide_length = peptide_length
#         self.job_num = job_num

#     @classmethod
#     def create_request(cls, data, meta_data):
#         return cls(data, meta_data.data_type, meta_data.peptide_length, meta_data.job_num)


# class LoadRequest:

#     def __init__(self, data_type, peptide_length=None, data_chunk=WholeRange()):
#         self.data_type = data_type
#         self.peptide_length = peptide_length
#         self.data_chunk = data_chunk
#         self.next_request = None

#     def add_request(self, request):
#         if self.next_request is None:
#             self.next_request = request
#         else:
#             self.next_request.add_request(request)

#     def get_next_request(self, result):
#         try:
#             self.next_request.update_data_chunk(result)
#         except AttributeError:
#             pass

#         return self.next_request

#     def update_data_chunk(self, data_chunk):
#         self.data_chunk = data_chunk

#     def format_result(self, func):
#         self.result = func(self.result)

# def create_repository():
#     import sys
#     import inspect
#     from collections import defaultdict

#     clsmembers = inspect.getmembers(sys.modules[__name__], lambda x: inspect.isclass(x) and 'DAO' in x.__name__)

#     daos = defaultdict(GenericDAO)
#     for _, instance in clsmembers:
#         try:
#             daos[instance.TYPE] = instance
#         except AttributeError:
#             continue

#     return Repository(daos)


# repository = create_repository()
