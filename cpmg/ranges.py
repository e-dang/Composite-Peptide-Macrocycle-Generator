import math


class Range:
    """
    Defines a range of numbers in interval [start, end).
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __str__(self):
        return f'start: {self.start}\nend: {self.end}\n'

    def __contains__(self, other):
        return (self.start < other.end and other.start < self.end) or (other.start < self.end and self.start < other.end)

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __len__(self):
        return self.end - self.start

    def __getitem__(self, idx):
        return self.start + idx


class WholeRange(Range):

    def __init__(self):
        super().__init__(float('-inf'), float('inf'))

    def __contains__(self, other):
        return True

    def __eq__(self, other):
        return isinstance(other, WholeRange)

    def __len__(self):
        return self.end

    def __getitem__(self, idx):
        return idx


class DiscreteDataChunk(Range):
    def __init__(self, indices):
        self.indices = sorted(indices)
        super().__init__(self.indices[0], self.indices[-1])

    def __str__(self):
        return str(self.indices)

    def __contains__(self, number):
        return number in self.indices

    def __eq__(self, other):
        return self.indices == other.indices

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        return self.indices[idx]


class BatchedDataChunk(Range):

    def __init__(self, total_num_data, total_num_jobs, job_num):
        self.calc_chunk_size(total_num_data, total_num_jobs)
        super().__init__(self.calc_start(job_num), self.calc_end(total_num_data, job_num))

    def __str__(self):
        return super().__str__() + f'chunk_size: {self.chunk_size}\n'

    def __len__(self):
        return self.chunk_size

    def calc_chunk_size(self, total_num_data, total_num_jobs):
        self.chunk_size = math.ceil(total_num_data / total_num_jobs)

    def calc_start(self, job_num):
        return self.chunk_size * (job_num - 1)

    def calc_end(self, total_num_data, job_num):
        end = self.chunk_size * job_num

        if end > total_num_data:
            end = total_num_data

        return end


class RangeKey:
    def __init__(self, range, peptide_length=None):
        self.range = range
        self.peptide_length = peptide_length

    def __str__(self):
        return self.range.__str__() + f'peptide_length: {self.peptide_length}\n'

    def __contains__(self, other):
        return other in self.range

    def __eq__(self, other):
        try:
            return self.range == other.range and self.peptide_length == self.peptide_length
        except AttributeError:
            return False

    def __len__(self):
        return len(self.range)

    def __getitem__(self, idx):
        return self.range[idx]

    @property
    def start(self):
        return self.range.start

    @property
    def end(self):
        return self.range.end


class IndexKey:
    def __init__(self, ids, peptide_length=None, index='id'):
        self.ids = sorted(ids)
        self.peptide_length = peptide_length
        self.index = index

    def __str__(self):
        return f'ids: {self.ids}\npeptide_length: {self.peptide_length}\n'

    def __contains__(self, _id):
        return _id in self.ids

    def __eq__(self, other):
        try:
            return self.ids == other.ids
        except AttributeError:
            return self.ids == other

    def __iter__(self):
        return iter(self.ids)

    def __next__(self):
        return next(self.ids)

    def __len__(self):
        return len(self.ids)

    def __getitem__(self, idx):
        return self.ids[idx]

    @property
    def start(self):
        return self.ids[0]

    @property
    def end(self):
        return self.ids[-1]


class Key:
    def __init__(self, identifiers, peptide_length=None, index='id'):
        if isinstance(identifiers, Range):
            self.key = RangeKey(identifiers, peptide_length=peptide_length)
        else:
            self.key = IndexKey(identifiers, peptide_length=peptide_length, index=index)

    def __str__(self):
        return self.key.__str__()

    def __contains__(self, identifier):
        return identifier in self.key

    def __eq__(self, other):
        return self.key == other

    def __iter__(self):
        return iter(self.key)

    def __next__(self):
        return next(self.key)

    def __len__(self):
        return len(self.key)

    def __getitem__(self, idx):
        return self.key[idx]

    def set_peptide_length(self, length):
        self.key.peptide_length = length

    @property
    def peptide_length(self):
        return self.key.peptide_length

    @property
    def index(self):
        try:
            return self.key.index
        except AttributeError:
            return None

    @property
    def start(self):
        return self.key.start

    @property
    def end(self):
        return self.key.end
