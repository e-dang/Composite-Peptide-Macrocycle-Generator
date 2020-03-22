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
