from time import perf_counter

from cpmg.exceptions import TimerAlreadyStarted, TimerNotStarted


class GlobalTimer:

    def __init__(self, time_allocation=None, buffer_time=300):  # buffer time is 5 minutes
        self.__start_time = None
        self.time_allocation = time_allocation or -1
        self.buffer_time = buffer_time

    @classmethod
    def instance(cls, time_allocation=None, buffer_time=300):
        if not hasattr(cls, '_instance'):
            cls._instance = cls(time_allocation=time_allocation, buffer_time=buffer_time)
        return cls._instance

    def has_started(self):
        return self.__start_time is not None

    def start(self):
        if self.__start_time is not None:
            raise TimerAlreadyStarted('Timer has already been started')

        self.__start_time = perf_counter()

    def elapsed_time(self):
        if self.__start_time is None:
            raise TimerNotStarted('Timer has not been started yet')

        return perf_counter() - self.__start_time

    def is_near_complete(self):
        if self.time_allocation < 0:
            return False

        return self.time_allocation - self.elapsed_time() <= self.buffer_time
