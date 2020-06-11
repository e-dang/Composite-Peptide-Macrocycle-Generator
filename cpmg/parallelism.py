import atexit
import mpi4py  # noqa
mpi4py.rc(initialize=False, finalize=False)  # noqa
from mpi4py import MPI  # noqa

from cpmg.exceptions import ParallelismAlreadySet

LEVEL_0 = 'single'
LEVEL_1 = 'multi'
LEVEL_2 = 'distributed'


def get_parallelism_strings():
    return (LEVEL_0, LEVEL_1, LEVEL_2)


class Parallelism:
    __LEVEL = None
    __LEVELS = get_parallelism_strings()

    @classmethod
    def get_level(cls):
        return cls.__LEVEL

    @classmethod
    def is_distributed(cls):
        return cls.__LEVEL == LEVEL_2

    @classmethod
    def set_level(cls, level):
        if cls.__LEVEL != None:
            raise ParallelismAlreadySet('Can not reset the parallelism level when it has already been set.')
        if level not in cls.__LEVELS:
            raise ValueError(f'Unrecognized parallelism option! Valid choices are {cls.__LEVELS}')

        cls.__LEVEL = level
        if level == LEVEL_2:
            MPI.Init_thread()
            atexit.register(MPI.Finalize)
