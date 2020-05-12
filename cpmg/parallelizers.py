import multiprocessing
import cpmg.config as config
import cpmg.utils as utils


class SingleProcess:
    STRING = 'single'

    def run(self, generator, data):

        for args in data:
            for result in generator.generate(*args):
                yield result


class Multiprocess:
    STRING = 'multi'

    def run(self, generator, data):

        with multiprocessing.Pool(processes=config.NUM_PROCS - 1, maxtasksperchild=config.TASKS_PER_CHILD) as pool:
            # try using callback function in data_handler.load() to calcualte the number of documents being loaded in order to calculate chunksize
            for result in pool.imap_unordered(generator.generate, data):
                for record in result:
                    yield record


get_all_parallelizer_strings = utils.get_module_strings(__name__)

create_parallelizer_from_string = utils.create_factory_function_closure(__name__, 'parallelizer')
