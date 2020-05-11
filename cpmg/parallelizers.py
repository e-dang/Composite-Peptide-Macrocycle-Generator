import multiprocessing
import cpmg.config as config
import cpmg.utils as utils


class Parallelizer:
    def __init__(self):
        self.results = []

    def save_result(self, result, data_handler):
        self.results.append(result)

        if len(self.results) == config.CAPACITY:
            self._flush(data_handler)

    def execute(self, generator, data_handler):
        self._execute(generator, data_handler)
        self._flush(data_handler)

    def _execute(self, generator, data_handler):
        pass

    def _flush(self, data_handler):
        data_handler.save(self.results)
        self.results = []


class SingleProcess(Parallelizer):
    STRING = 'single'

    def _execute(self, generator, data_handler):

        for args in data_handler.load():
            for result in generator.generate(args):
                self.save_result(result, data_handler)


class Multiprocess(Parallelizer):
    STRING = 'multi'

    def _execute(self, generator, data_handler):

        with multiprocessing.Pool(processes=config.NUM_PROCS - 1, maxtasksperchild=config.TASKS_PER_CHILD) as pool:
            # try using callback function in data_handler.load() to calcualte the number of documents being loaded in order to calculate chunksize
            for result in pool.imap_unordered(generator.generate, data_handler.load()):
                self.save_result(result, data_handler)


get_all_parallelizer_strings = utils.get_module_strings(__name__)

create_parallelizer_from_string = utils.create_factory_function_closure(__name__, 'parallelizer')
# def parallelizer_from_string(string):
#     parallelizer_map = {
#         SingleProcess.STRING: SingleProcess(),
#         Multiprocess.STRING: Multiprocess()
#     }

#     try:
#         return parallelizer_map[string]
#     except KeyError:
#         raise KeyError('Unrecognized parallelizer string!')
