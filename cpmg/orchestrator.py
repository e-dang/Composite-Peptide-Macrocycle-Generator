import multiprocessing
import mpi4py  # noqa
mpi4py.rc(initialize=False, finalize=False)  # noqa
from mpi4py import MPI  # noqa
from mpi4py.futures import MPICommExecutor  # noqa

import cpmg.data_handlers as h
import cpmg.generators as g
import cpmg.ranges as r
from cpmg.config import CAPACITY
import cpmg.config as config
from cpmg.parallelism import LEVEL_0, LEVEL_1, LEVEL_2
from cpmg.timer import GlobalTimer
from cpmg.exceptions import InvalidChunkSize


class ExecutionParameters:
    def __init__(self, command_line_args):
        command_line_args.pop('func')
        command_line_args.pop('time')
        command_line_args.pop('buffer_time')
        self.operation = command_line_args.pop('operation')
        self.parallelism = command_line_args.pop('parallelism')
        self.chunk_size = command_line_args.pop('chunk_size')
        self.operation_parameters = self.__extract_operation_parameters(command_line_args)

    def __extract_operation_parameters(self, command_line_args):
        if self.operation == g.SidechainModifier.STRING:
            operation_parameters = {
                'sidechain_key': r.Key(r.WholeRange()),
                'connection_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.MonomerGenerator.STRING:
            operation_parameters = {
                'sidechain_key': r.Key(r.WholeRange()),
                'backbone_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.PeptidePlanGenerator.STRING:
            operation_parameters = {
                'monomer_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.PeptideGenerator.STRING:
            operation_parameters = {
                'peptide_plan_key': self.__generate_length_dependent_key(command_line_args),
                'monomer_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.TemplatePeptideGenerator.STRING:
            operation_parameters = {
                'peptide_key': self.__generate_length_dependent_key(command_line_args),
                'template_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.MacrocycleGenerator.STRING:
            operation_parameters = {
                'template_peptide_key': self.__generate_length_dependent_key(command_line_args),
                'reaction_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.InterMolecularReactionGenerator.STRING:
            operation_parameters = {
                'sidechain_key': r.Key(r.WholeRange()),
                'monomer_key': r.Key(r.WholeRange()),
                'template_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.IntraMolecularReactionGenerator.STRING:
            operation_parameters = {
                'template_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.ConformerGenerator.STRING:
            operation_parameters = {
                'macrocycle_key': self.__generate_length_dependent_key(command_line_args)
            }

        operation_parameters.update(command_line_args)
        return operation_parameters

    def __generate_length_dependent_key(self, command_line_args):
        num_records = command_line_args.pop('num_records', None)
        peptide_length = command_line_args.pop('peptide_length', None)
        if num_records is None:
            return r.Key(r.WholeRange(), peptide_length=peptide_length)

        return r.Key(r.DiscreteDataChunk(list(range(num_records))), peptide_length=peptide_length)


class ResultBuffer:
    def __init__(self, saver, chunk_size=None):
        self.saver = saver
        self.buffer = []
        self.ids = []
        self.chunk_size = chunk_size or CAPACITY
        self._validate_chunk_size()

    def __len__(self):
        return len(self.buffer)

    def add(self, data):
        self.buffer.append(data)
        if len(self.buffer) >= self.chunk_size:
            self.flush()

    def flush(self):
        self.ids.extend(self.saver.save(self.buffer))
        self.buffer = []

    def _validate_chunk_size(self):
        if not isinstance(self.chunk_size, int) or self.chunk_size <= 0:
            raise InvalidChunkSize('The chunk size must be a positive integer')


class Orchestrator:

    def __init__(self, impl):
        self.impl = impl

    @classmethod
    def from_strings(cls, parallelism, operation, chunk_size=None):
        generator = g.create_generator_from_string(operation)
        handler = h.create_handler_from_string(operation)
        if parallelism == SingleProcessOrchestrator.STRING:
            return cls(SingleProcessOrchestrator(generator, handler, chunk_size=chunk_size))

        if parallelism == MultiProcessOrchestrator.STRING:
            return cls(MultiProcessOrchestrator(generator, handler, chunk_size=chunk_size))

        if parallelism == DistributedOrchestrator.STRING:
            return cls(DistributedOrchestrator(generator, handler, chunk_size=chunk_size))

    @classmethod
    def from_execution_parameters(cls, params):
        return cls.from_strings(params.parallelism, params.operation, chunk_size=params.chunk_size)

    def execute(self, **operation_parameters):
        return self.impl.execute(**operation_parameters)


class AbstractOrchestratorImpl:
    def __init__(self, generator, handler, chunk_size=None):
        self.generator = generator
        self.handler = handler
        self.result_buffer = ResultBuffer(handler, chunk_size=chunk_size)
        self.timer = GlobalTimer.instance()

    def execute(self, **operation_parameters):
        pass

    def _chunkify(self, data):
        buffer = []
        for record in data:
            buffer.append(record)
            if len(buffer) == self.result_buffer.chunk_size:
                yield buffer
                buffer = []
        yield buffer


class SingleProcessOrchestrator(AbstractOrchestratorImpl):
    STRING = LEVEL_0

    def execute(self, **operation_parameters):
        for chunk in self._chunkify(self.handler.load(**operation_parameters)):
            for args in chunk:
                for record in self.generator.generate(*args):
                    self.result_buffer.add(record)

                if self.timer.is_near_complete():
                    break

        self.result_buffer.flush()
        return self.result_buffer.ids


class MultiProcessOrchestrator(AbstractOrchestratorImpl):
    STRING = LEVEL_1

    def execute(self, **operation_parameters):

        with multiprocessing.Pool(processes=config.NUM_PROCS, maxtasksperchild=config.TASKS_PER_CHILD) as pool:
            for chunk in self._chunkify(self.handler.load(**operation_parameters)):
                future = pool.starmap_async(self.generator.generate, chunk)
                for result in future.get():
                    for record in result:
                        self.result_buffer.add(record)

                    if self.timer.is_near_complete():
                        break

                break

        self.result_buffer.flush()
        return self.result_buffer.ids


class DistributedOrchestrator(AbstractOrchestratorImpl):
    STRING = LEVEL_2

    def execute(self, **operation_parameters):
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        with MPICommExecutor() as executor:
            if executor is not None:
                for chunk in self._chunkify(self.handler.load(**operation_parameters)):
                    for result in executor.starmap(self.generator.generate, chunk,
                                                   chunksize=int(len(chunk) / size - 1), unordered=True):
                        for record in result:
                            self.result_buffer.add(record)
                            break

                        if self.timer.is_near_complete():
                            break

        if MPI.COMM_WORLD.Get_rank() == 0:
            self.result_buffer.flush()

        return self.result_buffer.ids
