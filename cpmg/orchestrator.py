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
import cpmg.utils as utils
from cpmg.parallelism import LEVEL_0, LEVEL_1, LEVEL_2


class ExecutionParameters:
    def __init__(self, command_line_args):
        command_line_args.pop('func')
        self.operation = command_line_args.pop('operation')
        self.parallelism = command_line_args.pop('parallelism')
        self.chunk_size = command_line_args.pop('chunk_size')
        self.operation_parameters = self._extract_operation_parameters(command_line_args)

    def _extract_operation_parameters(self, command_line_args):
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
                'peptide_plan_key': r.Key(r.WholeRange(), peptide_length=command_line_args.pop('peptide_length', None)),
                'monomer_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.TemplatePeptideGenerator.STRING:
            operation_parameters = {
                'peptide_key': r.Key(r.WholeRange(), peptide_length=command_line_args.pop('peptide_length', None)),
                'template_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.MacrocycleGenerator.STRING:
            operation_parameters = {
                'template_peptide_key': r.Key(r.WholeRange(), peptide_length=command_line_args.pop('peptide_length', None)),
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
                'macrocycle_key': r.Key(r.WholeRange(), peptide_length=command_line_args.pop('peptide_length', None)),
            }

        operation_parameters.update(command_line_args)
        return operation_parameters


class ResultBuffer:
    def __init__(self, saver, chunk_size=None):
        self.saver = saver
        self.chunk_size = chunk_size or CAPACITY
        self.buffer = []
        self.ids = []

    def __len__(self):
        return len(self.buffer)

    def add(self, data):
        self.buffer.append(data)
        if len(self.buffer) >= self.chunk_size:
            self.flush()

    def flush(self):
        self.ids.extend(self.saver.save(self.buffer))
        self.buffer = []


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

    def execute(self, **operation_parameters):
        pass


class SingleProcessOrchestrator(AbstractOrchestratorImpl):
    STRING = LEVEL_0

    def execute(self, **operation_parameters):

        for args in self.handler.load(**operation_parameters):
            for record in self.generator.generate(*args):
                self.result_buffer.add(record)

        self.result_buffer.flush()
        return self.result_buffer.ids


class MultiProcessOrchestrator(AbstractOrchestratorImpl):
    STRING = LEVEL_1

    def execute(self, **operation_parameters):

        with multiprocessing.Pool(processes=config.NUM_PROCS - 1, maxtasksperchild=config.TASKS_PER_CHILD) as pool:
            # try using callback function in data_handler.load() to calcualte the number of documents being loaded in order to calculate chunksize
            future = pool.starmap_async(self.generator.generate, self.handler.load(**operation_parameters))
            for result in future.get():
                for record in result:
                    self.result_buffer.add(record)

        self.result_buffer.flush()
        return self.result_buffer.ids


class DistributedOrchestrator(AbstractOrchestratorImpl):
    STRING = LEVEL_2

    def execute(self, **operation_parameters):
        with MPICommExecutor() as executor:
            if executor is not None:
                for args in self.handler.load(**operation_parameters):
                    future = executor.submit(self.generator.generate, *args)
                    for record in future.result():
                        self.result_buffer.add(record)

        if MPI.COMM_WORLD.Get_rank() == 0:
            self.result_buffer.flush()

        return self.result_buffer.ids
