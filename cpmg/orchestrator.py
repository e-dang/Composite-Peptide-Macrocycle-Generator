import cpmg.generators as g
import cpmg.parallelizers as p
import cpmg.data_handlers as h
import cpmg.ranges as r
from cpmg.config import CAPACITY


class ExecutionParameters:
    def __init__(self, command_line_args):
        command_line_args.pop('func')
        self.operation = command_line_args.pop('operation')
        self.parallelism = command_line_args.pop('parallelism')
        self.operation_parameters = self._extract_operation_parameters(command_line_args)

    def _extract_operation_parameters(self, command_line_args):
        if self.operation == g.PeptideGenerator.STRING:
            operation_parameters = {
                'peptide_plan_key': r.Key(r.WholeRange(), peptide_length=command_line_args.pop('peptide_length', None)),
                'monomer_key': r.Key(r.WholeRange())
            }
        elif self.operation == g.TemplatePeptideGenerator.STRING:
            operation_parameters = {
                'peptide_key': r.Key(r.WholeRange(), peptide_length=command_line_args.pop('peptide_length', None)),
                'template_key': r.Key(r.WholeRange())
            }
        else:
            operation_parameters = {
                'key': r.Key(r.WholeRange())
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


class Orchestractor:

    def __init__(self, generator, handler, parallelizer, chunk_size=None):
        self.generator = generator
        self.handler = handler
        self.parallelizer = parallelizer
        self.result_buffer = ResultBuffer(handler, chunk_size=chunk_size)

    @classmethod
    def from_execution_parameters(cls, params):
        generator = g.create_generator_from_string(params.operation)
        handler = h.create_handler_from_string(params.operation)
        parallelizer = p.create_parallelizer_from_string(params.parallelism)

        return cls(generator, handler, parallelizer)

    def execute(self, **operation_parameters):

        for record in self.parallelizer.run(self.generator, self.handler.load(**operation_parameters)):
            self.result_buffer.add(record)

        self.result_buffer.flush()
        return self.result_buffer.ids
