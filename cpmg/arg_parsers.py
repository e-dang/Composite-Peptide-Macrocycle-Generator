import argparse
import sys
import random

from rdkit import Chem

import cpmg.parallelism as p
import cpmg.models as m
import cpmg.generators as g
import cpmg.orchestrator as o
import cpmg.operations as ops
import cpmg.importers as i
import cpmg.repository as r
import cpmg.data_handlers as h
import cpmg.utils as utils
import cpmg.ranges as ranges
from cpmg.timer import GlobalTimer
from cpmg.initializer import CPMGInitializer
from cpmg.exporters import RegioSQMExporter


def check_positive(arg):
    int_value = int(arg)
    if int_value <= 0:
        raise argparse.ArgumentTypeError(f'{arg} is an invalid positive int value')
    return int_value


def add_length_option(parser):
    parser.add_argument('-l', '--length', '--peptide_length', type=int, choices=[3, 4, 5],
                        required=True, dest='peptide_length',
                        help='The length of the peptide in the desired records.')


def add_number_option(parser):
    parser.add_argument('-n', '--num', '--num_records', type=check_positive, dest='num_records',
                        help='The number of records to load or create for the given operation.')


def add_selection_options(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-i', '--ids', type=str, nargs='+',
                       help='The ids of the records that will be used in the operation.')
    group.add_argument('-s', '--smiles', type=str, nargs='+',
                       help='The SMILES strings of the records that will be used in the operation.')
    group.add_argument('-r', '--rand', '--random_sample', type=int, dest='rand',
                       help='Causes the specified number of records to be randomly selected for the operation.')
    # group.add_argument('--in_order', type=int,
    #                    help='Causes the specified number of records to be loaded in order for the operation.')
    # parser.add_argument('-t', '--type', '--selection_type', type=str, choices=['ids', 'smiles', 'rand', 'in_order'],
    #                     nargs='?', default='in_order', dest='selection_type',
    #                     help='The method to use to select the desired records. ids = list of ids, smiles = list of '
    #                     'SMILES strings, rand = random sampling (requires --num_records), in_order = loads the records '
    #                     'in the order they were stored (optionally takes --num_records). Defaults to in_order.')
    parser.add_argument('-n', '--num', '--num_records', type=check_positive,
                        required='--rand' in sys.argv or '-r' in sys.argv, dest='num_records',
                        help='The number of records to load. Invalid option if --selection_type is "ids" or "smiles". '
                        'Required if --selection_type is "rand". Optional if --selection_type is "in_order".')


def add_parallelism_options(parser):
    parser.add_argument('-p', '--parallelism', type=str, choices=p.get_parallelism_strings(), nargs='?',
                        default='single', help='The level of parallelism to use when executing the specified operation.')
    parser.add_argument('-c', '--chunk', '--chunk_size', type=int, default=None, dest='chunk_size',
                        help='The number of results to buffer before saving them. Default is specified in config.py.')


def add_time_options(parser):
    parser.add_argument('-t', '--time', type=int, default=None, const=None, nargs='?',
                        help='The total time allocated to the job in seconds. Defaults to -1 (unlimited time).')
    parser.add_argument('-b', '--buffer_time', type=int, default=300, const=300, nargs='?',
                        help='The amount of time before the maximum job time is reached that check pointing routines should be initiated. Defaults to 300 seconds.')


def add_file_options(parser):
    parser.add_argument('-f', '--file', '--filepath', type=str, dest='filepath',
                        help='The path to the file containing the data to import. Only valid for importing peptide plans.')


def add_repository_argument(parser, default=False):
    if default:
        parser.add_argument('repository', type=str, choices=r.get_all_repository_strings(), nargs='?',
                            default=r.CPMGRepository.STRING,
                            help='The repository that contains the desired types of records. Default repository is all.')
    else:
        parser.add_argument('repository', type=str,
                            help='The repository that contains the desired types of records.')


def add_record_type_subparser(parser):
    return parser.add_subparsers(dest='record_type', help='The type of record to use for the selected operation.')


def add_component_subparser(parser):
    return parser.add_subparsers(dest='component_type', help='The type of record(s) that are components of the selected type for the selected operation.')


def add_connection_parser(parser, selection=False):
    connection_parser = parser.add_parser(m.Connection.STRING)

    if selection:
        add_selection_options(connection_parser)

    return connection_parser


def add_backbone_parser(parser, selection=False):
    backbone_parser = parser.add_parser(m.Backbone.STRING)

    if selection:
        add_selection_options(backbone_parser)

    return backbone_parser


def add_template_parser(parser, selection=False):
    template_parser = parser.add_parser(m.Template.STRING)

    if selection:
        add_selection_options(template_parser)

    return template_parser


def add_sidechain_parser(parser, selection=False):
    sidechain_parser = parser.add_parser(g.SidechainModifier.STRING)

    if selection:
        add_selection_options(sidechain_parser)

    return sidechain_parser


def add_monomer_parser(parser, selection=False):
    monomer_parser = parser.add_parser(g.MonomerGenerator.STRING)

    if selection:
        add_selection_options(monomer_parser)

    return monomer_parser


def add_peptide_plan_parser(parser, selection=False, length=True, number=True):
    peptide_plan_parser = parser.add_parser(g.PeptidePlanGenerator.STRING)

    if selection:
        add_selection_options(peptide_plan_parser)

    if length:
        add_length_option(peptide_plan_parser)

    if number:
        add_number_option(peptide_plan_parser)

    return peptide_plan_parser


def add_peptide_parser(parser, selection=False, length=True):
    peptide_parser = parser.add_parser(g.PeptideGenerator.STRING)

    if selection:
        add_selection_options(peptide_parser)

    if length:
        add_length_option(peptide_parser)

    return peptide_parser


def add_template_peptide_parser(parser, selection=False, length=True):
    template_peptide_parser = parser.add_parser(g.TemplatePeptideGenerator.STRING)

    if selection:
        add_selection_options(template_peptide_parser)

    if length:
        add_length_option(template_peptide_parser)

    return template_peptide_parser


def add_macrocycle_parser(parser, selection=False, length=True):
    macrocycle_parser = parser.add_parser(g.MacrocycleGenerator.STRING)

    if selection:
        add_selection_options(macrocycle_parser)

    if length:
        add_length_option(macrocycle_parser)

    return macrocycle_parser


def add_conformer_parser(parser, selection=False, length=True):
    conformer_parser = parser.add_parser(g.ConformerGenerator.STRING)

    if selection:
        add_selection_options(conformer_parser)

    if length:
        add_length_option(conformer_parser)

    return conformer_parser


def add_inter_reaction_parser(parser, selection=False):
    inter_reaction = parser.add_parser(g.InterMolecularReactionGenerator.STRING)

    if selection:
        add_selection_options(inter_reaction)

    return inter_reaction


def add_intra_reaction_parser(parser, selection=False):
    intra_reaction = parser.add_parser(g.IntraMolecularReactionGenerator.STRING)

    if selection:
        add_selection_options(intra_reaction)

    return intra_reaction


def create_selection_key(command_line_args):
    peptide_length = command_line_args.peptide_length

    if command_line_args.ids is not None:
        return r.Key(command_line_args.ids, peptide_length=peptide_length)

    if command_line_args.smiles is not None:
        return r.Key(command_line_args.smiles, peptide_length=peptide_length, index='kekule')

    if command_line_args.rand is not None:
        repo = r.create_repository_from_string(command_line_args.record_type)
        indices = list(random.sample(range(repo.get_num_records()), command_line_args.num_records))
        return r.Key(indices, peptide_length=peptide_length)

    if command_line_args.num_records is not None:
        return r.Key(ranges.DiscreteDataChunk(list(range(command_line_args.num_records))), peptide_length=peptide_length)

    return r.Key(r.WholeRange(), peptide_length=peptide_length)


class AbstractParser:
    def __init__(self, args, func):
        self.args = args
        self.parser = argparse.ArgumentParser()
        self.parser.set_defaults(func=func)

    def execute(self):
        args = self.parser.parse_args(self.args)
        return args.func(args)


class GenerateParser(AbstractParser):
    def __init__(self, args):
        super().__init__(args, self.__execute)
        add_parallelism_options(self.parser)
        add_time_options(self.parser)

        subparser = add_record_type_subparser(self.parser)
        add_sidechain_parser(subparser)
        add_monomer_parser(subparser)
        add_peptide_plan_parser(subparser)
        add_peptide_parser(subparser)
        add_template_peptide_parser(subparser)
        add_macrocycle_parser(subparser)
        add_conformer_parser(subparser)
        add_inter_reaction_parser(subparser)
        add_intra_reaction_parser(subparser)

    @staticmethod
    def __execute(command_line_args):
        timer = GlobalTimer.instance(time_allocation=command_line_args.time, buffer_time=command_line_args.buffer_time)
        timer.start()
        if timer.is_near_complete():
            return []

        p.Parallelism.set_level(command_line_args.parallelism)
        params = o.ExecutionParameters(**vars(command_line_args))
        orchestrator = o.Orchestrator.from_execution_parameters(params)
        return orchestrator.execute(**params.operation_parameters)


class CalculateParser(AbstractParser):
    OPERATIONS = [ops.CalcMW.STRING, ops.CalcRB.STRING, ops.CalcTPSA.STRING]

    def __init__(self, args):
        super().__init__(args, self.__execute)
        self.parser.add_argument('operation', type=str, choices=self.OPERATIONS)
        add_parallelism_options(self.parser)
        add_time_options(self.parser)
        add_file_options(self.parser)

        subparser = add_record_type_subparser(self.parser)
        add_sidechain_parser(subparser, selection=True)
        add_monomer_parser(subparser, selection=True)
        add_peptide_parser(subparser, selection=True)
        add_template_peptide_parser(subparser, selection=True)
        add_macrocycle_parser(subparser, selection=True)
        add_conformer_parser(subparser, selection=True)
        # self.parser.add_argument('record_type', type=str, choices=[m.Connection.STRING, m.Template.STRING,
        #                                                            m.Backbone.STRING, m.Sidechain.STRING,
        #                                                            m.Monomer.STRING, m.Peptide.STRING,
        #                                                            m.TemplatePeptide.STRING, m.Macrocycle.STRING,
        #                                                            m.Conformer.STRING],
        #                          help='The type of record to perform the calculation on.')

    @staticmethod
    def __execute(command_line_args):
        timer = GlobalTimer.instance(time_allocation=command_line_args.time, buffer_time=command_line_args.buffer_time)
        timer.start()
        if timer.is_near_complete():
            return []

        p.Parallelism.set_level(command_line_args.parallelism)
        # params = o.ExecutionParameters(**vars(command_line_args))
        key = create_selection_key(command_line_args)
        generator = g.GeneratorWrapper(ops.get_calculation_from_string(command_line_args.operation))
        handler = h.DataHandlerWrapper(r.create_repository_from_string(
            command_line_args.record_type), utils.save_json, command_line_args.filepath)
        orchestrator = o.Orchestrator.from_objects(
            command_line_args.parallelism, generator, handler, chunk_size=command_line_args.chunk_size)
        return orchestrator.execute(key=key)


class ImportParser(AbstractParser):
    def __init__(self, args):
        super().__init__(args, self.__execute)
        self.parser.add_argument('record_type', type=str, choices=[m.Backbone.STRING, m.Connection.STRING,
                                                                   m.Template.STRING, m.Sidechain.STRING,
                                                                   m.Monomer.STRING, m.RegioSQMPrediction.STRING,
                                                                   m.pKaPrediction.STRING, m.PeptidePlan.STRING],
                                 help='The type of record to import.')
        add_file_options(self.parser)
        if m.PeptidePlan.STRING in args:
            add_length_option(self.parser)

    @staticmethod
    def __execute(command_line_args):
        command_line_args = vars(command_line_args)
        del command_line_args['func']
        record_type = command_line_args.pop('record_type')

        if record_type == m.Backbone.STRING:
            importer = i.BackboneImporter(i.JsonImporter())
        elif record_type == m.Connection.STRING:
            importer = i.ConnectionImporter(i.JsonImporter())
        elif record_type == m.Template.STRING:
            importer = i.TemplateImporter(i.JsonImporter())
        elif record_type == m.Sidechain.STRING:
            importer = i.SidechainImporter(i.JsonImporter())
        elif record_type == m.Monomer.STRING:
            importer = i.MonomerImporter(i.JsonImporter())
        elif record_type == m.RegioSQMPrediction.STRING:
            importer = i.RegioSQMPredictionImporter()
        elif record_type == m.pKaPrediction.STRING:
            importer = i.pKaPredictionImporter()
        elif record_type == m.PeptidePlan.STRING:
            importer = i.PeptidePlanImporter()

        return importer.import_data(**command_line_args)


class ExportParser(AbstractParser):
    def __init__(self, args):
        super().__init__(args, self.__execute)
        self.parser.add_argument('record_type', type=str, choices=[m.RegioSQMPrediction.STRING])

    @staticmethod
    def __execute(command_line_args):
        if command_line_args.record_type == m.RegioSQMPrediction.STRING:
            return RegioSQMExporter().export_regiosqm_smiles_file()


class InitializeParser(AbstractParser):
    RECORD_TYPES = [CPMGInitializer.ALL, CPMGInitializer.MOL_TYPE, CPMGInitializer.PREDICTION_TYPE]

    def __init__(self, args):
        super().__init__(args, self.__execute)
        self.parser.add_argument('record_types', type=str, choices=self.RECORD_TYPES, nargs='?',
                                 default=CPMGInitializer.ALL,
                                 help='The types of records to initialize in the repository.')

    @staticmethod
    def __execute(command_line_args):
        return CPMGInitializer().initialize(command_line_args.record_types)


class PrintParser(AbstractParser):
    def __init__(self, args):
        super().__init__(args, self.__execute)
        add_repository_argument(self.parser, default=True)

    @staticmethod
    def __execute(command_line_args):
        print(r.create_repository_from_string(command_line_args.repository))


class FindParser(AbstractParser):
    def __init__(self, args):
        super().__init__(args, self.__execute)
        add_repository_argument(self.parser)
        self.parser.add_argument('-a', '--all', action='store_true',
                                 help='Flag that causes all records in the specified repo to be loaded.')
        self.parser.add_argument('-k', '--kekule', type=str, nargs='+',
                                 help='Load and print the records corresponding to the specified kekule strings.')
        self.parser.add_argument('-f', '--filter', type=str,
                                 help='The filter criteria as a python dictionary, where the key is the field on the record and the value is the value of that field in desired record(s). File must be a .sdf file')
        self.parser.add_argument('-p', '--projection', type=str, nargs='+',
                                 help='Only print the specified field from the retrieved records.')
        self.parser.add_argument('-i', '--inactives', action='store_true',
                                 help='Flag that selects inactive records rather than active records in the given repository that meet the search criteria.')
        self.parser.add_argument('-c', '--chunk', '--chunk_size', type=int, const=10, default=10, nargs='?',
                                 dest='chunk_size',
                                 help='The number of records to display on screen at a single time.')

    @staticmethod
    def __execute(command_line_args):
        finder = ops.Finder()
        return finder.execute(command_line_args)


class RemoveParser(AbstractParser):
    def __init__(self, args):
        super().__init__(args, self.__execute)
        add_repository_argument(self.parser)
        self.parser.add_argument('-i', '--ids', type=str, nargs='+', help='The ids of the records to be removed.')

    @staticmethod
    def __execute(command_line_args):
        repo = r.create_repository_from_string(command_line_args.repository)
        return repo.remove_records(r.Key(command_line_args.ids))


class DeactivateParser(AbstractParser):
    def __init__(self, args):
        super().__init__(args, self.__execute)
        add_repository_argument(self.parser)

        file_deactivate = self.parser.add_argument_group('from_file')
        add_file_options(file_deactivate)
        file_deactivate.add_argument('-i', '--invert', action='store_true',
                                     help='Causes the structures in the file to remain active while all records not matched to structures in the file are deactivated.')

        completed_deactivate = self.parser.add_argument_group('completed')
        completed_deactivate.add_argument('-c', '--completed', action='store_true',
                                          help='Causes the completed records to be deactivated.')

    @staticmethod
    def __execute(command_line_args):
        repo = r.create_repository_from_string(command_line_args.repository)
        if command_line_args.filepath is not None:
            return DeactivateParser.__execute_deactivate(repo, command_line_args)

        if command_line_args.completed is not None:
            return repo.deactivate_completed()

    @staticmethod
    def __execute_deactivate(repo, command_line_args):
        repo_mols = {mol.kekule: mol._id for mol in repo.load()}

        file_mols = []
        reader = Chem.SDMolSupplier(command_line_args.filepath)
        for mol in reader:
            mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
            Chem.Kekulize(mol)
            file_mols.append(Chem.MolToSmiles(mol, kekuleSmiles=True))

        if command_line_args.invert:
            deactivated_mols = set(repo_mols.keys()).difference(file_mols)
        else:
            deactivated_mols = file_mols

        deactivated_ids = [repo_mols[kekule] for kekule in deactivated_mols]
        return repo.deactivate_records(r.Key(deactivated_ids))


class CPMGParser:
    GENERATE = 'generate'
    CALCULATE = 'calculate'
    IMPORT = 'import'
    EXPORT = 'export'
    INITIALIZE = 'initialize'
    PRINT = 'print'
    FIND = 'find'
    REMOVE = 'remove'
    DEACTIVATE = 'deactivate'
    OPERATIONS = [GENERATE, CALCULATE, IMPORT, EXPORT, INITIALIZE, PRINT, FIND, REMOVE, DEACTIVATE]

    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('operation', type=str, choices=self.OPERATIONS)

    def execute(self):
        args, extra = self.parser.parse_known_args()

        if args.operation == self.GENERATE:
            return GenerateParser(extra).execute()

        if args.operation == self.CALCULATE:
            return CalculateParser(extra).execute()

        if args.operation == self.IMPORT:
            return ImportParser(extra).execute()

        if args.operation == self.EXPORT:
            return ExportParser(extra).execute()

        if args.operation == self.INITIALIZE:
            return InitializeParser(extra).execute()

        if args.operation == self.PRINT:
            return PrintParser(extra).execute()

        if args.operation == self.FIND:
            return FindParser(extra).execute()

        if args.operation == self.REMOVE:
            return RemoveParser(extra).execute()

        if args.operation == self.DEACTIVATE:
            return DeactivateParser(extra).execute()
