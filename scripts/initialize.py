from argparse import ArgumentParser
from cpmg.initializer import CPMGInitializer
from cpmg.exporters import RegioSQMExporter
import cpmg.importers as i
import cpmg.models as m

IMPORT_OPERATION = 'import'
INIT_OPERATION = 'init'
EXPORT_OPERATION = 'export'


def execute(command_line_args):
    command_line_args.pop('func')
    operation = command_line_args.pop('operation')
    if operation == IMPORT_OPERATION:
        return execute_import(command_line_args)

    if operation == INIT_OPERATION:
        return execute_init()

    if operation == EXPORT_OPERATION:
        return execute_export()


def execute_import(command_line_args):
    model_type = command_line_args.pop('type')
    if model_type == m.Backbone.STRING:
        importer = i.BackboneImporter(i.JsonImporter())
    elif model_type == m.Connection.STRING:
        importer = i.ConnectionImporter(i.JsonImporter())
    elif model_type == m.Template.STRING:
        importer = i.TemplateImporter(i.JsonImporter())
    elif model_type == m.Sidechain.STRING:
        importer = i.SidechainImporter(i.JsonImporter())
    elif model_type == m.Monomer.STRING:
        importer = i.MonomerImporter(i.JsonImporter())
    elif model_type == m.RegioSQMPrediction.STRING:
        importer = i.RegioSQMPredictionImporter()
    elif model_type == m.pKaPrediction.STRING:
        importer = i.pKaPredictionImporter()
    elif model_type == m.PeptidePlan.STRING:
        importer = i.PeptidePlanImporter()

    return importer.import_data(**command_line_args)


def execute_init():
    return CPMGInitializer().initialize()


def execute_export():
    return RegioSQMExporter().export_regiosqm_smiles_file()


class InitializeArgParser:

    def __init__(self):
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest='operation')

        import_parser = subparsers.add_parser(IMPORT_OPERATION)
        import_parser.add_argument('type', choices=[m.Backbone.STRING, m.Connection.STRING, m.Template.STRING,
                                                    m.Sidechain.STRING, m.Monomer.STRING, m.RegioSQMPrediction.STRING,
                                                    m.pKaPrediction.STRING, m.PeptidePlan.STRING],
                                   help='The type of record to import.')
        import_parser.add_argument('-f', '--file', '--filepath', type=str, dest='filepath',
                                   help='The path to the file containing the data to import. Only valid for importing peptide plans.')
        import_parser.add_argument('-l', '--length', '--peptide_length', type=int, choices=[3, 4, 5], dest='peptide_length',
                                   help='The peptide length of the imported peptide plan.')
        import_parser.set_defaults(func=execute)

        init_parser = subparsers.add_parser(INIT_OPERATION)
        init_parser.set_defaults(func=execute)

        export_parser = subparsers.add_parser(EXPORT_OPERATION)
        export_parser.set_defaults(func=execute)

        args = parser.parse_args()
        self.return_val = args.func(vars(args))


if __name__ == "__main__":
    initializer = InitializeArgParser()
