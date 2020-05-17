from argparse import ArgumentParser
from cpmg.initializer import CPMGInitializer
from cpmg.exporters import RegioSQMExporter
import cpmg.importers as i
import cpmg.models as m


def import_type(model_type):
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

    return importer.import_data()


def initialize():
    return CPMGInitializer().initialize()


def export_data():
    return RegioSQMExporter().export_regiosqm_smiles_file()


class InitializeArgParser:

    def __init__(self):
        parser = ArgumentParser()
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--import', choices=[m.Backbone.STRING, m.Connection.STRING, m.Template.STRING,
                                                m.Sidechain.STRING, m.Monomer.STRING, m.RegioSQMPrediction.STRING,
                                                m.pKaPrediction.STRING], dest='load',
                           help='Choose which types of molecules to import')
        group.add_argument('--init', action='store_true', help='Automatically import all importable types')
        group.add_argument('--export', action='store_true',
                           help='Exports the RegioSQM SMILES file for generating predictions')

        self.args = parser.parse_args()

    def execute(self):
        if self.args.load is not None:
            return import_type(self.args.load)

        if self.args.init:
            return initialize()

        if self.args.export:
            return export_data()


if __name__ == "__main__":
    initializer = InitializeArgParser()
    initializer.execute()
