from cpmg.importers import create_importers, create_mol_importers, create_prediction_importers


class CPMGInitializer:

    def initialize_predictions_only(self):
        for importer in create_prediction_importers():
            importer.import_data()

        return True

    def initialize_mols_only(self):
        for importer in create_mol_importers():
            importer.import_data()

        return True

    def initialize(self):
        for importer in create_importers():
            importer.import_data()

        return True
