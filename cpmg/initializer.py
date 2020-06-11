from cpmg.importers import create_importers, create_mol_importers, create_prediction_importers


class CPMGInitializer:
    MOL_TYPE = 'mols'
    PREDICTION_TYPE = 'predictions'
    ALL = 'all'

    def initialize(self, record_types):
        if record_types == self.MOL_TYPE:
            return self._initialize_mols_only()

        if record_types == self.PREDICTION_TYPE:
            return self._initialize_predictions_only()

        if record_types == self.ALL:
            return self._initialize_all()

        raise ValueError('Unrecognized record type for initialization!')

    def _initialize_predictions_only(self):
        for importer in create_prediction_importers():
            importer.import_data()

        return True

    def _initialize_mols_only(self):
        for importer in create_mol_importers():
            importer.import_data()

        return True

    def _initialize_all(self):
        for importer in create_importers():
            importer.import_data()

        return True
