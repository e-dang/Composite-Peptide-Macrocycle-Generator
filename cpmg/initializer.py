from cpmg.importers import create_importers


class CPMGInitializer:

    def __init__(self):
        self.importers = create_importers()

    def initialize(self):
        for importer in self.importers:
            importer.import_data()

        return True
