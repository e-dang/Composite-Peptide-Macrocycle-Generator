from cpmg.repository import create_repository_initializer
from cpmg.importers import create_importers


class CPMGInitializer:

    def __init__(self):
        self.repo_initializer = create_repository_initializer()
        self.importers = create_importers()

    def initialize(self):
        self.repo_initializer.initialize()
        for importer in self.importers:
            importer.import_data()
