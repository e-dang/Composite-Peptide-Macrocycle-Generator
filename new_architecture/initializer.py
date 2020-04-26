from new_architecture.repository.repository import create_repository_initializer
from new_architecture.importers import create_importers


class CPMGInitializer:

    def __init__(self):
        self.repo_initializer = create_repository_initializer()
        self.importers = create_importers()

    def initialize(self):
        self.repo_initializer.initialize()
        for importer in self.importers:
            print(importer.loader.search_dir)
            importer.import_data()
