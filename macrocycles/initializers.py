from abc import ABC, abstractmethod

import project_io


class IDataInitializer(ABC):

    @abstractmethod
    def initialize(self):
        pass


class RecordInitializer(IDataInitializer):

    def __init__(self):

        self.initializers = [IDInitializer(), IndexInitializer()]

    def initialize(self):

        for initializer in self.initializers:
            initializer.initialize()


class IDInitializer(IDataInitializer):

    def initialize(self):

        project_io.get_id_io().save({'_id': 'id',
                                     'count': 0,
                                     'prefix': ''})


class IndexInitializer(IDataInitializer):

    def initialize(self):

        project_io.get_index_io().save({'_id': 'index',
                                        'index': 1})
