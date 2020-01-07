from abc import ABC, abstractmethod

import config
import project_io


class IDataInitializer(ABC):

    @abstractmethod
    def initialize(self):
        pass


class RecordInitializer(IDataInitializer):

    def __init__(self, data_format=config.DATA_FORMAT):

        if data_format == 'json':
            id_io = project_io.JsonIDIO()
            index_io = project_io.JsonIndexIO()
        elif data_format == 'mongo':
            id_io = project_io.MongoIDIO()
            index_io = project_io.MongoIndexIO()

        self.initializers = [IDInitializer(id_io), IndexInitializer(index_io)]

    def initialize(self):

        for initializer in self.initializers:
            initializer.initialize()


class IDInitializer(IDataInitializer):

    def __init__(self, id_io):

        self.saver = id_io

    def initialize(self):

        self.saver.save({'_id': 'id',
                         'count': 0,
                         'prefix': ''})


class IndexInitializer(IDataInitializer):

    def __init__(self, index_io):

        self.saver = index_io

    def initialize(self):
        self.saver.save({'_id': 'index',
                         'index': 1})
