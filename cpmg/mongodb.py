from cpmg.parallelism import Parallelism  # noqa
import mpi4py  # noqa
mpi4py.rc(initialize=False, finalize=False)  # noqa
from mpi4py import MPI  # noqa

import atexit
from subprocess import Popen

from pymongo import MongoClient
from pymongo.errors import ConnectionFailure

import cpmg.utils as utils
import cpmg.ranges as ranges
import cpmg.config as config


class DataBaseDaemon:
    __DAEMON = None

    @classmethod
    def start(cls):
        if cls.__DAEMON is None:
            if not Parallelism.is_distributed() or (Parallelism.is_distributed() and MPI.COMM_WORLD.Get_rank() == 0):
                cls.__DAEMON = Popen([config.MONGO_DB_EXECUTABLE, '--dbpath', config.MONGO_DB_DATA_PATH, '--logappend',
                                      '--logpath', config.MONGO_DB_LOG_PATH, '--port', config.MONGO_DB_DAEMON_PORT])
                atexit.register(cls.close)

    @classmethod
    def close(cls):
        cls.__DAEMON.terminate()


class MongoDataBase:

    def __init__(self):
        self.__client = None
        self.__database = None
        self.__collection = None

    def __repr__(self):
        return self.COLLECTION + ' - ' + str(self.get_num_records())  # pylint: disable=no-member

    def load(self, key):
        self.__make_connection()
        if key.peptide_length is None:
            return self.__load_no_pep_length(key)
        else:
            return self.__load_specific_pep_length(key)

    def save(self, data):
        self.__make_connection()
        return self.__collection.insert_many(utils.to_list(data), ordered=False).inserted_ids

    def get_num_records(self):
        self.__make_connection()
        return self.__collection.count_documents({})

    def remove_records(self, key):
        self.__make_connection()
        return self.__collection.delete({'_id': {'$in': utils.to_list(key)}}).deleted_count

    def mark_complete(self, key):
        self.__make_connection()
        return self.__collection.update_many({'_id': {'$in': utils.to_list(key)}}, {'$set': {'completed': True}})

    def deactivate_records(self, key):
        return self.mark_complete(key)

    def __make_connection(self):
        if None in (self.__client, self.__database, self.__collection):
            DataBaseDaemon.start()
            self.__connect_to_client()
            self.__database = self.__client[config.MONGO_DB_DATABASE]
            self.__collection = self.__database[self.COLLECTION]  # pylint: disable=no-member

    def __connect_to_client(self):
        from time import sleep
        for _ in range(10):  # try 10 times then fail if still not working
            try:
                self.__client = MongoClient(config.MONGO_DB_HOST, config.MONGO_DB_CLIENT_PORT)
                break
            except ConnectionFailure:
                sleep(1)

    def __load_no_pep_length(self, key):
        if isinstance(key.key, ranges.IndexKey):
            return self.__collection.find({'_id': {'$in': key}})
        else:
            return self.__collection.find({'completed': False})

    def __load_specific_pep_length(self, key):
        if isinstance(key.key, ranges.IndexKey):
            return self.__collection.find({'_id': {'$in': key}, 'length': key.peptide_length})
        else:
            return self.__collection.find({'completed': False, 'length': key.peptide_length})


class ConnectionMongoRepository(MongoDataBase):
    COLLECTION = 'connection'


class BackboneMongoRepository(MongoDataBase):
    COLLECTION = 'backbone'


class TemplateMongoRepository(MongoDataBase):
    COLLECTION = 'template'


class SidechainMongoRepository(MongoDataBase):
    COLLECTION = 'sidechain'


class MonomerMongoRepository(MongoDataBase):
    COLLECTION = 'monomer'


class PeptideMongoRepository(MongoDataBase):
    COLLECTION = 'peptide'


class TemplatePeptideMongoRepository(MongoDataBase):
    COLLECTION = 'template_peptide'


class MacrocycleMongoRepository(MongoDataBase):
    COLLECTION = 'macrocycle'


class ConformerMongoRepository(MongoDataBase):
    COLLECTION = 'conformer'


class ReactionMongoRepository(MongoDataBase):
    COLLECTION = 'reaction'


class RegioSQMMongoRepository(MongoDataBase):
    COLLECTION = 'regiosqm'


class pKaMongoRepository(MongoDataBase):
    COLLECTION = 'pka'


class PeptidePlanMongoRepository(MongoDataBase):
    COLLECTION = 'peptide_plan'


class MongoRepository:
    def __init__(self):
        self.connection_repo = ConnectionMongoRepository()
        self.backbone_repo = BackboneMongoRepository()
        self.template_repo = TemplateMongoRepository()
        self.sidechain_repo = SidechainMongoRepository()
        self.monomer_repo = MonomerMongoRepository()
        self.peptide_repo = PeptideMongoRepository()
        self.template_peptide_repo = TemplatePeptideMongoRepository()
        self.macrocycle_repo = MacrocycleMongoRepository()
        self.conformer_repo = ConformerMongoRepository()
        self.reaction_repo = ReactionMongoRepository()
        self.regiosqm_repo = RegioSQMMongoRepository()
        self.pka_repo = pKaMongoRepository()
        self.peptide_plan_repo = PeptidePlanMongoRepository()

    def __repr__(self):

        string = '/'
        for instance in self.__dict__.values():
            string += instance.__repr__()

        return string

    @classmethod
    def instance(cls):
        if not hasattr(cls, '_instance'):
            cls._instance = cls()
        return cls._instance
