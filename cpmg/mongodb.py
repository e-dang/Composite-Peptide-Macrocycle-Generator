# THESE MPI IMPORTS MUST COME FIRST
from cpmg.parallelism import Parallelism  # noqa
import mpi4py  # noqa
mpi4py.rc(initialize=False, finalize=False)  # noqa
from mpi4py import MPI  # noqa

import atexit
import json
import os
import socket
from subprocess import Popen

from bson import json_util
from pymongo import DESCENDING, ASCENDING, MongoClient
from pymongo.errors import BulkWriteError, ConnectionFailure

import cpmg.config as config
import cpmg.ranges as ranges
import cpmg.utils as utils

IP_ADDR = socket.gethostbyname(socket.gethostname())

MOL_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['kekule'],
        'properties': {
            'kekule': {
                'bsonType': 'string',
                'description': 'The atom kekule SMILES string of the molecule'
            }
        }
    }
}

PREDICTION_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['predictions', 'reacting_mol'],
        'properties': {
            'predictions': {
                'bsonType': ['array', 'object'],
                'description': 'The predictions'
            },
            'reacting_mol': {
                'bsonType': 'string',
                'description': 'The SMILES string of the reacting mol it is applicable to'
            }
        }
    }
}


class DataBaseDaemon:
    __DAEMON = None

    @classmethod
    def start(cls):
        if cls.__DAEMON is None:
            if not Parallelism.is_distributed() or (Parallelism.is_distributed() and MPI.COMM_WORLD.Get_rank() == 0):
                cls.__DAEMON = Popen([config.MONGO_DB_EXECUTABLE, '--dbpath', config.MONGO_DB_DATA_PATH, '--logappend',
                                      '--logpath', config.MONGO_DB_LOG_PATH, '--port', config.MONGO_DB_DAEMON_PORT,
                                      '--bind_ip', IP_ADDR])
                atexit.register(cls.close)

    @classmethod
    def close(cls):
        cls.__DAEMON.terminate()


# pylint: disable=no-member
class MongoDataBase:

    def __init__(self):
        self.__client = None
        self.__database = None
        self.__collection = None
        self.failed_instances = []

    def __repr__(self):
        return self.COLLECTION + ' - ' + str(self.get_num_records())

    def load(self, key):
        self.__make_connection()
        if key.peptide_length is None:
            return self.__load_no_pep_length(key)
        else:
            return self.__load_specific_pep_length(key)

    def save(self, data):
        self.__make_connection()
        try:
            return self.__collection.insert_many(utils.to_list(data), ordered=False).inserted_ids
        except BulkWriteError as err:
            self.failed_instances.extend([failed_instance['op'] for failed_instance in err.details['writeErrors']])
            print(
                f'Failed to write {len(self.failed_instances)} {self.COLLECTION}s! These are probably duplicates, but dumping to data directory anyway...', flush=True)
            utils.save_json(json.loads(json_util.dumps(self.failed_instances)),
                            utils.rotate_file(os.path.join(config.DATA_DIR, f'{self.COLLECTION}.json')))

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

    def activate_records(self, key):
        self.__make_connection()
        return self.__collection.update_many({'_id': {'$in': utils.to_list(key)}}, {'$set': {'completed': False}})

    def __make_connection(self):
        if None in (self.__client, self.__database, self.__collection):
            DataBaseDaemon.start()
            self.__connect_to_client()
            self.__database = self.__client[config.MONGO_DB_DATABASE]
            self.__setup_validator()
            self.__collection = self.__database[self.COLLECTION]

    def __connect_to_client(self):
        from time import sleep
        for _ in range(10):  # try 10 times then fail if still not working
            try:
                self.__client = MongoClient(IP_ADDR, config.MONGO_DB_CLIENT_PORT)
                break
            except ConnectionFailure:
                sleep(1)

    def __setup_validator(self):
        if self.COLLECTION not in self.__database.collection_names():
            collection = self.__database.create_collection(self.COLLECTION)
            self.__database.command({'collMod': self.COLLECTION,
                                     'validator': self.VALIDATOR,
                                     'validationLevel': self.VALIDATION_LEVEL})
            try:
                collection.create_index(self.INDEX, unique=True)
            except AttributeError:
                pass

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
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class BackboneMongoRepository(MongoDataBase):
    COLLECTION = 'backbone'
    VALIDATION_LEVEL = 'strict'
    INDEX = [('mapped_kekule', ASCENDING)]
    VALIDATOR = {
        '$jsonSchema': {
            'bsonType': 'object',
            'required': ['mapped_kekule'],
            'properties': {
                'kekule': {
                    'bsonType': 'string',
                    'description': 'The atom kekule SMILES string of the molecule'
                }
            }
        }
    }


class TemplateMongoRepository(MongoDataBase):
    COLLECTION = 'template'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class SidechainMongoRepository(MongoDataBase):
    COLLECTION = 'sidechain'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class MonomerMongoRepository(MongoDataBase):
    COLLECTION = 'monomer'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class PeptideMongoRepository(MongoDataBase):
    COLLECTION = 'peptide'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class TemplatePeptideMongoRepository(MongoDataBase):
    COLLECTION = 'template_peptide'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class MacrocycleMongoRepository(MongoDataBase):
    COLLECTION = 'macrocycle'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class ConformerMongoRepository(MongoDataBase):
    COLLECTION = 'conformer'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = MOL_VALIDATOR
    INDEX = [('kekule', ASCENDING)]


class ReactionMongoRepository(MongoDataBase):
    COLLECTION = 'reaction'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = REACTION_VALIDATOR = {
        '$jsonSchema': {
            'bsonType': 'object',
            'required': ['type', 'smarts'],
            'properties': {
                'type': {
                    'bsonType': 'string',
                    'description': 'The type of reaction it is'
                },
                'smarts': {
                    'bsonType': 'string',
                    'description': 'The SMARTS string that represents the reaction'
                }
            }
        }
    }


class RegioSQMMongoRepository(MongoDataBase):
    COLLECTION = 'regiosqm'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = PREDICTION_VALIDATOR


class pKaMongoRepository(MongoDataBase):
    COLLECTION = 'pka'
    VALIDATION_LEVEL = 'strict'
    VALIDATOR = PREDICTION_VALIDATOR


class PeptidePlanMongoRepository(MongoDataBase):
    COLLECTION = 'peptide_plan'
    VALIDATION_LEVEL = 'strict'
    INDEX = [('combination', DESCENDING), ('length', ASCENDING)]
    VALIDATOR = {
        '$jsonSchema': {
            'bsonType': 'object',
            'required': ['combination', 'length'],
            'properties': {
                'combination': {
                    'bsonType': 'string',
                    'description': 'The indices of the monomers to form into a peptide'
                },
                'length': {
                    'bsonType': 'int',
                    'description': 'The length of the resulting peptide. Note this is not always equal to the length of the combinations field'
                }
            }
        }
    }


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
