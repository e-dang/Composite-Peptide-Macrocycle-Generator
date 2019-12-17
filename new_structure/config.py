
import os
from collections import namedtuple
from pymongo import ASCENDING

################################################ Directories ################################################
PROJECT_DIR = os.environ['PROJECT_DIR']
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
LOG_DIR = os.path.join(PROJECT_DIR, 'logs')
TMP_DIR = os.path.join(PROJECT_DIR, 'tmp')

################################################ Globals ################################################
CAPACITY = 500000000

################################################ MongoDB Schema ################################################
HOST = 'localhost'
PORT = 27017
DATABASE = 'test'
# MongoSettings = namedtuple('MongoSettings', ['host', 'port', 'database'])
# MONGO_SETTINGS = MongoSettings(MONGO_HOST, MONGO_PORT, MONGO_DATABASE)

COL1 = 'molecules'
COL2 = 'reactions'
COL3 = 'filters'
COL4 = 'records'
COLLECTIONS = [COL1, COL2, COL3, COL4]

# used for storing molecules
COL1_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['_id', 'type', 'binary', 'kekule'],
        'properties': {
            '_id': {
                'bsonType': 'string',
                'description': 'The id of the molecule'
            },
            'type': {
                'enum': ['sidechain', 'monomer', 'peptide', 'template_peptide', 'macrocycle'],
                'description': 'The type of molecule or data contained in the document'
            },
            'binary': {
                'bsonType': ['binData', 'array'],
                'description': 'The binary representation of the molecule'
            },
            'kekule': {
                'bsonType': ['string', 'array'],
                'description': 'The kekule SMILES string of the molecule'
            }
        }
    }
}

# used for storing reactions
COL2_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['_id', 'type', 'binary', 'smarts'],
        'properties': {
            '_id': {
                'bsonType': 'string',
                'description': 'The id of the molecule'
            },
            'type': {
                'enum': ['friedel_crafts', 'tsuji_trost', 'pictet_spangler', 'pyrrolo_indolene',
                         'template_pictet_spangler'],
                'description': 'The type of molecule or data contained in the document'
            },
            'binary': {
                'bsonType': ['binData', 'array'],
                'description': 'The binary representation of the molecule'
            },
            'smarts': {
                'bsonType': 'string',
                'description': 'The atom mapped SMARTS/SMILES string'
            }
        }
    }
}

# used for storing filtering results
COL3_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['_id', 'type'],
        'properties': {
            '_id': {
                'bsonType': 'string',
                'description': 'The id of the molecule'
            },
            'type': {
                'enum': ['regiosqm', 'pka'],
                'description': 'The filter that has been applied to the molecule'
            }
        }
    }
}

# used for recording metadata such as used _ids and last insertions
COL4_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['_id'],
        'properties': {
            'type': {
                'bsonType': 'string',
                'description': 'The type of molecule this counter applies to'
            }
            # 'count': {
            #     'bsonType': 'int',
            #     'maximum': 25,
            #     'minimum': 0,
            #     'description': 'Determines the letter to be appended to prefix'
            # },
            # 'prefix': {
            #     'bsonType': 'string',
            #     'description': 'The prefix to the ID determined by count'
            # }
        }
    }
}
VALIDATORS = [COL1_VALIDATOR, COL2_VALIDATOR, COL3_VALIDATOR, COL4_VALIDATOR]

COL1_INDEX = [[('kekule', ASCENDING)]]
COL2_INDEX = [[('smarts', ASCENDING)]]
# COL3_INDEX = [('_id', ASCENDING)]
# COL4_INDEX = [('_id', ASCENDING)]
COL3_INDEX = None
COL4_INDEX = None
INDICES = [COL1_INDEX, COL2_INDEX, COL3_INDEX, COL4_INDEX]
