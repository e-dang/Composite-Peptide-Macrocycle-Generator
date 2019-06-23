"""
Written by Eric Dang.
github: https://github.com/e-dang
email: edang830@gmail.com
"""

import os
from collections import namedtuple
from pymongo import ASCENDING


####################################### project constants #######################################
PROJECT_DIR = os.environ['PROJECT_DIR']
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
LOG_DIR = os.path.join(PROJECT_DIR, 'logs')

####################################### MongoDataBase constants #######################################

# default connection settings
MONGO_HOST = 'localhost'
MONGO_PORT = 27017
MONGO_DATABASE = 'test'
MongoSettings = namedtuple('MongoSettings', ['host', 'port', 'database'])
MONGO_SETTINGS = MongoSettings(MONGO_HOST, MONGO_PORT, MONGO_DATABASE)

# default collection settings
COL1 = 'molecules'
COL2 = 'reactions'
COL3 = 'filters'
COL4 = 'counter'
Collections = namedtuple('Collections', ['mols', 'rxns', 'filters', 'counter'])
COLLECTIONS = Collections(COL1, COL2, COL3, COL4)

# used for storing molecules
COL1_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['ID', 'type', 'smiles', 'kekule'],
        'properties': {
            'ID': {
                'bsonType': 'string',
                'description': 'The ID of the molecule'
            },
            'type': {
                'enum': ['parent_side_chain', 'side_chain', 'backbone', 'monomer', 'peptide', 'template', 'tp_hybrid', 'candidate', 'conformer'],
                'description': 'The type of molecule or data contained in the document'
            },
            'smiles': {
                'bsonType': 'string',
                'description': 'The SMILES string of the molecule'
            },
            'kekule': {
                'bsonType': 'string',
                'description': 'The kekule SMILES string of the molecule'
            }
        }
    }
}

# used for storing reactions
COL2_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['ID', 'type'],
        'properties': {
            'ID': {
                'bsonType': 'string',
                'description': 'The ID of the molecule'
            },
            'type': {
                'enum': ['side_chain', 'template', 'reaction'],
                'description': 'The type of molecule or data contained in the document'
            }
        }
    }
}

# used for storing filtering results
COL3_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['ID', 'type'],
        'properties': {
            'ID': {
                'bsonType': 'string',
                'description': 'The ID of the molecule'
            },
            'type': {
                'enum': ['regiosqm', 'indole', 'manual'],
                'description': 'The filter that has been applied to the molecule'
            }
        }
    }
}

# used for calculating new IDs
COL4_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['type', 'count', 'prefix'],
        'properties': {
            'type': {
                'bsonType': 'string',
                'description': 'The type of molecule this counter applies to'
            },
            'count': {
                'bsonType': 'int',
                'maximum': 25,
                'minimum': 0,
                'description': 'Determines the letter to be appended to prefix'
            },
            'prefix': {
                'bsonType': 'string',
                'description': 'The prefix to the ID determined by count'
            }
        }
    }
}
Validators = namedtuple('Validators', ['mols', 'rxns', 'filters', 'counter'])
VALIDATORS = Validators(COL1_VALIDATOR, COL2_VALIDATOR, COL3_VALIDATOR, COL4_VALIDATOR)

COL1_INDEX = [[('ID', ASCENDING)], [('smiles', ASCENDING)]]
COL2_INDEX = [('ID', ASCENDING)]
COL3_INDEX = [('ID', ASCENDING)]
COL4_INDEX = [('ID', ASCENDING)]
Index = namedtuple('Index', ['mols', 'rxns', 'filters', 'counter'])
INDEX = Index(COL1_INDEX, COL2_INDEX, COL3_INDEX, COL4_INDEX)

####################################### DataInitializer #######################################
DI_INPUT_DIR = 'chemdraw'

####################################### SideChainModifier #######################################
HETERO_MAP_NUM = 1
CONN_MAP_NUM = 2
SCM_INPUT_DIR = 'pre_monomer'
SCM_OUTPUT_DIR = 'pre_monomer'
SCM_INPUT_COL = [COL1]
SCM_DOC_TYPE = ['parent_side_chain']
SCM_OUTPUT_COL = COL1
Connections = namedtuple('Connections', ['smarts', 'mod_array'])
CONNECTIONS = [Connections(con, mod) for con, mod in [(f'[CH3:{CONN_MAP_NUM}]', [0, 3]),
                                                      (f'[CH3][CH2:{CONN_MAP_NUM}]', [1, 3]),
                                                      (f'[CH3][CH2][CH2:{CONN_MAP_NUM}]', [2, 3])]]

####################################### MonomerGenerator #######################################
BB_MAP_NUM = 1
SC_MAP_NUM = 2
STEREOCHEMISTRY = ['CW', 'CCW']
MG_INPUT_DIR = 'pre_monomer'
MG_OUTPUT_DIR = 'monomers'
MG_INPUT_COL = [COL1, COL1]
MG_DOC_TYPE = ['side_chain', 'backbone']
MG_OUTPUT_COL = COL1
