
import os
from collections import namedtuple

from pymongo import ASCENDING

################################################ Directories ################################################
PROJECT_DIR = os.environ['PROJECT_DIR']
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
# DATA_DIR = os.path.join('/u/scratch/e/ericdang', 'data')
LOG_DIR = os.path.join(PROJECT_DIR, 'logs')
TMP_DIR = os.path.join(PROJECT_DIR, 'tmp')
#TMP_DIR = os.environ['TMPDIR']

################################################ Miscellaneous ################################################
CAPACITY = 1000000
DATA_FORMAT = 'json'
NUM_PROCS = None
TASKS_PER_CHILD = None

################################################ Filter Constants ################################################
PKA_CUTOFF = 13.5
MAX_MW = 1200
MAX_ROTATABLE_BONDS = 10
MAX_TPSA = 200

############################################ ConformerGenerator Parameters #############################################
REPEATS_PER_CUT = 1
NUM_CONFS_EMBED = 3
NUM_CONFS_GENETIC = 5
NUM_CONFS_ROTAMER_SEARCH = 15
FORCE_FIELD = 'MMFF94s'
DIELECTRIC = 1.0
SCORE = 'energy'
MIN_RMSD = 0.5
ENERGY_DIFF = 5
EMBED_PARAMS = None
SMALL_ANGLE_GRAN = 5
LARGE_ANGLE_GRAN = 15
CLASH_THRESHOLD = 0.9
DISTANCE_INTERVAL = [1.0, 2.5]
NUM_THREADS = 0
MAX_ITERS = 1000
EXTRA_ITERS = 50
MIN_MACRO_RING_SIZE = 10
ConformerArgs = namedtuple('ConformerArgs', 'repeats_per_cut num_confs_embed num_confs_genetic num_confs_rotamer_search force_field \
                           dielectric score min_rmsd energy_diff embed_params small_angle_gran large_angle_gran \
                           clash_threshold distance_interval num_threads max_iters extra_iters min_macro_ring_size')
CONFORMER_ARGS = ConformerArgs(REPEATS_PER_CUT, NUM_CONFS_EMBED, NUM_CONFS_GENETIC, NUM_CONFS_ROTAMER_SEARCH, FORCE_FIELD, DIELECTRIC,
                               SCORE, MIN_RMSD, ENERGY_DIFF, EMBED_PARAMS, SMALL_ANGLE_GRAN, LARGE_ANGLE_GRAN,
                               CLASH_THRESHOLD, DISTANCE_INTERVAL, NUM_THREADS, MAX_ITERS, EXTRA_ITERS, MIN_MACRO_RING_SIZE)

########################################## EbejerConformerGenerator Parameters #########################################
D_MIN = 0.5
NUM_CONFS = 50

################################################ MongoDB Schema ################################################
HOST = 'localhost'
PORT = 27017
DATABASE = 'test'

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
                         'template_pictet_spangler', 'unmasked_aldehyde_cyclization'],
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

# used for storing filters
COL3_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['type'],
        'properties': {
            'type': {
                'enum': ['regiosqm', 'pka'],
                'description': 'The filter that has been applied to the molecule'
            }
        }
    }
}

# used for recording metadata such as used _ids and indices
COL4_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['_id'],
        'properties': {
            'type': {
                'bsonType': 'string',
                'description': 'The type of molecule this counter applies to'
            }
        }
    }
}
VALIDATORS = [COL1_VALIDATOR, COL2_VALIDATOR, COL3_VALIDATOR, COL4_VALIDATOR]

COL1_INDEX = [[('kekule', ASCENDING)]]
COL2_INDEX = [[('smarts', ASCENDING)]]
COL3_INDEX = None
COL4_INDEX = None
INDICES = [COL1_INDEX, COL2_INDEX, COL3_INDEX, COL4_INDEX]
