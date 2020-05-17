
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
STDOUT_DIR = os.path.join(PROJECT_DIR, 'output')

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


##################################################### Filepaths ######################################################
PEPTIDE_PLAN_FP = os.path.join(DATA_DIR, 'generated', 'peptide_plan.txt')
CONFORMER_PLAN_FP = os.path.join(DATA_DIR, 'generated', 'conformer_plan.txt')
REGIOSQM_SMILES_FP = os.path.join(DATA_DIR, 'external', 'regiosqm_mols.smiles')
RAW_REGIOSQM_RESULT_FP = os.path.join(DATA_DIR, 'external', 'regiosqm_results_nm_3.csv')
RAW_PKAS_FP = os.path.join(DATA_DIR, 'external', 'pkas.txt')
RAW_SIDECHAINS_FP = os.path.join(DATA_DIR, 'chemdraw', 'sidechains.sdf')
RAW_MONOMERS_FP = os.path.join(DATA_DIR, 'chemdraw', 'monomers.sdf')
IDS_FP = os.path.join(DATA_DIR, 'generated', 'ids.json')
INDEX_FP = os.path.join(DATA_DIR, 'generated', 'index.json')
SIDECHAINS_FP = os.path.join(DATA_DIR, 'generated', 'sidechains.json')
MONOMERS_FP = os.path.join(DATA_DIR, 'generated', 'monomers.json')
PEPTIDES_FP = os.path.join(DATA_DIR, 'generated', 'peptides.json')
TP_HYBRIDS_FP = os.path.join(DATA_DIR, 'generated', 'template_peptides.json')
MACROCYCLES_FP = os.path.join(DATA_DIR, 'generated', 'macrocycles.json')
CONFORMERS_FP = os.path.join(DATA_DIR, 'generated', 'conformers.json')
EBEJER_FP = os.path.join(DATA_DIR, 'generated', 'ebejer.json')
REACTIONS_FP = os.path.join(DATA_DIR, 'generated', 'reactions.json')
REGIOSQM_FP = os.path.join(DATA_DIR, 'generated', 'regiosqm.json')
PKAS_FP = os.path.join(DATA_DIR, 'generated', 'pka.json')
MW_FP = os.path.join(DATA_DIR, 'generated', 'mw.json')
RB_FP = os.path.join(DATA_DIR, 'generated', 'rb.json')
TPSA_FP = os.path.join(DATA_DIR, 'generated', 'tpsa.json')
FILEPATHS = {
    'peptide_plan': PEPTIDE_PLAN_FP,
    'conformer_plan': CONFORMER_PLAN_FP,
    'regiosqm_smiles': REGIOSQM_SMILES_FP,
    'raw_regiosqm_results': RAW_REGIOSQM_RESULT_FP,
    'raw_pkas': RAW_PKAS_FP,
    'raw_sidechains': RAW_SIDECHAINS_FP,
    'raw_monomers': RAW_MONOMERS_FP,
    'ids': IDS_FP,
    'index': INDEX_FP,
    'sidechains': SIDECHAINS_FP,
    'monomers': MONOMERS_FP,
    'peptides': PEPTIDES_FP,
    'tp_hybrids': TP_HYBRIDS_FP,
    'macrocycles': MACROCYCLES_FP,
    'conformers': CONFORMERS_FP,
    'ebjer': EBEJER_FP,
    'reactions': REACTIONS_FP,
    'regiosqm': REGIOSQM_FP,
    'pkas': PKAS_FP,
    'mw': MW_FP,
    'rb': RB_FP,
    'tpsa': TPSA_FP
}

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
