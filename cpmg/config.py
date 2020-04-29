import os
from collections import namedtuple

################################################ Directories ################################################
PROJECT_DIR = os.environ['PROJECT_DIR']
DATA_DIR = os.path.join(PROJECT_DIR, 'data')
# DATA_DIR = os.path.join('/u/scratch/e/ericdang', 'data')
LOG_DIR = os.path.join(PROJECT_DIR, 'logs')
TMP_DIR = os.path.join(PROJECT_DIR, 'tmp')
#TMP_DIR = os.environ['TMPDIR']
STDOUT_DIR = os.path.join(PROJECT_DIR, 'output')
TEST_DIR = os.path.join(PROJECT_DIR, 'tests')
IMPORT_DIR = os.path.join(DATA_DIR, 'imports')

################################################ Miscellaneous ################################################
CAPACITY = 1000000
DATA_FORMAT = 'hdf5'
NUM_PROCS = None
TASKS_PER_CHILD = None

##################################################### HDF5Params ######################################################
HDF5_FILEPATH = os.path.join(DATA_DIR, 'macrocycle_lib.hdf5')
COMPRESSION = 'gzip'
COMPRESSION_OPTS = 9

##################################################### RegioSQM ######################################################
REGIOSQM_SMILES_FILEPATH = os.path.join(DATA_DIR, 'external', 'regiosqm_smiles.smiles')

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
