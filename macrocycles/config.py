
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

################################################ Map Numbers ################################################
# SideChainGenerator
CONN_MAP_NUM = 1
PSC_MAP_NUM = 2

# MonomerGenerator
BB_MAP_NUM = 1
SC_MAP_NUM = 2

# PeptideGenerator
MONO_NITROGEN_MAP_NUM = 1
PEP_CARBON_MAP_NUM = 2

# TPHybridGenerator
TEMP_CARBON_MAP_NUM = 1
PEP_NITROGEN_MAP_NUM = 2

# ReactionGenerator
TEMP_WILDCARD_MAP_NUM = TEMP_CARBON_MAP_NUM
TEMP_EAS_MAP_NUM = 2
SC_WILDCARD_MAP_NUM = 3
SC_EAS_MAP_NUM = 4

################################################ MongoDB Schema ################################################
MONGO_HOST = 'localhost'
MONGO_PORT = 27017
MONGO_DATABASE = 'tests'
MongoSettings = namedtuple('MongoSettings', ['host', 'port', 'database'])
MONGO_SETTINGS = MongoSettings(MONGO_HOST, MONGO_PORT, MONGO_DATABASE)

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
                'enum': ['parent_side_chain', 'connection', 'side_chain', 'backbone', 'monomer', 'peptide', 'template', 'tp_hybrid', 'macrocycle', 'conformer'],
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
                'enum': ['template', 'friedel_crafts', 'tsuji_trost', 'pictet_spangler'],
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
                'enum': ['regiosqm', 'indole', 'manual'],
                'description': 'The filter that has been applied to the molecule'
            }
        }
    }
}

# used for recording metadata such as used _ids and last insertions
COL4_VALIDATOR = {
    '$jsonSchema': {
        'bsonType': 'object',
        'required': ['type'],
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
####################################### DataInitializer #######################################
DI_DEFAULTS = {
    'inputs': {
        'fp_psc': os.path.join(DATA_DIR, 'chemdraw', 'pre_monomer', 'side_chains_likely1.sdf'),
        'fp_connections': os.path.join(DATA_DIR, 'chemdraw', 'pre_monomer', 'connections.sdf'),
        'fp_backbones': os.path.join(DATA_DIR, 'chemdraw', 'pre_monomer', 'backbones.sdf'),
        'fp_monomers': os.path.join(DATA_DIR, 'chemdraw', 'monomers', 'modified_prolines.sdf'),
        'fp_templates': os.path.join(DATA_DIR, 'chemdraw', 'templates.sdf')
    },
    'outputs': {
        'col_psc': COL1,
        'fp_psc': os.path.join(DATA_DIR, 'starter', 'parent_side_chains.json'),
        'col_connections': COL1,
        'fp_connections': os.path.join(DATA_DIR, 'starter', 'connections.json'),
        'col_backbones': COL1,
        'fp_backbones': os.path.join(DATA_DIR, 'starter', 'backbones.json'),
        'col_monomers': COL1,
        'fp_monomers': os.path.join(DATA_DIR, 'starter', 'monomers.json'),
        'col_templates': COL1,
        'fp_templates': os.path.join(DATA_DIR, 'starter', 'templates.json'),
        'col_rtemplates': COL2,
        'fp_rtemplates': os.path.join(DATA_DIR, 'starter', 'rtemplates.json')
    }
}

####################################### SideChainGenerator #######################################
SCM_DEFAULTS = {
    'inputs': {
        'col_psc': COL1,
        'col_connections': COL1,
    },
    'outputs': {
        'col_side_chains': COL1
    }
}

####################################### MonomerGenerator #######################################
MG_DEFAULTS = {
    'inputs': {
        'col_backbones': COL1,
        'col_side_chains': COL1
    },
    'outputs': {
        'col_monomers': COL1
    }
}

####################################### PeptideGenerator #######################################
PG_DEFAULTS = {
    'inputs': {
        'col_monomers': COL1
    },
    'outputs': {
        'col_peptides': COL1
    }
}

####################################### TPHybridGenerator #######################################
TPHG_DEFAULTS = {
    'inputs': {
        'col_peptides': COL1,
        'col_templates': COL1
    },
    'outputs': {
        'col_tp_hybrids': COL1
    }
}

####################################### MacrocycleGenerator #######################################
MCG_DEFAULTS = {
    'inputs': {
        'col_tp_hyrbids': COL1,
        'col_reactions': COL2
    },
    'outputs': {
        'col_macrocycles': COL1,
        'json_macrocycles': os.path.join(DATA_DIR, 'macrocycles.json')
    }
}

####################################### ConformerGenerator #######################################
CG_DEFAULTS = {
    'inputs': {
        'col_macrocycles': COL1,
        'json_macrocycles': MCG_DEFAULTS['outputs']['json_macrocycles']
    },
    'outputs': {
        'col_conformers': COL1,
        'json_conformers': os.path.join(DATA_DIR, 'conformers', 'conformers.json'),
        'fp_conformers': os.path.join(DATA_DIR, 'conformers'),
        'tmp_molecule': os.path.join(TMP_DIR, 'conf_macrocycle.sdf'),
        'tmp_genetic_results': os.path.join(TMP_DIR, 'genetic_results.sdf')
    }
}

####################################### ReactionGenerator #######################################
RG_DEFAULTS = {
    'inputs': {
        'col_templates': COL2,
        'col_side_chains': COL1,
    },
    'outputs': {
        'col_reactions': COL2
    },
    'template_dependent_rxns': ['PictetSpangler']
}

####################################### RegioSQMFilter #######################################
RSQM_DEFAULTS = {
    'inputs': {
        'col_filters': COL3,
        'col_macrocycles': COL1
    },
    'outputs': {
        'col_filtered': COL1
    }
}

####################################### DEFAULT #######################################
DEFAULTS = {
    'DataInitializer': DI_DEFAULTS,
    'SideChainGenerator': SCM_DEFAULTS,
    'MonomerGenerator': MG_DEFAULTS,
    'PeptideGenerator': PG_DEFAULTS,
    'TPHybridGenerator': TPHG_DEFAULTS,
    'MacrocycleGenerator': MCG_DEFAULTS,
    'ConformerGenerator': CG_DEFAULTS,
    'ReactionGenerator': RG_DEFAULTS,
    'RegioSQMFilter': RSQM_DEFAULTS
}
