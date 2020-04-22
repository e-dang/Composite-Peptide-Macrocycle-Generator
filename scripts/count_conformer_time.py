import macrocycles.project_io as project_io
import macrocycles.utils as utils
import numpy as np
import glob
import json
from bson import json_util

PEPTIDE_LENGTH = 3
conformer_fp = project_io.JsonConformerIO.FILEPATH
filepaths = glob.glob(utils.attach_file_num(conformer_fp, PEPTIDE_LENGTH, '*'))
time = []
for filepath in filepaths:
    with open(filepath, 'r') as file:
        for doc in json_util.loads(json_util.dumps(json.load(file))):
            time.append(doc['time'])

print(f'Average Time peptide length {PEPTIDE_LENGTH}: {np.average(time)}')
print(f'Num conformers {len(time)}')
