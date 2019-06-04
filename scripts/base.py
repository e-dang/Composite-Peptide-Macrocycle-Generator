from pathlib import Path
from utils import Database
import json


class Base():

    def __init__(self, project_dir=Path(__file__).resolve().parents[1]):

        self.project_dir = project_dir

    def write_json(self, data, fp_out):

        try:
            with open(fp_out, 'w') as f:
                json.dump(data, f)
        except OSError:
            print('File does not exist:', fp_out)
            return False
        except TypeError:
            print('Object contains invalid types')
            return False

        return True

    def write_txt(self, data, fp_out):
        with open(fp_out, 'w') as f:
            for dp in data:
                f.write(dp + '\n')

        return True
