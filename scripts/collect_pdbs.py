import os
import shutil
import json
from bson import json_util


def read_text_file(filepath):
    with open(filepath, 'r') as file:
        for line in file.readlines():
            yield line.strip('\n')


def read_json_file(filepath):
    with open(filepath, 'r') as file:
        for doc in json_util.loads(json_util.dumps(json.load(file))):
            yield doc


def copy_pdbs_to_current_dir(desired_pdbs_filepath, pdbs_filepath, destination_dir):
    selected_pdbs = list(read_text_file(desired_pdbs_filepath))
    for filename in os.listdir(pdbs_filepath):
        if filename in selected_pdbs:
            shutil.copyfile(os.path.join(pdbs_filepath, filename), os.path.join(destination_dir, filename))


if __name__ == "__main__":
    copy_pdbs_to_current_dir('./pdbs.txt', '/u/home/e/ericdang/project-pharran/conformers',
                             os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pdbs'))
