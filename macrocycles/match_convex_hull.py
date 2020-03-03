import glob
import json
import shutil
import os

import macrocycles.config as config
import macrocycles.exceptions as exceptions


def load_pdb_filepaths(glob_path):
    return list(glob.glob(glob_path))


def load_idxs(filepath):
    idxs = []
    with open(filepath, 'r') as file:
        for line in file.readlines():
            idxs.append(int(line.split(',')[0]))

    return idxs


def match_and_move(idxs, filepaths):
    for idx in idxs:
        path = filepaths[idx]
        parent, file = os.path.split(path)
        parent, _ = os.path.split(parent)
        shutil.copyfile(filepaths[idx], os.path.join(parent, 'convex_hull', file))


filepaths = load_pdb_filepaths('/u/home/e/ericdang/project-pharran/conformers/*.pdb')
idxs = load_idxs(os.path.join(config.DATA_DIR, 'external', 'PCA_convex_hull.txt'))
match_and_move(idxs, filepaths)
