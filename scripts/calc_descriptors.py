import json
import os
import glob
from bson import json_util
from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing

import macrocycles.config as config


def get_filepaths(glob_filepath):
    for filepath in glob.glob(glob_filepath):
        yield filepath


# def open_files(filepaths):
#     for filepath in filepaths:
#         with open(filepath, 'r') as file:
#             for doc in json_util.loads(json_util.dumps(json.load(file))):
#                 yield Chem.Mol(doc['binary'])


def open_files(filepaths):
    for filepath in filepaths:
        yield Chem.MolFromPDBFile(filepath)


def calc_rb(mol):
    return AllChem.CalcNumRotatableBonds(mol)


def calc_mw(mol):
    return AllChem.CalcExactMolWt(mol)


def calc_tpsa(mol):
    return AllChem.CalcTPSA(mol, includeSandP=True)


def run(data, func):
    results = []
    with multiprocessing.Pool(processes=config.NUM_PROCS) as pool:
        for value in pool.imap_unordered(func, data):
            results.append(value)

    return results


def save_results(results, name):
    with open(os.path.join(config.DATA_DIR, 'generated', f'{name}.json'), 'w') as file:
        json.dump(results, file)


if __name__ == "__main__":
    filepaths = get_filepaths(os.path.join('./pdbs', '*.pdb'))
    # filepaths = get_filepaths(os.path.join(config.DATA_DIR, 'generated', 'conformers*'))
    data = open_files(filepaths)
    func, name = calc_mw, 'mw'
    results = run(data, func)
    # print(results)
    save_results(results, name)
