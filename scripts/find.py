import json
import os
import glob


def find(filepath, ids):
    filepaths = glob.glob(filepath)
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            for doc in json.load(file):
                if doc['_id'] in ids:
                    print(doc['kekule'])


if __name__ == "__main__":
    find('/u/home/e/ericdang/project-pharran/conformers*', ('FJGB2KW0ER0dv80m28s2'))
