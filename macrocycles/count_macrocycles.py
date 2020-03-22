import os
import macrocycles.config as config
import multiprocessing
import json
from bson import json_util


GENERATED_DIR = os.path.join(config.DATA_DIR, 'generated')


def count_macrocycles(filepath):
    with open(filepath, 'r') as file:
        count = len(list(json_util.loads(json_util.dumps(json.load(file)))))

    return filepath, count


def output(counts):
    with open(os.path.join(GENERATED_DIR, 'counts.json'), 'w') as file:
        json.dump(counts, file)


def run():

    counts = {}
    filepaths = map(lambda x: os.path.join(GENERATED_DIR, x), filter(
        lambda x: 'macrocycles' in x, os.listdir(GENERATED_DIR)))
    with multiprocessing.Pool() as pool:
        result = pool.map_async(count_macrocycles, filepaths)
        for filepath, count in result.get():
            counts[filepath] = count

    output(counts)


if __name__ == "__main__":
    run()
