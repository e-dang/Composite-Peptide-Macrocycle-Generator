import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
import glob
import os
import json
from collections import Counter, defaultdict
import macrocycles.config as config
import numpy as np
from scipy.stats import iqr


def open_files(glob_filepath):
    for filepath in glob.glob(glob_filepath):
        with open(filepath, 'r') as file:
            print(filepath)
            yield json.load(file), filepath


def parse_data(data):
    for doc, filepath in data:
        rmsd = [0.0000] + doc['rmsd']
        ring_rmsd = [0.00000] + doc['ring_rmsd']
        yield doc['energies'], rmsd, ring_rmsd, filepath


def get_unique_vals(values):
    s = set(val for val in values)
    return s


def plot(tuple_data):

    for energy, rmsd, ring_rmsd, filepath in tuple_data:
        energy = np.round(np.array(energy), decimals=1)
        rmsd = np.round(np.array(rmsd), decimals=1)
        energy_rmsd = list(zip(energy, rmsd))
        counter = Counter(energy_rmsd)

        partitions = defaultdict(list)
        for tup, count in counter.items():
            partitions[tup[0]].append((tup, count))

        merged_counts = {}
        for partition in partitions.values():
            partition_dict = {}
            for (energy, rmsd), counts in partition:
                if len(partition_dict) != 0:
                    for (h_energy, h_rmsd), h_counts in list(partition_dict.items()):
                        if abs(rmsd - h_rmsd) <= config.MIN_RMSD:
                            if rmsd < h_rmsd:
                                partition_dict[(energy, rmsd)] = h_counts + counts
                                partition_dict.pop((h_energy, h_rmsd))
                            else:
                                partition_dict[(h_energy, h_rmsd)] = h_counts + counts
                            break
                    else:
                        partition_dict[(energy, rmsd)] = counts
                else:
                    partition_dict[(energy, rmsd)] = counts
            merged_counts.update(partition_dict)

        for tup, count in list(merged_counts.items()):
            if count == 1:
                merged_counts.pop(tup)
            else:
                print(tup, count)
        print('\n')
        plt.figure(figsize=(200, 6))
        labels = [x.strip('(').strip(')') for x in list(map(str, merged_counts.keys()))]
        plt.bar(list(range(1, len(merged_counts) + 1)), merged_counts.values(),
                tick_label=labels)
        plt.xticks(rotation='vertical', fontsize=5)
        plt.title('Number of Times The Same Conformer is Found by ConfBuster++')
        plt.xlabel('Conformer with (Energy, RMSD) Tuple ($kcal \cdot mol^{-1}, Å$)')
        plt.ylabel('Count')
        plt.xlim(-0.5, len(merged_counts) - .5)
        _, basename = os.path.split(filepath)
        basename, ext = os.path.splitext(basename)
        print(basename)
        plt.show()
        # plt.savefig(f'repeats_{basename}.png')
        plt.cla()
        plt.clf()


if __name__ == "__main__":
    data = open_files(os.path.join(config.PROJECT_DIR, 'scripts', 'data', 'energy', '*.json'))
    tuple_data = parse_data(data)
    plot(tuple_data)
