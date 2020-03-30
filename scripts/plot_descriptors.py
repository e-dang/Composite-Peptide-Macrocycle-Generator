import json
import matplotlib.pyplot as plt
import os
import macrocycles.config as config
from scipy.stats import iqr
import numpy as np


def open_file(filepath):
    with open(filepath, 'r') as file:
        for val in json.load(file):
            yield val


def plot_data(data, title, x_label, save=False):
    data = list(data)
    maximum = max(data)
    minimum = min(data)
    bin_width = 2 * iqr(data) * (len(data) ** (-1 / 3))
    bin_width = bin_width if bin_width > 1 else 1
    bins = int(round((maximum - minimum) / bin_width))
    print(bins)
    _, _, patches = plt.hist(data, bins=bins, color='royalblue', ec='black')
    # _, _, patches = plt.hist(data, bins=bins, color='royalblue', ec='black', align='left') # rb
    for patch in patches:
        plt.setp(patch, alpha=1, lw=0.1)
    # plt.xticks(np.arange(minimum, maximum))
    plt.xlabel(x_label)
    plt.ylabel('Counts')
    plt.title(title)
    if save:
        path = os.path.join(config.PROJECT_DIR, 'scripts', 'figs', 'hist_rb.png')
        plt.savefig(path, dpi=2000)
    else:
        plt.show()


if __name__ == "__main__":
    data = open_file(os.path.join(config.PROJECT_DIR, 'scripts', 'data', 'rb.json'))
    # plot_data(data, 'Molecular Weight Distribution', 'Molecular Weight ($g \cdot mol^{-1}$)', True)
    plot_data(data, 'Rotatable Bonds Distribution', 'Number of Rotatable Bonds', False)
    # plot_data(data, 'TPSA Distribution', 'TPSA ($Å^{2}$)', False)
