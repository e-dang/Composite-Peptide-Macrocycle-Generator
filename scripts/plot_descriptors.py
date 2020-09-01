import matplotlib as mpl
import json
from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import cpmg.config as config
from scipy.stats import iqr
import numpy as np
from numpy.random import rand
from collections import Counter


def open_file(filepath):
    with open(filepath, 'r') as file:
        for val in json.load(file):
            yield val


mpl.rcParams.update({'font.size': 15})
cm = plt.cm.get_cmap('inferno')
cm = cm.reversed()
COLOR = 'inf'


class MidpointNormalize(mpl.colors.Normalize):
    # class from the mpl docs:
    # https://matplotlib.org/users/colormapnorms.html

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def plot_mw(data, save=False):
    data = list(data)
    maximum = max(data)
    minimum = min(data)
    bin_width = np.round(2 * iqr(data) * (len(data) ** (-1 / 3)))
    bin_width = bin_width if bin_width > 1 else 1
    bins = int(round((maximum - minimum) / bin_width))
    vmax = 0.057 + 0.01
    normalizer = MidpointNormalize(vmin=0, vmax=vmax, midpoint=vmax / 2)

    fig = plt.figure(figsize=(15, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.03])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    #
    n, bins, patches = ax1.hist(data, bins=bins, color=COLOR, ec='black')
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm, orientation='vertical', norm=normalizer)
    cb1.set_ticks(cb1.get_ticks(), update_ticks=True)
    cb1.set_ticklabels(list(map(lambda x: str(round(x * 100)), cb1.get_ticks())))
    print(cb1.get_ticks())
    norm = np.array(n) / np.sum(n)
    print(norm, np.max(norm), np.sum(norm))
    norm = normalizer(norm)
    for i, patch in enumerate(patches):
        plt.setp(patch, alpha=1, lw=0.1)
        plt.setp(patch, 'facecolor', cm(norm[i]))

    # n, bins, patches = plt.hist(data, bins=bins, color=COLOR, ec='black')
    # for n_i, patch in zip(n, patches):
    #     plt.setp(patch, alpha=1, lw=0.1)
    #     plt.setp(patch, 'facecolor', cm(float(n_i) / float(max(bins))))
    # plt.xticks(np.arange(minimum, maximum + bin_width, step=bin_width))
    ax1.set_xticks(bins)
    ax1.set_xticklabels(map(lambda x: int(round(x)), bins))
    for tick in ax1.get_xticklabels():
        tick.set_rotation(45)
    ax1.set_xlabel('Molecular Weight ($g/mol$)')
    ax1.set_ylabel('Counts')
    ax1.set_title('Molecular Weight Distribution')
    ax2.set_ylabel('% Total')
    # ax1.tick_params(axis='x',which='minor',direction='out',bottom=True,length=5)
    plt.tight_layout()
    if save:
        path = os.path.join(config.PROJECT_DIR, 'images', f'hist_mw_{COLOR}_new.png')
        # plt.savefig(path, dpi=2000)
        fig.savefig(path, dpi=1000)
    else:
        plt.show()


def plot_tpsa(data, save=False):
    data = list(data)
    maximum = max(data)
    minimum = min(data)
    bin_width = 2 * iqr(data) * (len(data) ** (-1 / 3))
    bin_width = bin_width if bin_width > 1 else 1
    bins = int(round((maximum - minimum) / bin_width))

    # 0.073 is maximum normalized value
    vmax = 0.073 + 0.01
    normalizer = MidpointNormalize(vmin=0, vmax=vmax, midpoint=vmax / 2)
    fig = plt.figure(figsize=(15, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.03])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    n, bins, patches = ax1.hist(data, bins=bins, color=COLOR, ec='black')
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm, orientation='vertical', norm=normalizer)
    cb1.set_ticks(cb1.get_ticks(), update_ticks=True)
    cb1.set_ticklabels(list(map(lambda x: str(round(x * 100)), cb1.get_ticks())))
    print(cb1.get_ticks())
    norm = np.array(n) / np.sum(n)
    print(norm, np.max(norm), np.sum(norm))
    norm = normalizer(norm)
    for i, patch in enumerate(patches):
        plt.setp(patch, alpha=1, lw=0.1)
        plt.setp(patch, 'facecolor', cm(norm[i]))

    # n, bins, patches = plt.hist(data, bins=bins, color=COLOR, ec='black')
    # for n_i, patch in zip(n, patches):
    #     plt.setp(patch, alpha=1, lw=0.1)
    #     x = (n_i / max(n)) / 2
    #     plt.setp(patch, 'facecolor', cm(x))
    # plt.xticks(bins, rotation='vertical')
    # plt.xlabel('TPSA ($Å^{2}$)')
    # plt.ylabel('Counts')
    # plt.title('TPSA Distribution')

    ax1.set_xticks(bins)
    ax1.set_xlabel('TPSA ($Å^{2}$)')
    ax1.set_xticklabels(map(lambda x: int(round(x)), bins))
    for tick in ax1.get_xticklabels():
        tick.set_rotation(45)
    ax1.set_ylabel('Counts')
    ax1.set_title('TPSA Distribution')
    ax2.set_ylabel('% Total')
    plt.tight_layout()
    if save:
        path = os.path.join(config.PROJECT_DIR, 'images', f'hist_tpsa_{COLOR}_new.png')
        # plt.savefig(path, dpi=2000)
        fig.savefig(path, dpi=1000)
    else:
        plt.show()


def plot_rb(data, save=False):
    data = list(data)
    # maximum = max(data)
    # minimum = min(data)
    # bin_width = 2 * iqr(data) * (len(data) ** (-1 / 3))
    # bin_width = bin_width if bin_width > 1 else 1
    # bins = int(round((maximum - minimum) / bin_width))
    counter = Counter(data)
    bins = counter.keys()
    counts = np.array(list(counter.values()))
    norm = counts / np.sum(counts)
    normalizer = MidpointNormalize(vmin=0, vmax=0.2, midpoint=0.1)

    fig = plt.figure(figsize=(15, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.03])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    # n, bins, patches = ax1.hist(data, bins=bins, color=COLOR, ec='black', align='left')
    my_data = 5 * rand(5)
    bars = ax1.bar(counter.keys(), counter.values(), color=cm(normalizer(norm)), ec='black')
    # ax1.pcolormesh(counter.keys(), counter.values(), norm=MidpointNormalize(midpoint=0))
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm, orientation='vertical', norm=normalizer)
    cb1.set_ticks(cb1.get_ticks(), update_ticks=True)
    cb1.set_ticklabels(list(map(lambda x: str(round(x * 100)), cb1.get_ticks())))
    # cb1.set_clim((0, 0.5))
    # n, bins, patches = plt.hist(data, bins=bins, color=COLOR, ec='black', align='left')  # rb
    # for n_i, patch in zip(n, patches):
    #     plt.setp(patch, alpha=1, lw=0.1)
    #     x = (n_i / max(n)) / 2
    #     plt.setp(patch, 'facecolor', cm(x))
    ax1.set_xticks(list(bins))
    ax1.set_xlabel('Number of Rotatable Bonds')
    ax1.set_ylabel('Counts')
    ax2.set_ylabel('% Total')
    ax1.set_title('Rotatable Bonds Distribution')
    plt.tight_layout()
    if save:
        path = os.path.join(config.PROJECT_DIR, 'images', f'hist_rb_{COLOR}_new.png')
        fig.savefig(path, dpi=1000)
    else:
        plt.show()


if __name__ == "__main__":
    mw = open_file(os.path.join(config.PROJECT_DIR, 'scripts', 'data', 'mw.json'))
    tpsa = open_file(os.path.join(config.PROJECT_DIR, 'scripts', 'data', 'tpsa.json'))
    rb = open_file(os.path.join(config.PROJECT_DIR, 'scripts', 'data', 'rb.json'))
    # plot_mw(mw, True)
    plot_tpsa(tpsa, True)
    plot_rb(rb, True)
