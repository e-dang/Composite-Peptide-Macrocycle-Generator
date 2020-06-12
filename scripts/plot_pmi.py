import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import cpmg.config as config
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import numpy as np
from sklearn.neighbors import KernelDensity

cm = plt.cm.get_cmap('inferno')
# cm = cm.reversed()
# cm = mpl.colors.ListedColormap(name='inferno')
COLOR = 'purple'


def open_file(filepath):
    with open(filepath, 'r') as file:
        for line in file.readlines():
            yield line.strip('\n').split(',')


def parse_data(lines):
    for line in lines:
        yield np.array(list(map(float, line[:-1])))


def plot(data, title, save=False):
    data = np.array(list(data))
    # print(data)
    z = kde(data)
    # z_x = z[:, 0] / np.max(z[:, 0])
    # z_y = z[:, 1] / np.max(z[:, 1])
    # z = np.vstack([z_x, z_y])
    z = z / np.max(z)
    # print(z)
    norm = kde_color_weights(data)
    idx = z.argsort()
    x, y, z = data[:, 0][idx], data[:, 1][idx], z[idx]
    # for i in z:
    #     print(z)
    x_min, x_max = 0, 1
    y_min, y_max = 0.5, 1
    x_edge_width = 0.1
    y_edge_width = 0.1
    fig, ax = plt.subplots(1, 1)
    plot_triangle_overlay(ax)
    sc = ax.scatter(x, y, c=z, cmap=cm, s=0.1, marker='o', facecolors='full', linewidths=0, norm=norm)  # 1 mil
    # sc = ax.scatter(x, y, c=z, cmap=cm, s=4, marker='o', facecolors='full', linewidths=0, norm=norm)  # 10 k
    # sc = ax.scatter(data[:, 0], data[:, 1], c=z, cmap=cm, s=4, marker='o', facecolors='full', linewidths=0, norm=norm)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(sc, ax=ax, ticks=[-0.15, (1.0 - 0.15) / 2.0, 1], shrink=0.85)  # fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels(['Low Density', 'Medium Density', 'High Density'])
    ax.set_xlim((x_min - x_edge_width, x_max + x_edge_width))
    ax.set_ylim((y_min, y_max + y_edge_width))
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    # ax.grid(alpha=0.25)
    ax.set_xlabel('i1/i3')
    ax.set_ylabel('i2/i3')
    ax.set_title(title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_bounds(0.4, 1.0)
    ax.spines['bottom'].set_bounds(-0.1, 1.0)
    if save:
        # path = os.path.join(config.PROJECT_DIR, 'scripts', 'figs', 'pmi_10k_inf.png')
        path = os.path.join(config.PROJECT_DIR, 'scripts', 'figs', 'pmi_mil_inf.png')
        plt.savefig(path, dpi=2000)
    else:
        plt.show()


def kde_color_weights(data):
    # wt = len(data) / np.sum(data[:, 0])
    norm = mplc.Normalize(vmin=-0.15, vmax=1)
    return norm


def kde(data):
    kd = KernelDensity(kernel='tophat', bandwidth=0.02).fit(data)
    return kd.score_samples(data)


def plot_triangle_overlay(axes):
    top_left = [0, 1]
    top_right = [1, 1]
    bottom_middle = [0.5, 0.5]
    points = [top_left, top_right, bottom_middle]
    polygon = plt.Polygon(points, fill=False)
    axes.add_patch(polygon)


if __name__ == "__main__":
    # lines = open_file(os.path.join(config.PROJECT_DIR, 'scripts', 'data', 'pmi10k.csv'))
    lines = open_file(os.path.join(config.PROJECT_DIR, 'scripts', 'data', 'pmimil.csv'))
    coords = parse_data(lines)
    # plot(coords, 'PMI $10^3$ Macrocycles', True)
    plot(coords, 'PMI $10^6$ Macrocycles', True)
