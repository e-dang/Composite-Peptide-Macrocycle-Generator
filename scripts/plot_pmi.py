import matplotlib.pyplot as plt
import macrocycles.config as config
import os
import numpy as np


def open_file(filepath):
    with open(filepath, 'r') as file:
        for line in file.readlines():
            yield line.strip('\n').split(',')


def parse_data(lines):
    for line in lines:
        yield np.array(list(map(float, line[:-1])))


def plot(data, title, save=False):
    data = np.array(list(data))
    x_min, x_max = 0, 1
    y_min, y_max = 0.5, 1
    x_edge_width = 0.1
    y_edge_width = 0.1
    fig, ax = plt.subplots(1, 1)
    plot_triangle_overlay(ax)
    ax.scatter(data[:, 0], data[:, 1], s=0.2, c='b', marker='o', facecolors='full', linewidths=0)
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
        # path = os.path.join(config.PROJECT_DIR, 'scripts', 'figs', 'pmi_10k_p.png')
        path = os.path.join(config.PROJECT_DIR, 'scripts', 'figs', 'pmi_mil_b.png')
        plt.savefig(path, dpi=2000)
    else:
        plt.show()


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
