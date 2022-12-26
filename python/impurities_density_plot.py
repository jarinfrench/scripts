#! /usr/bin/env python3.9

import matplotlib.pyplot as plt
from scipy.stats import kde
import numpy as np
import argparse, sys

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file', description = 'Plot a density plot of the specified atom type from the LAMMPS dump file')
parser.add_argument("file", help = "The LAMMPS dump file containing the atom positions and types")
parser.add_argument('-t', '--type', type = int, nargs = '*', default = [2], help = "The atom type to plot")
parser.add_argument('-H', '--hist', action = 'store_true', help = "Flag to show a 2D histogram rather than a scatter plot")

args = parser.parse_args()

xs = []
ys = []

with open(args.file) as f:
    for _ in range(8):
        next(f)
    data_types = f.readline().split()[2:]
    type_index = data_types.index("type") if "type" in data_types else -1
    x_index = data_types.index("x") if "x" in data_types else -1
    y_index = data_types.index("y") if "y" in data_types else -1

    if x_index == -1:
        x_index = data_types.index("xu") if "xu" in data_types else -1
    if y_index == -1:
        y_index = data_types.index("yu") if "yu" in data_types else -1

    if -1 in [type_index, x_index, y_index]:
        print("One or more of type, x, or y index is missing")
        sys.exit(1)
    for line in f:
        sline = line.split()
        if int(sline[type_index]) in args.type:
            xs.append(float(sline[x_index]))
            ys.append(float(sline[y_index]))

# evaluate a Gaussian KDE (kernel density estimate) on a grid of 20x20 bins over data extents
# See https://python-graph-gallery.com/86-avoid-overlapping-in-scatterplot-with-2d-density
xs = np.array(xs)
ys = np.array(ys)
k = kde.gaussian_kde([xs,ys])
xi, yi = np.mgrid[xs.min():xs.max():20*1j, ys.min():ys.max():20*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))

fig,ax = plt.subplots()
# TODO: look into density plots that take into account PBCs
if args.hist:
    ax.set_title(f'2D histogram of atom type(s) {args.type}')
    h=ax.hist2d(xs,ys, bins=50, cmap=plt.cm.BuGn_r)
    fig.colorbar(h[3], ax=ax)
    # ax.pcolormesh(xi,yi,zi.reshape(xi.shape), shading = 'auto', cmap = plt.cm.BuGn_r)
else:
    ax.set_title(f'Scatter plot of atom type(s) {args.type}')
    ax.scatter(xs,ys)
    plt.legend()

plt.show()
