#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from sys import argv, exit
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file', description = "Python script to plot phase field results")
parser.add_argument('file', help = 'Input file containing the simulation details')
parser.add_argument('-i','--ic', default = 'initial_structure.txt', help = 'File containing the initial condition.  Defaults to \'initial_structure.txt\'')
parser.add_argument('--initial', action = 'store_true', help = "Flag to only show the IC")
parser.add_argument('-c','--cutoff', type = int, help = "The maximum timestep to show - default is the last timestep in the simulation")
parser.add_argument('--end', action = 'store_true', help = "Flag to only show the final time step")
parser.add_argument('-s','--save', action = 'store_true', help = "Flag to save the generated images")
parser.add_argument('--no-display', action = 'store_true', help = "Flag to not display the figures")

args = parser.parse_args()

file_data = np.loadtxt(args.file, dtype = str, usecols = 2, skiprows = 1)
grid_size = int(file_data[0])
dx = float(file_data[2])
num_grains = int(file_data[4])
num_steps = int(file_data[5])
n_step = int(file_data[6])
dt = float(file_data[7])
x_max = grid_size * dx
y_max = grid_size * dx

if args.cutoff is None:
    args.cutoff = num_steps

if args.no_display is None:
    while num_steps / n_step > 10:
        n_step += n_step

if args.end:
    n_step = args.cutoff

y, x = np.mgrid[slice(0, y_max, dx), slice(0,x_max, dx)]
z = np.loadtxt(args.ic) - 1

plt.figure()
c = plt.pcolormesh(x,y,z, cmap = 'jet', vmin = 0, vmax = num_grains - 1)
plt.colorbar(c, ax = plt.gca())
plt.title("Initial Condition")
plt.axis('square')

if args.save:
    plt.savefig('initial_structure.png')

if args.initial:
    plt.show()
    exit(0)

for i in tqdm(range(n_step, args.cutoff + 1, n_step), ncols = 80):
    file = 'out_' + str(i) + '.txt'
    z = np.loadtxt(file)
    plt.figure()
    c = plt.pcolormesh(x,y,z, cmap = 'gray', vmin = 0, vmax = 1)
    plt.colorbar(c, ax = plt.gca())
    plt.title("Timestep {} ({}s)".format(i, dt * i))
    plt.axis('square')

    if args.save:
        plt.savefig('out.' + str(i).zfill(len(str(num_steps))) + '.png')

    if args.no_display:
        plt.close()

if not args.no_display:
    plt.show()
