#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from sys import argv, exit
from tqdm import tqdm
from math import copysign
import argparse

sign = lambda x: copysign(1, x)

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file[s]', description = "A script to plot grain sizes over time.")
parser.add_argument('files', nargs = '+', help = "The data files to plot") # NOTE: will be output from the phase_field.cpp script as (generally) 'grain_size_distribution.txt'
parser.add_argument('-g','--grain-ids', type = int, nargs = '+', help = "The specific grain numbers to plot (indexed from 0)")
parser.add_argument('--dt', type = float, nargs = '+', help = "Optional timesteps for each file (multiplied by the step number to get a time value)")
parser.add_argument('--names', nargs = '+', help = "The legend names for each file specified, defaults to the file name")
parser.add_argument('-i','--ic', help = "The (optional) initial condition image to be embedded in the figure")
parser.add_argument('--normalize', action = 'store_true', help = "Normalize the area of each grain to it's initial value (i.e. all grain areas start at 1.0)")
parser.add_argument('--no-display', action = 'store_true', help = "Flag to not display the figure")
parser.add_argument('--growth-table', choices = ["end", "max"], default = argparse.SUPPRESS, help = "Shows a table showing the growth of each specified grain at the end of the simulation (end), or the maximum difference (max)")
parser.add_argument('-s','--save', help = "Filename to save the figure as (saves as a png image)")

args, unknown = parser.parse_known_args()

if args.growth_table is not None:
    from prettytable import PrettyTable
    pt = PrettyTable()

if args.grain_ids is not None:
    grain_ids = [i + 1 for i in args.grain_ids]

if args.dt is not None:
    if not len(args.dt) == len(args.files):
        print("The number of dt's specified must be equal to the number of files")
        exit(10)

if args.names is None:
    args.names = []
    for i in args.files:
        args.names.append(i)
else:
    if not len(args.names) == len(args.files):
        print("The number of legend names specified must be equal to the number of files")

plot_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'tab:orange', 'tab:purple', 'tab:brown',
               'tab:pink', 'tab:olive']
marker_styles = ['o', 'v', '^', 's', 'p', '*', 'P', 'x', '8', 'd']
line_styles = ['-', '--', ':']

all_data = [None for i in range(len(args.files))]
color_index = 0 # index for current color
marker_index = 0 # index for current marker style
plot_grains = [] # the list of plots labeled with 'Grain {i}' (for identifying color trends)
plot_sets = [] # the list of ghost plots for labeling marker styles
gb_lengths = []

for i,file in enumerate(args.files):
    all_data[i] = np.loadtxt(file, delimiter = ' ')

all_data = np.array(all_data)
if args.growth_table is not None:
    field_names = ["File"] + ["Grain {}".format(grain_num) if args.grain_ids is not None else "Grain {}".format(grain_num) for grain_num in range(len(all_data[0][0][1:]))]
    pt.field_names = field_names
    pt.align["File"] = "l"

for i in tqdm(range(len(args.files)), ncols = 80):
    step = [j[0] for j in all_data[i]] # step is the first value in each row
    data = np.array([j[1:] for j in all_data[i]]) # the remaining columns are the sizes for each grain
    if args.normalize:
        data = np.array([[data[j,k] / data[0,k] for k in range(len(data[0]))] for j in range(len(data))])
    mark_every = int(len(step) / 20) # keep 20 markers to avoid clutter
    plt.xlabel("Step", fontsize = 16)
    plt.title("Area of all grains", fontsize = 20)
    if args.dt is not None:
        plt.xlabel("Time (Arbitrary Units)", fontsize = 16)
        step = [args.dt[i] * x for x in step]

    if args.grain_ids is not None:
        if len(args.grain_ids) > 1:
            plt.title("Area of grains {}".format(','.join([str(a) for a in args.grain_ids])), fontsize = 20)
        else:
            plt.title("Area of grain {}".format(','.join([str(a) for a in args.grain_ids])), fontsize = 20)
        data = [j[grain_ids] for j in all_data[i]]

    if args.growth_table is not None:
        if len(args.names[i]) > 40:
            row = [args.names[i][0:15] + "..."]
        else:
            row = [args.names[i]]
        if args.growth_table == "end":
            row += ["{0:.3f}%".format((data[-1, j] - data[0, j]) / abs(data[0, j]) * 100) for j in range(len(data[0]))]
        else:
            # Complicated list comprehension that does the following:
            # max(range(len(list)), key = list.__getitem__) - returns the index of the list that has the max value
            # sign(value) - returns the sign (+1 or -1) of the value
            # max(abs(list) / abs(data[0][0])) - returns the maximum percent difference for each grain
            # NOTE: Still not quite correct. - e.g. no negative values shown, when there is clearly some negative growth
            row += ["{0:.3f}%".format(sign(data[max(range(len(abs(data[:,j] - data[0,j]))), key = data[:,0].__getitem__),j] - data[0,j]) * max(abs(data[:,j] - data[0,j])) / abs(data[0][j]) * 100.0) for j in range(len(data[0]))]
        pt.add_row(row)

    for j in range(len(data[0])): # each 'row' of data has the number of grains being plotted
        if i == 0:
            grain_num = args.grain_ids[j] if args.grain_ids is not None else j
            plot_grains.append(plt.plot(step, [k[j] for k in data], plot_colors[color_index] + marker_styles[marker_index] + '-', label = "Grain {}".format(grain_num), markevery = mark_every, lw = 2))
        else:
            plt.plot(step, [k[j] for k in data], plot_colors[color_index] + marker_styles[marker_index] + '-', markevery = mark_every, lw = 2)
        color_index += 1
    plot_sets.append(plt.plot([], [], 'k-' + marker_styles[marker_index], label = '{}'.format(args.names[i][0:15] + "..."), markevery = mark_every, lw = 2))
    marker_index += 1
    if marker_index >= len(marker_styles):
        marker_index = 0
    color_index = 0

# flatten the arrays
plot_sets = [item for sublist in plot_sets for item in sublist]
plot_grains = [item for sublist in plot_grains for item in sublist]
file_legend = plt.legend(handles = plot_sets, bbox_to_anchor = (1.04, 1), loc = 'upper left', fontsize = 14) # legend for the marker styles
plt.gca().add_artist(file_legend)
plt.legend(handles = plot_grains, bbox_to_anchor = (1.04, 0), loc = 'lower left', markerscale = 0, fontsize = 14) # legend for the grain colors
plt.ylabel("Area (Arbitrary Units)", fontsize = 16)
plt.gca().tick_params(axis='x', labelsize = 14)
plt.gca().tick_params(axis='y', labelsize = 14)
plt.tight_layout()
fig = plt.gcf()
fig.canvas.set_window_title('Grain Size Over Time')
fig.set_size_inches(16.0,9.0, forward = True)
if args.ic is not None: # embeds the initial condition image on the right side of the plot for easy comparison
    image = plt.imread(args.ic)
    newax = fig.add_axes([0.70, 0.25, 0.3, 0.5], anchor = 'E') # FIXME: These values need to change for every plot... is there a way to automatically determine them?
    newax.imshow(image)
    newax.axis('off')

if args.save is not None:
    fig.savefig(args.save, format = 'png')

if args.growth_table is not None:
    print(pt)

if not args.no_display:
    plt.show()
