#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from sys import argv, exit
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file[s]', description = "A script to plot grain sizes over time.")
parser.add_argument('files', nargs = '+', help = "The data files to plot")
parser.add_argument('-g','--grain-ids', type = int, nargs = '+', help = "The specific grain numbers to plot (indexed from 0)")
parser.add_argument('--dt', type = float, nargs = '+', help = "Optional timesteps for each file (multiplied by the step number to get a time value)")
parser.add_argument('--names', nargs = '+', help = "The legend names for each file specified, defaults to the file name")
parser.add_argument('-i','--ic', help = "The (optional) initial condition image to be embedded in the figure")
parser.add_argument('-s','--save', help = "Filename to save the figure as (saves as a png image)")
parser.add_argument('--no-display', action = 'store_true', help = "Flag to not display the figure")

args = parser.parse_args()

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

all_data = [None for i in range(len(args.files))]
color_index = 0 # index for current color
marker_index = 0 # index for current marker style
plot_grains = [] # the list of plots labeled with 'Grain {i}' (for identifying color trends)
plot_sets = [] # the list of ghost plots for labeling marker styles

for i,file in enumerate(args.files):
    all_data[i] = np.loadtxt(file, delimiter = ' ')

all_data = np.array(all_data)
for i in tqdm(range(len(args.files))):
    step = [j[0] for j in all_data[i]] # step is the first value in each row
    data = [j[1:] for j in all_data[i]] # the remaining lines are the sizes for each grain
    mark_every = int(len(step) / 20) # keep 20 markers to avoid clutter
    plt.xlabel("Step")
    plt.title("Area of all grains")
    if args.dt is not None:
        plt.xlabel("Time")
        step = [args.dt[i] * x for x in step]

    if args.grain_ids is not None:
        if len(args.grain_ids) > 1:
            plt.title("Area of grains {}".format(','.join([str(a) for a in args.grain_ids])))
        else:
            plt.title("Area of grain {}".format(','.join([str(a) for a in args.grain_ids])))
        data = [j[grain_ids] for j in all_data[i]]

    for j in range(len(data[0])): # each 'row' of data has the number of grains being plotted
        if i == 0:
            grain_num = args.grain_ids[j] if args.grain_ids is not None else j
            plot_grains.append(plt.plot(step, [k[j] for k in data], plot_colors[color_index] + marker_styles[marker_index] + '-', label = "Grain {}".format(grain_num), markevery = mark_every))
        else:
            plt.plot(step, [k[j] for k in data], plot_colors[color_index] + marker_styles[marker_index] + '-', markevery = mark_every)
        color_index += 1
    plot_sets.append(plt.plot([], [], 'k-' + marker_styles[marker_index], label = '{}'.format(args.names[i]), markevery = mark_every))
    marker_index += 1
    if marker_index >= len(marker_styles):
        marker_index = 0
    color_index = 0

# flatten the arrays
plot_sets = [item for sublist in plot_sets for item in sublist]
plot_grains = [item for sublist in plot_grains for item in sublist]
file_legend = plt.legend(handles = plot_sets, bbox_to_anchor = (1.04, 1), loc = 'upper left') # legend for the marker styles
plt.gca().add_artist(file_legend)
plt.legend(handles = plot_grains, bbox_to_anchor = (1.04, 0), loc = 'lower left', markerscale = 0) # legend for the grain colors
plt.ylabel("Area")
plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(16.0,9.0, forward = True)
if args.ic is not None: # embeds the initial condition image on the right side of the plot for easy comparison
    image = plt.imread(args.ic)
    newax = fig.add_axes([0.70, 0.25, 0.3, 0.5], anchor = 'E') # FIXME: These values need to change for every plot... is there a way to automatically determine them?
    newax.imshow(image)
    newax.axis('off')

if args.save is not None:
    fig.savefig(args.save, format = 'png')

if not args.no_display:
    plt.show()
