#! /usr/bin/env python3

import matplotlib.pyplot as plt
from sys import argv, exit
from os.path import basename, splitext
import csv
import argparse

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file[s]', description = "A simple script to plot two-dimensional data.")
parser.add_argument('files', nargs = '+', help = "The csv file(s) to plot")
parser.add_argument('-x', '--abscissa', type = int, help = "The column used for the x axis (0-based indexing). Default: 0", default = 0)
parser.add_argument('-y', '--ordinate', type = int, help = "The column used for the y axis (0-based indexing). Default: 1", default = 1)
parser.add_argument('--x-label', help = "The label of the x axis")
parser.add_argument('--y-label', help = "The label of the y axis")
parser.add_argument('--title', help = "The title of the plot")
parser.add_argument('-l', '--labels', nargs = '+', help = "The labels to place as legend entries. Default is the filename without the extension")
parser.add_argument('-p', '--plot-type', choices = ['s','l'], default = 's', help = "(s)catter plot or (l)ine plot")
parser.add_argument('-d','--delimiter', default = ',', help = "The delimiter between data values (default: \',\')")
parser.add_argument('-m', '--minimum', action = 'store_true', help = "Print the minimum data value of the plot")
parser.add_argument('-M', '--maximum', action = 'store_true', help = "Print the maximum data value of the plot")
parser.add_argument('-n', '--norm', choices = ['x','y','xy'], help = "Normalize the x, y, or both by the respective maximums")
parser.add_argument('--add-origin', action = 'store_true', help = "Adds the point (0,0) to the data")
parser.add_argument('--no-legend', action = 'store_true', help = "Flag to not put a legend on the plot")
parser.add_argument('--no-show', action = 'store_true', help = "Flag to not show the plot")
parser.add_argument('--abs', choices = ['x','y','xy'], help = "Take the absolute value of the x, y, or both data")
parser.add_argument('--legend-columns', type = int, default = 1, help = "The number of columns to split the legend into (default: 1)")
parser.add_argument('-s','--save', help = "Saves the plot with the specified file name")
args = parser.parse_args()

if not args.x_label:
    args.x_label = " "
if not args.y_label:
    args.y_label = " "
if not args.title:
    args.title = " "
if not args.labels:
    tmp = []
    for i in range(len(args.files)):
        tmp += [splitext(basename(args.files[i]))[0]]
    args.labels = tmp
if not len(args.labels) == len(args.files):
    # Using the file names for the unlabeled files
    tmp = args.labels
    for i in range(len(tmp), len(args.files)):
        tmp += [splitext(basename(args.files[i]))[0]]
    args.labels = tmp

plot_style = ['bo', 'go', 'ro', 'co', 'mo', 'yo', 'ko',
              'b1', 'g1', 'r1', 'c1', 'm1', 'y1', 'k1',
              'bv', 'gv', 'rv', 'cv', 'mv', 'yv', 'kv',
              'b^', 'g^', 'r^', 'c^', 'm^', 'y^', 'k^',
              'b<', 'g<', 'r<', 'c<', 'm<', 'y<', 'k<',
              'b>', 'g>', 'r>', 'c>', 'm>', 'y>', 'k>',
              'b2', 'g2', 'r2', 'c2', 'm2', 'y2', 'k2',
              'b3', 'g3', 'r3', 'c3', 'm3', 'y3', 'k3',
              'b4', 'g4', 'r4', 'c4', 'm4', 'y4', 'k4',
              'b8', 'g8', 'r8', 'c8', 'm8', 'y8', 'k8',
              'bs', 'gs', 'rs', 'cs', 'ms', 'ys', 'ks',
              'bp', 'gp', 'rp', 'cp', 'mp', 'yp', 'kp',
              'bP', 'gP', 'rP', 'cP', 'mP', 'yP', 'kP',
              'b*', 'g*', 'r*', 'c*', 'm*', 'y*', 'k*',
              'bh', 'gh', 'rh', 'ch', 'mh', 'yh', 'kh',
              'bH', 'gH', 'rH', 'cH', 'mH', 'yH', 'kH',
              'b+', 'g+', 'r+', 'c+', 'm+', 'y+', 'k+',
              'bx', 'gx', 'rx', 'cx', 'mx', 'yx', 'kx',
              'bX', 'gX', 'rX', 'cX', 'mX', 'yX', 'kX',
              'bD', 'gD', 'rD', 'cD', 'mD', 'yD', 'kD',
              'bd', 'gd', 'rd', 'cd', 'md', 'yd', 'kd',
              'b|', 'g|', 'r|', 'c|', 'm|', 'y|', 'k|',
              'b_', 'g_', 'r_', 'c_', 'm_', 'y_', 'k_']
x_max = 0.0
x_min = 1.0e20
y_max = 0.0
data_mins = []
data_maxs = []
data = []
for i in range(len(args.files)):
    if args.add_origin:
        x_data = [0.0]
        y_data = [0.0]
    else:
        x_data = []
        y_data = []
    fin = open(args.files[i])
    reader = csv.reader(fin, delimiter = args.delimiter)

    for j,line in enumerate(reader):
        if not (line) or line[0].startswith("#"):
            continue
        else:
            if len(line) == 1:
                x_data.append(j)
                y_data.append(float(line[0]))
            else:
                x_data.append(float(line[args.abscissa]))
                y_data.append(float(line[args.ordinate]))


    # Sort the data
    xy1 = sorted(zip(x_data, y_data))
    x_data = [x for x,y in xy1]
    y_data = [y for x,y in xy1]

    x_max = max(max(x_data), x_max)
    x_min = min(min(x_data), x_min)
    if args.norm:
        if 'x' in args.norm:
            x_data = [ii / x_max for ii in x_data]
        if 'y' in args.norm:
            y_max = max(max(y_data), y_max)
            y_data = [ii / y_max for ii in y_data]

    if args.abs:
        if 'x' in args.abs:
            x_data = [abs(ii) for ii in x_data]
        if 'y' in args.abs:
            y_data = [abs(ii) for ii in y_data]

    # recalculate these, as they may have changed
    x_max = max(max(x_data), x_max)
    x_min = min(min(x_data), x_min)

    # Plot the results
    if args.plot_type == 'l':
        plot_style[i] += '-'
    plt.plot(x_data, y_data, plot_style[i], label=args.labels[i], markersize = 0.5)

    if args.maximum:
        data_maxs.append(max(xy1, key = lambda item: item[1]))
    if args.minimum:
        data_mins.append(min(xy1, key = lambda item: item[1]))

if args.maximum:
    for i,item in enumerate(data_maxs):
        print(f"Max value for {args.files[i]}: {item}")
if args.minimum:
    for i,item in enumerate(data_mins):
        print(f"Min value for {args.files[i]}: {item}")


plt.xlabel(args.x_label)
plt.ylabel(args.y_label)
plt.xlim(x_min, x_max)
plt.title(args.title)
if not args.no_legend:
    plt.legend(loc='best', ncol = args.legend_columns, markerscale = 10)
if args.save:
    plt.savefig(args.save)
if not args.no_show:
    plt.show()
