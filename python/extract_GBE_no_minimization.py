#! /usr/bin/env python3

from natsort import natsorted
from math import ceil
from sys import exit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons # look at https://stackoverflow.com/questions/6697259/interactive-matplotlib-plot-with-two-sliders to determine how to use sliders.
import numpy as np
import os, glob, csv, argparse
from myModules import *

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] timestep',
    description = "Calculate the approximate energy of a system that has not been minimized")
parser.add_argument('timestep', type = int, help = "The timestep at which the temperature has stabilized")
parser.add_argument('--end', type = int, help = "The timestep at which the stable temperature ends (default: last timestep)")

args = parser.parse_args()

# determine how many data sets we need\
steady = ceil(args.timestep / 100) * 100
files = natsorted(glob.glob("parsed_output_*.csv"))

average_list = []
n_list = []
for file in files:
    with open(file) as fopen:
        reader = csv.reader(fopen)
        data = []
        for row in reader:
            if row[0] == "Step":
                data_types = row;
                step_index = find(data_types, "Step")
                eng_index = find(data_types, "PotEng")
            if row[0] == "N =":
                n_list.append(int(row[1]))
            elif row[0].startswith(tuple('0123456789')):
                if args.end is not None:
                    if args.end < int(row[step_index]):
                        break
                data.append(row)

        x_data = [float(val[step_index]) for val in data]
        y_data = [float(val[eng_index]) for val in data]
        plt.figure(1)
        plt.plot(x_data,y_data)
        plt.xlabel("Step")
        plt.ylabel("Potential Energy")

        if args.end is not None:
            x_data = [float(val[step_index]) for val in data if float(val[step_index]) >= steady and float(val[step_index]) <= args.end]
            y_data = [float(val[eng_index]) for val in data if float(val[step_index]) >= steady and float(val[step_index]) <= args.end]
        else:
            x_data = [float(val[step_index]) for val in data if float(val[step_index]) >= steady]
            y_data = [float(val[eng_index]) for val in data if float(val[step_index]) >= steady]

        plt.figure(2)
        plt.plot(x_data,y_data)
        plt.xlabel("Step")
        plt.ylabel("Potential Energy")

        average_eng = np.mean(y_data)
        average_list.append([file, average_eng])

plt.show()
average_list = natsorted(average_list)
output = []
with open("GB_total_energies.csv", 'w') as fout:
    csvwriter = csv.writer(fout, delimiter = ',')
    for i in range(0, len(average_list)):

        csvwriter.writerow([average_list[i][0], average_list[i][1], n_list[i]])
