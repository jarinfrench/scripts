#! /usr/bin/env python

from __future__ import division,print_function
from natsort import natsorted
from math import ceil
from sys import exit
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons # look at https://stackoverflow.com/questions/6697259/interactive-matplotlib-plot-with-two-sliders to determine how to use sliders.
import numpy as np
import os
import glob
import csv

def find(lst, elem):
    for i, x in enumerate(lst):
        if x == elem:
            return i
    return None

curr_dir = os.getcwd() # gets the current working directory

# determine how many data sets we need
num_data_sets = int(input("Please enter the number of data sets to check: "))
base_file_name = "parsed_output"
element = raw_input("Please enter the element: ").title()
num_cells = int(input("Please enter the number of cells in the x direction: "))
steady = input("Please enter the timestep where the temperatures stabilizes: ")
steady = ceil(steady / 100) * 100
lattice_param = float(input("Please enter the lattice parameter at 0 K: "))

average_list = []

for i in range(num_data_sets):
    os.chdir(curr_dir + "/dir_" + str(i+1))
    files_list = natsorted(glob.glob(base_file_name + "_*.csv"))
    k = 0
    if files_list[0] == base_file_name + "_0.csv":
        files_checked = files_list[1:]
    else:
        files_checked = files_list[0:]
    for j in files_checked:
        with open(j) as fopen:
            reader = csv.reader(fopen)
            data = []
            for row in reader:
                if row[0]== "Step":
                    data_types = row;
                elif row[0].startswith(tuple('0123456789')):
                    data.append(row)
            step_index = find(data_types, "Step")
            temp_index = find(data_types, "Temp")
            lx_index = find(data_types, "Lx")
            x_data = [float(val[step_index]) for val in data]
            y_data = [float(val[temp_index]) for val in data]
            plt.figure(1)
            plt.plot(x_data,y_data)
            plt.xlabel("Step")
            plt.ylabel("Temperature")
            #plt.show(block = False)
            # Note that steps are in increments of 100, so we need to apply a ceiling to the nearest multiple of 100
            #plt.clf() # clear the figure
            x_data = [float(val[step_index]) for val in data if float(val[step_index]) >= steady]
            y_data = [float(val[lx_index]) for val in data if float(val[step_index]) >= steady]
            plt.figure(2)
            plt.plot(x_data,y_data)
            plt.xlabel("Step")
            plt.ylabel("Lx")
            #plt.axvline(x=steady, color='k')
            #plt.ylim(0,max(y_data)*1.2)
            #plt.xlim(0,max(x_data))
            #plt.show(block = False)
            #correct = raw_input("Enter y|Y if the plot is roughly straight after the vertical line, n|N otherwise: ")
            #plt.cla()
            #if any((c in ['n', 'N']) for c in correct):
            #    print("Error!")
            #    exit(3)
            average_lx = np.mean(y_data)
            average_list.append([j, average_lx])


    os.chdir(curr_dir)
plt.figure(2)
plt.axvline(x = steady, color = 'k')
plt.axhline(y = lattice_param * num_cells, color = 'k')
plt.show()
average_list = natsorted(average_list)
output = []
with open(element+"_a0_data.csv", 'w') as fout:
    csvwriter = csv.writer(fout, delimiter = ',')
    for i in range(0, len(average_list), num_data_sets):
        data_range = average_list[i:i+num_data_sets]
        mean = np.mean([j[1] for j in data_range]) / num_cells
        output.append(mean)
        csvwriter.writerow([average_list[i][0], mean])
