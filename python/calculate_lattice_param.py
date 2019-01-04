#! /usr/bin/env python3

from natsort import natsorted
from math import ceil
from sys import exit
import matplotlib.pyplot as plt
import numpy as np
import os, glob, csv, argparse

def find(lst, elem):
    for i, x in enumerate(lst):
        if x == elem:
            return i
    return None

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] num_data_sets element', description = "Calculate the temperature-based lattice parameter based on multiple sets of data")
parser.add_argument('num_data_sets', type = int, help = "The number of subdirectories from the current directory containing the data.")
parser.add_argument('element', help = "The element/molecule being examined")
parser.add_argument('--non-cubic', action = 'store_true', help = "Flag specifying the system is not cubic (three separate lattice parameters)")
parser.add_argument('-t', '--stabilized-at', type = int, help = "The timestep at which the temperature stabilizes", default = 2000)
parser.add_argument('-b', '--file-name-base', help = "The base name of the files containing the data.  This program assumes a format of <file_name_base>_*.csv", default = "parsed_output")

args = parser.parse_args()

num_cells_x = int(input("Please enter the number of cells in the x direction: "))
if args.non_cubic:
    num_cells_y = int(input("Please enter the number of cells in the y direction: "))
    num_cells_z = int(input("Please enter the number of cells in the z direction: "))

lattice_param_x = float(input("Please enter the x lattice parameter at its lowest stable temperature (or 0 K): "))
if args.non_cubic:
    lattice_param_y = float(input("Please enter the y lattice parameter at its lowest stable temperature (or 0 K): "))
    lattice_param_z = float(input("Please enter the z lattice parameter at its lowest stable temperature (or 0 K): "))

steady = ceil(args.stabilized_at / 100) * 100 # determines the nearest 100 based on the input.  Used for averaging

curr_dir = os.getcwd() # gets the current working directory

average_list = []

for i in range(args.num_data_sets):
    os.chdir(curr_dir + "/dir_" + str(i+1))
    files_list = natsorted(glob.glob(args.file_name_base + "_*.csv"))
    k = 0
    if files_list[0] == args.file_name_base + "_0.csv":
        files_checked = files_list[1:] # we ignore the 0K case, because the simulations are run differently
    else:
        files_checked = files_list[:]

    for j in files_checked:
        with open(j) as fopen:
            reader = csv.reader(fopen)
            data = []
            for row in reader:
                if not row:
                    continue

                if row[0] == "Step":
                    data_types = row
                elif row[0].startswith(tuple('0123456789')): # note that the actual data starts with a number based on the file formatting
                    data.append(row)

            step_index = find(data_types, "Step")
            temp_index = find(data_types, "Temp")
            lx_index = find(data_types, "Lx")
            if args.non_cubic:
                ly_index = find(data_types, "Ly")
                lz_index = find(data_types, "Lz")
                if ly_index is None or lz_index is None:
                    print("Unable to continue: unknown box size")
                    exit(4)

            x_data = [float(val[step_index]) for val in data]
            y_data = [float(val[temp_index]) for val in data]

            plt.figure(1)
            plt.plot(x_data,y_data)
            plt.xlabel("Step")
            plt.ylabel("Temperature")
            plt.title("Temperature Analysis")

            x_data = [float(val[step_index]) for val in data if float(val[step_index]) >= steady]
            y_x_data = [float(val[lx_index]) for val in data if float(val[step_index]) >= steady]
            if args.non_cubic:
                y_y_data = [float(val[ly_index]) for val in data if float(val[step_index]) >= steady]
                y_z_data = [float(val[lz_index]) for val in data if float(val[step_index]) >= steady]

            plt.figure(2)
            plt.plot(x_data, y_x_data)
            plt.xlabel("Step")
            plt.ylabel("Lx")
            if args.non_cubic:
                plt.figure(3)
                plt.plot(x_data, y_y_data)
                plt.xlabel("Step")
                plt.ylabel("Lx")
                plt.figure(4)
                plt.plot(x_data, y_z_data)
                plt.xlabel("Step")
                plt.ylabel("Lx")

            average_lx = np.mean(y_x_data)
            if args.non_cubic:
                average_ly = np.mean(y_y_data)
                average_lz = np.mean(y_z_data)
                average_list.append([j, average_lx, average_ly, average_lz])
            else:
                average_list.append([j, average_lx])

    os.chdir(curr_dir)
plt.figure(2)
plt.axvline(x = steady, color = 'k')
plt.axhline(y = lattice_param_x * num_cells_x, color = 'k')
if args.non_cubic:
    plt.figure(3)
    plt.axvline(x = steady, color = 'k')
    plt.axhline(y = lattice_param_y * num_cells_y, color = 'k')
    plt.figure(4)
    plt.axvline(x = steady, color = 'k')
    plt.axhline(y = lattice_param_z * num_cells_z, color = 'k')

plt.show()
average_list = natsorted(average_list)
output = []

with open(args.element + "_A0_data.csv", 'w') as fout:
    csvwriter = csv.writer(fout, delimiter = ',')
    for i in range(0, len(average_list), args.num_data_sets):
        data_range = average_list[i:(i + args.num_data_sets)]
        mean_x = np.mean([j[1] for j in data_range]) / num_cells_x
        if args.non_cubic:
            mean_y = np.mean([j[2] for j in data_range]) / num_cells_y
            mean_z = np.mean([j[3] for j in data_range]) / num_cells_z
            output.append([mean_x, mean_y, mean_z])
            csvwriter.writerow([average_list[i][0], mean_x, mean_y, mean_z])
        else:
            output.append(mean_x)
            csvwriter.writerow([average_list[i][0], mean_x])
