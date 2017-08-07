#! /usr/bin/env python

from __future__ import division, print_function
from sys import argv
from collections import defaultdict
import matplotlib.pyplot as plt
import csv

# This script is should be utilized AFTER using the parse_lammps_output.cpp
# script.  This script reads in the data from that file, and then prompts
# the user as to which plots to show.

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs) > 1)

def depth(data):
    try:
        depth = len(data) # depth is the length of the list
    except:
        depth = 1 # or it's one if it's an integer or float
    return depth

def determine_label(label_set, label_dict, unit_dict):
    l = []
    for i in range(depth(label_set)):
        if label_set[i] in label_dict:
            l.append(unit_dict.get(label_dict.get(label_set[i])))
    if not l or l.count(l[0]) == len(l):
        if depth(l) == 1:
            return l
        else:
            return l[0]
    else:
        return "Unknown Units"


if len(argv) == 1:
    filename = input("Please enter the filename of the parsed output data in quotes: ")
else:
    filename = argv[1]

# Now that we have the filename, lets open it
fopen = open(filename, 'r')
reader = csv.reader(fopen)
line1 = next(reader) # Get the LAMMPS version

if (line1[0] != "LAMMPS version:"):
    print("Unable to determine LAMMPS version.")
    version = "Unknown"
else:
    version = line1[1]

line2 = next(reader) # Get the box size

# Double check that we have what we expect.  NOTE that this assumes a 3D system
if (line2[0] != "Lx Ly Lz ="):
    print("Unable to determine box size")
    Lx = -1.0
    Ly = -1.0
    Lz = -1.0
else:
    Lx = line2[1]
    Ly = line2[2]
    Lz = line2[3]

line3 = next(reader) # Get the number of atoms

# Again, double check that we have the right input.
if (line3[0] != "N ="):
    print("Unable to determine the number of atoms.")
    N = -1
else:
    N = int(line3[1])

line4 = next(reader) # get the units style

# Double check for correctness
if (line4[0] != "Unit style:"):
    print("Unable to determine unit style.")
    unit_style = "Unknown"
else:
    unit_style = line4[1]

if unit_style == "metal":
    unit_labels = {"mass": "g/mol", "distance": r"$\AA$", "time": "ps",
                   "energy": "eV", "velocity": r"$\AA$/ps", "force": r"eV/$\AA$",
                   "torque": "eV", "temperature": "K", "pressure": "bars",
                   "dynamic viscosity": "Poise", "charge": r"Multiples of e$^-1$ charge",
                   "dipole": r"charge*$\AA$", "electric field": r"V/$\AA$",
                   "volume": r"$\AA^3$",
                   "density": r"g/cm$^3$", "none": " ", "cpu": "s"}
elif unit_style == "si":
    unit_labels = {"mass": "kg", "distance": "m", "time": "s",
                   "energy":"J", "velocity": "m/s", "force": "N",
                   "torque": "N*m", "temperature": "K", "pressure": r"N/m$^2$",
                   "dynamic viscosity": r"N*s/m$^2$", "charge": "C",
                   "dipole": "C*m", "electric field": "V/m",
                   "volume": r"m$^3$",
                   "density": r"kg/m$^3$", "none": " ", "cpu": "s"}
else:
    print("Please see the LAMMPS manual for relevant units.")
    unit_labels = "Not implemented yet"

# Map from labels to units.
lammps_thermo = {"Step": "none", "Elapsed": "time", "Elaplong": "time",
                 "Dt": "time", "Time": "time", "CPU": "cpu", "Atoms": "none",
                 "Temp": "temperature", "Press": "pressure", "PotEng": "energy",
                 "KinEng": "energy", "TotEng": "energy", "Enthalpy": "none",
                 "Pxx": "pressure", "Pyy": "pressure", "Pzz": "pressure",
                 "Pxy": "pressure", "Pxz": "pressure", "Pyz": "pressure",
                 "Fmax": "force", "Fnorm": "force", "E_per_atom": "energy",
                 "Lx": "distance", "Ly": "distance", "Lz": "distance",
                 "Volume": "volume"}

# Now we get the labels
labels = next(reader)
data = []
for i in range(len(labels)): # for the total number of labels we have
    data.append([])

# Now we read through the data until either the end of the file, or until we
# come to another label set. <-- Second label set NOT IMPLEMENTED YET (TODO)
for row in reader:
    for i in range(len(labels)):
        try:
            data[i].append(float(row[i]))
        except:
            print("Second label set not implemented yet.")
            break;

# This removes duplicated values from the FIRST data set (the assumption is that
# the first column of data contains the step number)
for dup in sorted(list_duplicates(data[0]), reverse = True):
    for i in range(len(labels)):
        del data[i][dup[1][1]]

# Calculate the energy per atom
if "TotEng" in labels and N > 1:
    data.append([])
    totEngIndex = labels.index("TotEng")
    for i in range(len(data[totEngIndex])):
        data[len(labels)].append(data[totEngIndex][i] / N)
    labels.append("E_per_atom")

print("There are %d data sets:"%len(labels))
for i in range(len(labels)):
    print("%d - %s"%(i + 1, labels[i]))

print("Please enter the data you would like plotted against each other.")
plotted = input("Use tuple format (ex: [[x_data,y_data], [x_data2,y_data2]]): ")

data_len = []
for i in range(len(plotted)): # now we find the lengths of each data set being plotted
    data_len.append([i,[depth(plotted[i][0]), depth(plotted[i][1])]])

for i in range(len(plotted)):
    try:
        x = [j - 1 for j in plotted[i][0]] # index of x data and label
    except:
        x = plotted[i][0] - 1
    try:
        y = [j - 1 for j in plotted[i][1]] # index of y data and label
    except:
        y = plotted[i][1] - 1

    if not depth(y) == 1:
        title_label_main = "["
        for j in range(len(y)):
            title_label_main += labels[y[j]]
            if not j + 1 == len(y):
                title_label_main += ", "
            else:
                title_label_main += "] + vs "
    else:
        title_label_main = labels[y] + " vs "
    if not depth(x) == 1:
        title_label_main += "["
        for j in range(len(x)):
            title_label_main += labels[x[j]]
            if not j + 1 == len(x):
                title_label_main += ", "
            else:
                title_label_main += "]"
    else:
        title_label_main += labels[x]

    title_label_right = "LAMMPS version: " + version + "; N = " + str(N)
    if depth(x) > 1:
        labels_x = []
        for j in range(len(x)):
            labels_x .append(labels[x[j]])
        x_label = determine_label(labels_x, lammps_thermo, unit_labels)
    else:
        x_label = labels[x] + " (" + determine_label([labels[x]], lammps_thermo, unit_labels)[0] + ")"

    if depth(y) > 1:
        labels_y = []
        for j in range(len(y)):
            labels_y .append(labels[y[j]])
        y_label = determine_label(labels_y, lammps_thermo, unit_labels)
    else:
        y_label = labels[y] + " (" + determine_label([labels[y]], lammps_thermo, unit_labels)[0] + ")"
    plt.figure(i+1)

    if depth(x) > 1 and depth(y) > 1:
        for j in range(depth(x)):
            for m in range(depth(y)):
                plt.plot(data[x[j]], data[y[m]], label = labels[y[m]])
    elif depth(x) > 1 and depth(y) == 1:
        for j in range(depth(x)):
            plt.plot(data[x[j]], data[y], label = labels[y])
    elif depth(x) == 1 and depth(y) > 1:
        for j in range(depth(y)):
            plt.plot(data[x], data[y[j]], label = labels[y[j]])
    else:
        plt.plot(data[x], data[y], label=labels[y])
    plt.title(title_label_main + "\n")
    plt.suptitle("\n\n" + title_label_right, fontsize=8)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if depth(x) > 1 or depth(y) > 1:
        plt.legend(loc='best')

plt.show()
