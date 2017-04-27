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

# Now we get the labels
labels = next(reader)
data = []
for i in range(len(labels)): # for the total number of labels we have
    data.append([])

# Now we read through the data until either the end of the file, or until we
# come to another label set. <-- Second label set NOT IMPLEMENTED YET (TODO)
for row in reader:
    for i in range(len(labels)):
        data[i].append(float(row[i]))

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

for i in range(len(plotted)):
    x = plotted[i][0] - 1
    y = plotted[i][1] - 1
    title_label = labels[x] + " vs " + labels[y]
    x_label = labels[x]
    y_label = labels[y]
    plt.figure(i+1)
    plt.plot(data[x], data[y])
    plt.title(title_label)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

plt.show()
