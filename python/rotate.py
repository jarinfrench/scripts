#! /usr/bin/env python

from __future__ import division,print_function
from sys import argv
from math import sqrt


if len(argv) == 1:
    filename = raw_input("Please enter the filename to be processed: ")
else:
    filename = argv[1]

#filename2 = filename[:filename.find("_marked")] + "_flipped.dat"
filename2 = "flipped_data.dat"

f2 = open(filename2,"w")
atoms = []

with open(filename, "r") as f:
    for line in f:
        data = line.split()
        if len(data) != 6:
            continue
        if line.startswith("ITEM"):
            continue
        atoms.append([int(data[0]), int(data[1]), float(data[2]), float(data[3]), float(data[4]), float(data[5])])

xmin = min([x[2] for x in atoms])
ymin = min([y[3] for y in atoms])
zmin = min([z[4] for z in atoms])
xmax = max([x[2] for x in atoms])
ymax = max([y[3] for y in atoms])
zmax = max([z[4] for z in atoms])

xmid = (xmin + xmax) / 2.0
ymid = (ymin + ymax) / 2.0
zmid = (zmin + zmax) / 2.0

for i in range(len(atoms)):
    ytemp = atoms[i][3] - ymid
    ztemp = atoms[i][4] - zmid
    # rotate
    ytemp = -ytemp
    ztemp = -ztemp
    # assign back to the atom
    atoms[i][3] = ytemp + ymid
    atoms[i][4] = ztemp + zmid
    f2.write("%d %d %8.6f %8.6f %8.6f %d\n"%(atoms[i][0],atoms[i][1],atoms[i][2],atoms[i][3],atoms[i][4],atoms[i][5]))

f2.close()
