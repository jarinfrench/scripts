#! /usr/bin/env python3

import math
import matplotlib.pyplot as plt
from sys import argv, exit
import argparse

parser = argparse.ArgumentParser(usage = "%(prog)s [-h] file radius",
    description = "Plot the distance between each atom and the center of the cylindrical grain boundary structure")
parser.add_argument("file", help = "The structure file to examine")
parser.add_argument("radius", type = float, default = 100.0, help = "The radius of the cylindrical grain")

args = parser.parse_args()

with open(args.file) as f:
    header = f.readline()
    line = f.readline() # blank line
    line = f.readline().split()
    N = int(line[0])
    line = f.readline() # atom types
    line = f.readline().split() # x bounds
    xlo,xhi = float(line[0]),float(line[1])
    line = f.readline().split() # y bounds
    ylo,yhi = float(line[0]),float(line[1])
    line = f.readline().split() # z bounds
    zlo,zhi = float(line[0]),float(line[1])
    line = f.readline() # another blank line
    line = f.readline() # Atoms
    line = f.readline() # yet another blank line

    atoms = [None] * N
    for line in f:
        id = int(line.split()[0])
        atoms[id - 1] = line.split()

center = [(xhi-xlo)/2.0,(yhi-ylo)/2.0]
distances = [None] * N
for index,value in enumerate(atoms,0):
    if not int(atoms[index][1]) == 1:
        continue
    dist = math.sqrt((float(atoms[index][3]) - center[0])**2 + (float(atoms[index][4]) - center[1])**2)
    distances[index] = dist - args.radius

x = [int(x[0]) for x in atoms]
plt.scatter(x,distances)
plt.show()
