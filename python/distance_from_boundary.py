#! /usr/bin/env python3

import math
from sys import argv, exit
import argparse

# Taken from http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

parser = argparse.ArgumentParser(usage = "%(prog)s [-h] file radius",
    description = "Calculate the distance between each atom and the center of the cylindrical grain boundary structure")
parser.add_argument("file", help = "The structure file to examine")
parser.add_argument('-o','--output', help = "Name of the output file to write to", default = 'distances.txt')
parser.add_argument('-r','--radius', type = float, default = 100.0, help = "The radius of the cylindrical grain")
#TODO: add optional arguments for atom types to examine/ignore

args = parser.parse_args()

with open(args.file) as f:
    header = f.readline()
    reading_atoms = False

    while not reading_atoms:
        line = f.readline()
        ls = line.split()
        if "xlo" in line:
            xlo,xhi = float(ls[0]), float(ls[1])
        if "ylo" in line:
            ylo,yhi = float(ls[0]), float(ls[1])
        if "zlo" in line:
            zlo,zhi = float(ls[0]), float(ls[1])
        if "atoms" in line:
            N = int(ls[0])
        if "Atoms" in line:
            reading_atoms = True
            #TODO: add parsing rules for LAMMPS output formats (e.g. atomic, charge, etc.) as well as custom formats (e.g. 'id type q xu yu zu')

    atoms = [None] * N
    i = 0
    if reading_atoms:
        for line in f:
            if not line or line == "\n":
                continue
            ls = line.split()
            id = int(ls[0])
            atoms[i] = [float(k) if j > 1 else int(k) for j,k in enumerate(ls)]
            atoms[i][2] -= xlo
            atoms[i][3] -= ylo
            atoms[i][4] -= zlo
            i += 1

center = [(xhi-xlo)/2.0,(yhi-ylo)/2.0]
distances = [None] * N
for index,value in enumerate(atoms,0):
    # TODO: use the parsing rules established above (e.g. line 50 TODO) to calculate the (2D) distance
    dist = math.sqrt((float(atoms[index][2]) - center[0])**2 + (float(atoms[index][3]) - center[1])**2)
    distances[index] = dist - args.radius

with open(args.output, 'w') as f:
    for i,j in zip(atoms,distances):
        f.write(' '.join(str(k) for k in flatten([i,j])))
        f.write("\n")
