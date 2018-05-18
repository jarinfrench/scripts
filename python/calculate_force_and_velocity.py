#! /usr/bin/env python

from __future__ import division, print_function
from sys import exit, argv
from numpy import polyfit
import math, itertools, argparse

def threes(iterator):
    "s ->(s0,s1,s2),(s1,s2,s3),(s2,s3,s4), ..."
    a,b,c = itertools.tee(iterator,3)
    next(b,None)
    next(c,None)
    next(c,None)
    return zip(a,b,c)

def calculateLatticeParam(T, potential = 0):
    dbFile="/home/jarinf/projects/scripts/lattice_params.db"
    data = []
    with open(dbFile,'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            else:
                data.append(line)
    if potential == 0:
        print("There are %d fits:" %len(data))
        for i in range(len(data)):
            print("  {num} - {name}".format(num = i + 1, name = data[i].split()[0]))
        potential = int(input("Please specify the fit to use: "))

    while potential > len(data) or potential < 1:
        print("Invalid number. Additional potentials should be added to the database file.")
        potential = int(input("Please specify the fit to use: "))


    name, T1, yInt, slope, T2, yInt2, linC, paraC = data[potential - 1].split()
    print("Using {name} potential".format(name = name))

    if T < 0 or T > float(T2):
        print("Temperature out of fitted range.")
        exit(3)
    elif T >= 0 and T <= float(T1):
        return float(yInt) + float(slope) * T
    elif T > float(T1) and T <= float(T2):
        return float(yInt2) + float(linC) * T + float(paraC) * T**2
    else:
        print("Error calculating lattice parameter")


def calculateR (N,a0,Lz):
    return math.sqrt(N*a0**3/(4*math.pi*Lz))

def calculateForce(r):
    if r == 0:
        return 0
    else:
        return 1.6/r

def calculateVelocity(ns,ts,a0,Lz):
    coeff = math.sqrt(a0**3/(4*math.pi*Lz))
    if not len(ns) == 3 or not len(ts) == 3:
        print("Error calculating velocity of grain boundary")
    sqrtN = [math.sqrt(n) for n in ns]
    fit = polyfit(ts, sqrtN,1)
    return coeff * fit[0]


parser = argparse.ArgumentParser(description="Calculates the velocity and forces for an assumed cylindrical grain boundary given a data file in the format <timestep> <n grain 1> <n grain 2>")
parser.add_argument('t', metavar = 'T', type = float, help = "Temperature of the simulation")
parser.add_argument('l', metavar = 'Lz', type = float, help = "Thickness of the grain")
parser.add_argument('a', metavar = 'a0', type = float, help = "Lattice parameter at 0 K")
parser.add_argument('-p', '--potential', type = int, help = "Number of the potential to use from the database file", default = 0)

args=parser.parse_args()

# May change this to be an argument
dataFile = "data.txt" # data file containing the timestep, and number of atoms in each grain
dataOutfile = "force_velocity_data.txt"
a0=calculateLatticeParam(args.t, args.potential)

data = []
with open(dataFile,'r') as f:
    for _ in range(3):
        next(f)
    for line in f:
        data.append(line)

vel = []
force = []
r = []

for n, i in enumerate(threes(data)):
    ts  = [int(i[0].split()[0]) * 0.002, int(i[1].split()[0]) * 0.002, int(i[2].split()[0]) * 0.002]
    n1s = [float(i[0].split()[1]), float(i[1].split()[1]), float(i[2].split()[1])]
    n2s = [float(i[0].split()[2]), float(i[1].split()[2]), float(i[2].split()[2])]

    if n == 0:
        fit1 = polyfit(ts,n1s,1)
        fit2 = polyfit(ts,n2s,1)
        if fit1[0] > fit2[0]:
            use1 = False
        else:
            use1 = True

    if use1:
        r.append(calculateR(n1s[1], args.a, args.l))
        vel.append(calculateVelocity(n1s, ts, args.a, args.l))
    else:
        r.append(calculateR(n2s[1], args.a, args.l))
        vel.append(calculateVelocity(n2s, ts, args.a, args.l))

    force.append(calculateForce(r[n]))

with open(dataOutfile, 'w') as f:
    f.write("# Grain radius, Force on the grain, and instantaneous velocity of the boundary.\n")
    f.write("# Note that a negative velocity indicates grain shrinking.\n")
    for i in range(len(force)):
        f.write("{rad},{force},{vel}\n".format(rad=r[i], force=force[i], vel=vel[i]))
