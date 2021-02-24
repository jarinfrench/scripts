#! /usr/bin/env python3

from sys import argv, exit
import os
import argparse

def check_range(val):
    try:
        value = float(val)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err))

    if value < 0 or value > 3300:
        message = "Temperature out of fitted range. Expected 0 <= T <= 3300, got value = {}".format(value)
        raise argparse.ArgumentTypeError(message)

    return value

def printFits(fits):
    print("There are {l} fits:".format(l = len(fits)))
    for i,fit in enumerate(fits):
        print("  {} - {}".format(i + 1, fits[i][0]))

parser = argparse.ArgumentParser(usage = '%(prog)s temperature concentration potential', description = "Get the lattice parameter based on T given the lattice parameter database file")
parser.add_argument('T', type = float, help = "Desired temperature")
parser.add_argument('c', type = float, help = "Desired concentration")
parser.add_argument('p', type = int, help = "The desired potential to use")
parser.add_argument('--database', help = "Database file location", default = os.environ['HOME']+'/projects/scripts/lattice_params_2.db')

args = parser.parse_args()

fits = []

with open(args.database) as f:
    for _ in range(5):
        next(f)
    for line in f:
        fits.append(line.split())

if args.p > len(fits) or args.p - 1 < 0:
    print("Invalid potential specifier.")
    printFits(fits)
    exit(2)


if args.c >= float(fits[args.p - 1][1]) and args.c <= float(fits[args.p - 1][2]) and args.T >= float(fits[args.p - 1][3]) and args.T <= float(fits[args.p - 1][4]):
    x3 = float(fits[args.p - 1][5])
    x2 = float(fits[args.p - 1][6])
    x1 = float(fits[args.p - 1][7])
    y3 = float(fits[args.p - 1][8])
    y2 = float(fits[args.p - 1][9])
    y1 = float(fits[args.p - 1][10])
    xy = float(fits[args.p - 1][11])
    xy2 = float(fits[args.p - 1][12])
    x2y = float(fits[args.p - 1][13])
    z0 = float(fits[args.p - 1][14])
else:
    print("Error! Temperature or concentration not in fitted range.")
    exit(10)

a0 = (x3 * args.c**3 + y3 * args.T**3 +
     x2 * args.c**2 + y2 * args.T**2 +
     x1 * args.c    + y1 * args.T    +
     xy * args.c * args.T +
     xy2 * args.c * args.T**2 +
     x2y * args.c**2 * args.T + z0)
print(u"The lattice parameter using the {} potential at {:.1f} K and concentration {:.2f} is {:.5f} Angstroms".format(fits[args.p - 1][0], args.T, args.c, a0))
