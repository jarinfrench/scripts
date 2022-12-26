#! /usr/bin/env python3

from sys import argv, exit
import os
import argparse

def printFits(fits):
    print("There are {l} fits:".format(l = len(fits)))
    for i,fit in enumerate(fits):
        print("  {} - {}".format(i + 1, fits[i][0]))

parser = argparse.ArgumentParser(usage = '%(prog)s temperature potential', description = "Get the lattice parameter based on T given the lattice parameter database file")
parser.add_argument('T', type = float, help = "Desired temperature")
parser.add_argument('p', type = int, help = "The desired potential to use")
parser.add_argument('--database', help = "Database file location", default = os.environ['HOME']+'/projects/scripts/lattice_params.db')

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


if args.T >= float(fits[args.p - 1][1]) and args.T <= float(fits[args.p - 1][2]):
    A = float(fits[args.p - 1][3]);
    B = float(fits[args.p - 1][4]);
    C = 0.0;
elif args.T > float(fits[args.p - 1][2]) and args.T <= float(fits[args.p - 1][5]):
    A = float(fits[args.p - 1][6]);
    B = float(fits[args.p - 1][7]);
    C = float(fits[args.p - 1][8]);
else:
    print("Error! Temperature not in fitted range.")

print(u"The lattice parameter using the {} potential at {:.1f} K is {:.5f} Angstroms".format(fits[args.p - 1][0], args.T, A + B * args.T + C * args.T * args.T))
