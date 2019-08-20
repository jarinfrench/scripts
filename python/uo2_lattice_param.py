#! /usr/bin/env python3

from sys import argv, exit
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

parser = argparse.ArgumentParser(usage = '%(prog)s T', description = "Get the lattice parameter of UO2 based on T (using the Basak potential)")
parser.add_argument('T', type = float, help = "Desired temperature")

args = parser.parse_args()

if args.T > 0 and args.T <= 1000.0:
    A = 5.45297739;
    B = 5.89383E-5;
    C = 0.0;
elif args.T > 1000.0 and args.T <= 3300.0:
    A = 5.44364863;
    B = 6.18966E-5;
    C = 5.27784E-9;

print("The lattice parameter at {:.1f} K is {:.5f} Angstroms"%(args.T, A + B * args.T + C * args.T * args.T))
