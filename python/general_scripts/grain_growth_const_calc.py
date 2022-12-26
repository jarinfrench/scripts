#! /usr/bin/env python3

import argparse, sys
from math import log

def positiveFloat(num):
    try:
        num = float(num)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid float value {num}")
    if num <= 0:
        raise argparse.ArgumentTypeError("Value must be positive")
    else:
        return num


parser = argparse.ArgumentParser(usage = '%(prog)s [-h] d0 df t [options]', description = "Calculate the grain growth constant for a given initial and final grain size")
parser.add_argument('d0', type = positiveFloat, help = "Initial grain size (in microns)")
parser.add_argument('df', type = positiveFloat, help = "Final grain size (in microns)")
parser.add_argument('t', type = positiveFloat, help = "Time between start and end times (in hours)")
parser.add_argument('-g', '--gamma', type = positiveFloat, default = 1.0, help = "Grain boundary energy in J/m^2. Default = 1.0 J/m^2")
parser.add_argument('--ln', action = 'store_true', help = "Flag to show the natural log value of the calculated mobility")
parser.add_argument('--log', action = 'store_true', help = "Flag to show the log10 value of the calculated mobility")

args = parser.parse_args()

if args.df < args.d0:
    raise ValueError("Final grain size must be larger than initial grain size")

convert = 1e-12/3600 # converts microns^2/hr to m^2/s
K = (args.df**2 - args.d0**2)/args.t * convert
M = K / args.gamma

print(f"Mobility = {M} m^4/Js")
if args.ln:
    print(f"ln(M) = {log(M)}")
if args.log:
    print(f"log(M) = {log(M,10)}")
