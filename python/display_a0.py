#! /usr/bin/env python3

import sympy as sp
from thefuzz import fuzz
import argparse
import os

sp.init_printing(use_unicode = True)

def rlen(lst):
    return range(len(lst))

def printFits(fits):
    print("There are {l} fits:".format(l = len(fits)))
    for i,fit in enumerate(fits):
        print("  {} - {}".format(i + 1, fits[i][0]))

def isfloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

parser = argparse.ArgumentParser(usage = '%(prog)s potential', description = "Display the lattice parameter as a function of T and impurity concentration")
parser.add_argument('p', nargs = '+', help = "The desired potential from the database (integer or string)") # can either be a number or a string
parser.add_argument('-T', '--temperature', type = float, help = "Temperature")
parser.add_argument('-c', '--concentration', type = float, help = "Concentration")
parser.add_argument('--database', help = "Database file location", default = os.environ['HOME']+'/projects/scripts/lattice_params.db')
parser.add_argument('--precision', type = int, default = 5, help = "The precision to use for float values (default: 5)")

args = parser.parse_args()

fits = []
potential_index = None
args.p = '_'.join(args.p)

with open(args.database) as f:
    for _ in range(5):
        next(f)
    for line in f:
        fits.append([sp.Float(float(i),args.precision) if isfloat(i) else i for i in line.split()])

if args.p.isdecimal() and (int(args.p) > len(fits) or int(args.p) - 1 < 0):
    print("Invalid potential specifier.")
    printFits(fits)
    exit(2)
elif args.p.isdecimal():
    potential_index = int(args.p) - 1
    print(f"Using potential {fits[potential_index][0]}")
else:
    distances = [fuzz.partial_ratio(args.p, fits[i][0]) for i in range(len(fits))]
    potential_index = max(range(len(distances)), key=distances.__getitem__)
    print(f"Using potential {fits[potential_index][0]}")

fit = fits[potential_index]
T, c = sp.symbols('T c')


expr = c**3 * fit[5] + c**2 * fit[6] + c * fit[7] + T**3 * fit[8] + T**2 * fit[9] + T * fit[10] + c*T * fit[11] + c * T**2 * fit[12] + c**2 * T * fit[13] + fit[14]
c_domain = sp.Rel(fit[1], c, '<=') & sp.Rel(c, fit[2], '<=')
T_domain = sp.Rel(fit[3], T, '<=') & sp.Rel(T, fit[4], '<=')

if args.temperature and not T_domain.subs(T, args.temperature):
    print(f"Temperature {args.temperature} is not in the expected temperature domain:\n  {T_domain}")
    exit(1)
elif args.temperature:
    expr_T_sub = expr.subs(T, args.temperature)
    print(f"\nLattice parameter as a function of concentration at T = {args.temperature}:\n  {expr_T_sub}")

if args.concentration and not c_domain.subs(c, args.concentration):
    print(f"Concentration {args.concentration} is not in the expected concentration domain:\n  {c_domain}")
    exit(1)
elif args.concentration:
    expr_c_sub = expr.subs(c, args.concentration)
    print(f"\nLattice parameter as a function of Temperature at c = {args.concentration}:\n  {expr_c_sub}")


if args.concentration and args.temperature:
    expr_T_c_subs = expr.evalf(5, subs={T: args.temperature, c: args.concentration})
    print(f"\nLattice parameter at T = {args.temperature} and c = {args.concentration}:\n  {expr_T_c_subs}")

print(f"\nLattice parameter as a function of temperature (T) and concentration (c):\n  {expr}")
print(f"    Temperature domain: {fit[3]} <= T <= {fit[4]}")
print(f"    Concentration domain: {fit[1]} <= c <= {fit[2]}")
