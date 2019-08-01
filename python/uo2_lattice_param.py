#! /usr/bin/env python

from __future__ import division, print_function
from sys import argv, exit

if len(argv) == 1:
    T = float(input("Please enter the temperature in Kelvin: "))
else:
    T = float(argv[1])

if T > 0 and T <= 1000.0:
    A = 5.45297739;
    B = 5.89383E-5;
    C = 0.0;
elif T > 1000.0 and T <= 3300.0:
    A = 5.44364863;
    B = 6.18966E-5;
    C = 5.27784E-9;
else:
    print("Temperature out of fitted range (0 K - 3300 K).")
    exit(10)

print("The lattice parameter at %2.1f K is %7.5f Angstroms"%(T, A + B * T + C * T * T))
