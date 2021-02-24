#! /usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import math
import numpy as np

def periodCalc(m,k,ratio):
    try:
        tmp = len(ratio)
    except TypeError:
        ratio = [ratio]
    if any([i > 1 or i < 0 for i in ratio]):
        raise ValueError("number outside of valid range")
    if m <= 0 or k <= 0:
        raise ValueError("invalid value")
    convert = 1.0364364e-28 # converts from g eV /mol*Angstrom^2 to seconds^2
    return [2 * math.pi * np.sqrt(r * (1 - r) * m / k * convert) for r in ratio]

parser = argparse.ArgumentParser(usage = '%(prog)s [-h]', description = "Calculator to determine the best mass split ratio for the core-shell model")
parser.add_argument('-m', '--mass', nargs = '+', type = float, help = "The total mass(es) of each particle")
parser.add_argument('-k', '--spring', nargs = '+', type = float, help = "The spring constant(s) of each core-shell system")
parser.add_argument('--elements', nargs = '+', help = "The names of the elements")
parser.add_argument('--plot', action = 'store_true', help = "Plot the period of each particle")

args = parser.parse_args()

if not args.mass:
    args.mass = []
    tmp = 1
    n = 1
    while tmp > 0:
        try:
            tmp = float(input(f"Enter the mass of particle {n} (<=0 to quit): "))
            if tmp <= 0:
                break
            args.mass.append(tmp)
            n += 1
        except ValueError:
            print("Input must be a number")

if not args.spring:
    args.spring = []
    tmp = 1
    n = 0
    while len(args.spring) < len(args.mass):
        try:
            tmp = float(input(f"Enter the spring constant for particle {n+1} (<=0 to quit): "))
            if tmp <= 0:
                break
            args.spring.append(tmp)
            n += 1
        except ValueError:
            print("Input must be a number")

if not args.elements or len(args.elements) != len(args.mass):
    if not args.elements:
        args.elements = []
    while len(args.elements) != len(args.mass):
        args.elements.append("None")

even_split_vals = [0] * len(args.mass)
for i,mk in enumerate(zip(args.mass, args.spring)):
    even_split_vals[i] = periodCalc(mk[0], mk[1], 0.5)[0]

min_val_index = np.argmin(even_split_vals) # get the index of the minimum value: all other values are calculated based on this one

splits = [0] * len(args.mass)
m_min = args.mass[min_val_index]
k_min = args.spring[min_val_index]

for i,mk in enumerate(zip(args.mass, args.spring)):
    if i == min_val_index:
        splits[i] = 0.5
    else:
        splits[i] = 0.5 - 0.5 * np.sqrt(1 - (m_min * args.spring[i]) / (k_min * args.mass[i]))

print(f"For a period of {min(even_split_vals):.4e} seconds ({min(even_split_vals) * 1e12:.4f} ps):")
print(f"Recommended timestep: {np.floor(min(even_split_vals)*1e11*1000)/1000:.3f} ps")
for i,val_el in enumerate(zip(splits, args.elements)):
    val = val_el[0]
    el = i + 1 if val_el[1] == "None" else val_el[1]
    print(f"The ideal value for shell mass of element {el} is {val:.3f} (core = {(1-val)*args.mass[i]:.3f}, shell = {val*args.mass[i]:.3f})")

if args.plot:
    for m,k,el in zip(args.mass, args.spring, args.elements):
        x_vals = np.arange(0, 1.0001, 0.0001)
        y_vals = periodCalc(m,k,x_vals)
        plt.plot(x_vals, y_vals, label = el)
    plt.legend(loc='best')
    plt.show()
