#! /usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import math
import numpy as np
from myModules import depth
from functools import partial

element_mass_dict = {
"H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011,
"N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
"Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
"K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
"Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
"Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
"Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95,
"Tc": 98, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41,
"In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29,
"Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
"Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50,
"Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49,
"Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
"Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209,
"At": 210, "Rn": 222, "Fr": 223, "Ra": 226, "Ac": 227, "Th": 232.04, "Pa": 231.04,
"U": 238.03, "Np": 237, "Pu": 244, "Am": 243, "Cm": 247, "Cf": 251, "Es": 252, "Fm": 257,
"Md": 258, "No": 259, "Lr": 266, "Rf": 267, "Db": 268, "Sg": 269, "Bh": 270, "Hs": 277,
"Mt": 278, "Ds": 281, "Rg": 282, "Cn": 285, "Nh": 286, "Fl": 289, "Mc": 290, "Lv": 293,
"Ts": 294, "Og": 294
}

def range_type(astr, min, max):
    value = float(astr)
    if not min < value <= max:
        raise argparse.ArgumentTypeError(f'Value not in range {min}, {max}')
    return value

def periodCalc(m,k,ratio):
    length = depth(ratio)
    if length == 1:
        ratio_list = [ratio]
    else:
        ratio_list = ratio
    if any([i > 1 or i < 0 for i in ratio_list]):
        raise ValueError("number outside of valid range")
    if m <= 0 or k <= 0:
        raise ValueError("invalid value")
    convert = 1.0364364e-28 # converts from g eV /mol*Angstrom^2 to seconds^2
    if length == 1:
        return 2 * math.pi * np.sqrt(ratio * (1 - ratio) * m / k * convert)
    else:
        return [2 * math.pi * np.sqrt(r * (1 - r) * m / k * convert) for r in ratio]

parser = argparse.ArgumentParser(usage = '%(prog)s [-h]', description = "Calculator to determine the best mass split ratio for the core-shell model")
parser.add_argument('-m', '--mass', nargs = '+', type = float, help = "The total mass(es) of each particle")
parser.add_argument('-k', '--spring', nargs = '+', type = float, help = "The spring constant(s) of each core-shell system")
parser.add_argument('-s', '--splits', nargs = '+', type = partial(range_type, min = 0, max = 0.5), help = "The percent of the mass assigned to the shell (0 < value <= 0.5)")
parser.add_argument('--elements', nargs = '+', help = "The names of the elements. Used to determine mass based on the most stable isotope if no masses are specified")
parser.add_argument('--plot', action = 'store_true', help = "Plot the period of each particle")

args = parser.parse_args()

if not args.mass:
    if args.elements:
        mass_tmp = []
        for el in args.elements:
            try:
                mass_tmp.append(element_mass_dict[el])
            except KeyError:
                prompt_anyway = True
                break
        else:
            args.mass = mass_tmp
            prompt_anyway = False
    else:
        prompt_anyway = True

    if prompt_anyway:
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
    even_split_vals[i] = periodCalc(mk[0], mk[1], 0.5)

min_val_index = np.argmin(even_split_vals) # get the index of the minimum value: all other values are calculated based on this one

splits = [0] * len(args.mass)
periods = [0] * len(args.mass)
m_min = args.mass[min_val_index]
k_min = args.spring[min_val_index]

for i in range(len(args.mass)):
    if i == min_val_index:
        splits[i] = 0.5
    else:
        splits[i] = 0.5 - 0.5 * np.sqrt(1 - (m_min * args.spring[i]) / (k_min * args.mass[i]))
    periods[i] = periodCalc(args.mass[i], args.spring[i], splits[i])

print(f"For a period of {min(even_split_vals):.4e} seconds ({min(even_split_vals) * 1e12:.4f} ps):")
print(f"Recommended timestep: {np.floor(min(even_split_vals)*1e15)/10000:.4f} ps")
for i,val_el in enumerate(zip(splits, args.elements)):
    val = val_el[0]
    el = i + 1 if val_el[1] == "None" else val_el[1]
    print(f"The ideal value for shell mass of element {el} is {val:.3f} (core = {(1-val)*args.mass[i]:.3f}, shell = {val*args.mass[i]:.3f}) (T = {periods[i]*1e15:.3f} fs)")

if args.splits:
    for val in args.splits:
        for i in range(len(args.mass)):
            el = i + 1 if args.elements[i] == "None" else args.elements[i]
            print(f"For the value of {val}, {el} has a period of {periodCalc(args.mass[i], args.spring[i], val)*1e15:.3f} fs")

if args.plot:
    for m,k,el in zip(args.mass, args.spring, args.elements):
        x_vals = np.arange(0, 1.0001, 0.0001)
        y_vals = periodCalc(m,k,x_vals)
        plt.plot(x_vals, [1e15*i for i in y_vals], label = el)
    plt.xlabel("Core Mass Percentage")
    plt.ylabel("Period (fs)")
    plt.legend(loc='best')
    plt.show()
