#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

try:
    data_initial = np.loadtxt("rdf_initial.txt", dtype = {
        "names": ("separation", "1-1", "1-2", "1-3", "2-2", "2-3", "3-3"),
        "formats": ("float", "float", "float", "float", "float", "float", "float")},
        skiprows = 2)
    data_final = np.loadtxt("rdf_final.txt", dtype = {
        "names": ("separation", "1-1", "1-2", "1-3", "2-2", "2-3", "3-3"),
        "formats": ("float", "float", "float", "float", "float", "float", "float")},
        skiprows = 2)
except ValueError:
    data_initial = np.loadtxt("rdf_initial.txt", dtype = {
        "names": ("separation", "1-1", "1-2", "2-2"),
        "formats": ("float", "float", "float", "float")},
        skiprows = 2)
    data_final = np.loadtxt("rdf_final.txt", dtype = {
        "names": ("separation", "1-1", "1-2", "2-2"),
        "formats": ("float", "float", "float", "float")},
        skiprows = 2)

for field in data_initial.dtype.names[1:]:
    fig, ax = plt.subplots()
    ax.plot(data_initial["separation"], data_initial[field], color = '#1b9e77', label = "Initial")
    ax.plot(data_final["separation"], data_final[field], '--', color = '#d95f02', label = "Final")
    ax.set_xlabel("Separation ($\AA$)")
    ax.set_ylabel(f"g(r) [{field}]")
    plt.legend(loc = 'best')
    fig.savefig(f"rdf_{field}_comparison.png")
    plt.close()
