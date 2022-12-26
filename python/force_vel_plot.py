#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("force_velocity_data.txt", dtype = {
                "names": ("Time", "N", "r", "F", "vel"),
                "formats": ("float", "int", "float", "float", "float"),
                },
                skiprows = 2)

fig, ax = plt.subplots()
ax.scatter(-data["vel"], data["F"], color = '#1b9e77')
ax.set_xlabel("Velocity (m/s)")
ax.set_ylabel("Force (GPa)")

plt.savefig("force_vs_velocity.png")
