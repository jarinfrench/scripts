#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt(
    "force_velocity_data.txt",
    dtype={
        "names": ("Time", "N", "r", "F", "vel"),
        "formats": ("float", "int", "float", "float", "float"),
    },
    skiprows=2,
)

fig, ax = plt.subplots()
ax.scatter(1 / data["r"], -data["vel"], color="#1b9e77")
ax.set_xlabel("Inverse Radius ($\305^{-1}$)")
ax.set_ylabel("Velocity (m/s)")

plt.savefig("force_vs_velocity.png")
