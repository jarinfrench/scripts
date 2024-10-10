#! /usr/bin/env python3

import xorbit.numpy as np
import os
import matplotlib.pyplot as plt

offset = 0.0
files = [
    f
    for f in os.listdir(".")
    if os.path.isfile(f) and f.startswith("cluster_results_type") and f.endswith("txt")
]
types = set([i.split("_")[2] for i in files])
ddtypes = {
    "names": ("id", "size", "CoM X", "CoM Y", "CoM Z"),
    "formats": ("int", "int", "float", "float", "float"),
}


def getData(imp_type):
    if imp_type not in ["type2", "type3"]:
        print(f"Unrecognized atom type: {imp_type}")
        exit()
    return (
        np.loadtxt(
            f"cluster_results_{imp_type}_d4.4_initial.txt", dtype=ddtypes, skiprows=2
        ),
        np.loadtxt(
            f"cluster_results_{imp_type}_d4.4_final.txt", dtype=ddtypes, skiprows=2
        ),
    )


for t in types:
    data_initial, data_final = getData(t)
    min_max = max(5, max(max(data_initial["size"]), max(data_final["size"])))
    h = {i: 0 for i in range(1, min_max + 1)}  # initial results
    h2 = {i + offset: 0 for i in range(1, min_max + 1)}  # final results
    for s in data_initial["size"]:
        h[s] += 1
    for s in data_final["size"]:
        h2[s + offset] += 1

    for i in range(1, 3):
        h.pop(i)
        h2.pop(i + offset)
    fig, ax = plt.subplots()
    ax.hist(
        [list(h.keys()), list(h2.keys())],
        bins=list(h.keys()),
        weights=[list(h.values()), list(h2.values())],
        color=["#1b9e77", "#d95f02"],
        label=["Initial", "Final"],
    )  # , alpha = 0.7)

    ax.set_xticks([i for i in range(1, min_max + 1)])
    ax.set_xlabel("Cluster size")
    ax.set_ylabel("Count")
    ax.set_ylim(0, max(1.2 * max(max(list(h.values())), max(list(h2.values()))), 1))
    ax.legend(loc="best")
    plt.savefig(f"{t}_cluster_histogram_d4.4.png")
