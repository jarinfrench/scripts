#! /usr/bin/env python
"""Quick script to summarize the 3D MSD/diffusion data for defect diffusion in UMo"""
import itertools

import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from numpy.lib.recfunctions import append_fields

# cs = [
#     str(i) + "at%"
#     for i in [
#         1,
#         3,
#         5,
#         7,
#         10,
#         12,
#         14,
#         17,
#         19,
#         21,
#         25,
#         30,
#         35,
#         40,
#         45,
#         50,
#         55,
#         60,
#         65,
#         70,
#         75,
#         80,
#         85,
#         90,
#         95,
#         100,
#     ]
# ]
cs = [f"{i:.2f}at%Xe" for i in [0.25, 0.50, 0.75, 1.00]]
# cMos = [str(i) + "at%Mo" for i in [1, 3, 5]]
# cXes = [i + "at%Xe" for i in ["0.25", "0.50", "0.75"]]
Ts = ["T3000"]
dirs = ["dir_1", "dir_2", "dir_3"]

dirs_to_analyze = list(itertools.product(cs, Ts, dirs))

with open("diffusion_results_UXe.txt", "w", encoding="utf-8") as diff_file:
    diff_file.write(
        "#'concentration (at%)' 'temperature (K)' 'run' 'D_xyz (m^2/s)' 'std_err'\n"
    )
    for c, T, d in dirs_to_analyze:
        msd = np.recfromcsv(
            "/".join([c, T, d, "MSD_UXe.dat"]),
            delimiter=" ",
            dtype=[
                ("timestep", int),
                ("msd_x", float),
                ("msd_y", float),
                ("msd_z", float),
                ("msd_xyz", float),
            ],
            skip_header=2,
            names=None,
        )
        time = msd["timestep"] * 0.002
        msd = append_fields(msd, "time", time, usemask=False)
        x = msd["time"]
        y = msd["msd_xyz"]
        X = sm.add_constant(x)
        model = sm.OLS(y, X)
        result = model.fit()
        m = result.params[1]
        m_err = result.bse[1]
        plt.plot(x, y, ".", label="data")
        plt.plot(x, result.fittedvalues, "r-", label="fit")
        plt.legend()
        plt.savefig("/".join([c, T, d, "MSD_plot.png"]))
        plt.clf()
        diff_file.write(
            f"{c.strip('at%')} {T.strip('T')} {d.strip('dir_')} {m/6*1e-8} {m_err}\n"
        )
