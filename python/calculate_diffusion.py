# /usr/bin/env python3

import os
import numpy as np
from numpy.lib.recfunctions import append_fields
from sklearn.linear_model import LinearRegression
import argparse

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file [file ...] [option(s)]', description = "Calculate the 1D, 2D, and 3D diffusion coefficients for a given MSD file")
parser.add_argument('file', nargs = '+', help = "The MSD file(s) to analyze")
parser.add_argument('-F', '--dump-freq', type = int, default = 10000, help = "How often dump files were written (default: 10000)")
parser.add_argument('--dt', type = float, default = 0.002, help = "Timestep size (default: 0.002 ps)")
parser.add_argument('--post-growth-diff', action = 'store_true', help = "Calculate the diffusion coefficients after the grain disappears separately")
parser.add_argument('--show', nargs = '+', choices = ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'], default = ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'], help = "Show the specified diffusion coefficients (default: all). Order matters")

args = parser.parse_args()

mobility_file = "slope_calc.txt"
diff_file = open("diffusion_results.txt", 'w')
if not os.path.exists(mobility_file):
    print("Error: slope_calc.txt not found. Use the mobility_plots.py script to generate it.")
    exit()

# if os.path.exists("MSD_U.dat"):
#     MSD_file = "MSD_U.dat"
#     header_skip = 3
# elif os.path.exists("A_MSD_U.dat"):
#     MSD_file = "A_MSD_U.dat"
#     header_skip = 2
# else:
#     print("Error: no MSD file for U found")
#     exit()

# Taken from https://www.codingem.com/how-to-read-the-last-line-of-a-file-in-python/
with open(mobility_file, "rb") as f:
    try:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
    except OSError:
        f.seek(0)
    idx_start, idx_end = f.readline().decode().split()[-1].split(':')
    idx_start = int(idx_start) - 1
    idx_end = int(idx_end) - 1

time_start = idx_start * args.dump_freq * args.dt
time_end = idx_end * args.dump_freq * args.dt
# time_start = (idx_start) * 10000 * 0.002 if "Basak" in os.getcwd().split('/') else (idx_start) * 20000 * 0.002
# time_end = (idx_end - 1) * 10000 * 0.002 if "Basak" in os.getcwd().split('/') else (idx_end) * 20000 * 0.002

for MSD_file in args.file:
    msd = np.recfromcsv(MSD_file, delimiter = " ", dtype = [('timestep', int), ('msd_x', float), ('msd_y', float), ('msd_z', float),('msd_xyz',float)], skip_header = 2, names = None)
    xy = msd['msd_x'] + msd['msd_y']
    xz = msd['msd_x'] + msd['msd_z']
    yz = msd['msd_y'] + msd['msd_z']
    time = msd['timestep'] * args.dt
    msd = append_fields(msd, 'time', time, usemask = False)
    msd = append_fields(msd, 'msd_xy', xy, usemask = False)
    msd = append_fields(msd, 'msd_xz', xz, usemask = False)
    msd = append_fields(msd, 'msd_yz', yz, usemask = False)

    diff_file.write(f"MSD for {MSD_file}:\n")
    active_grain_growth_time = (msd['time'] >= time_start) & (msd['time'] <= time_end)
    for i in args.show:
        factor = 2 if i in ['x','y','z'] else 4 if i in ['xy', 'xz', 'yz'] else 6
        x_aggt = msd['time'][active_grain_growth_time]
        y_aggt = msd['msd_' + i][active_grain_growth_time]
        model_aggt = LinearRegression().fit(x_aggt[:,None],y_aggt)
        r_sq_aggt = model_aggt.score(x_aggt[:,None],y_aggt)
        m_aggt = model_aggt.coef_[0]
        diff_file.write(f"  GB Diffusion for {i} = {m_aggt/factor*1e-8} m^2/s, r_sq = {r_sq_aggt}\n")

    if args.post_growth_diff:
        after_grain_disappears = msd['time'] > time_end
        if np.sum(after_grain_disappears) <= 2:
            diff_file.write("\n  Grain did not disappear entirely, no post-growth diffusion calculated\n")
            continue
        for i in args.show:
            factor = 2 if i in ['x','y','z'] else 4 if i in ['xy', 'xz', 'yz'] else 6
            x_agd = msd['time'][after_grain_disappears]
            y_agd = msd['msd_' + i][after_grain_disappears]
            model_agd = LinearRegression().fit(x_agd[:,None],y_agd)
            r_sq_agd = model_agd.score(x_agd[:,None],y_agd)
            m_agd = model_agd.coef_[0]
            diff_file.write(f"  Bulk diffusion for {i} = {m_agd/factor*1e-8} m^2/s, r_sq = {r_sq_agd}\n")
