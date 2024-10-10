# coding: utf-8
import argparse
import glob
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from myModules import promptForContinue
from numpy.lib.recfunctions import append_fields

# mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['axes.linewidth'] = 2
# mpl.rcParams['xtick.labelsize'] = 'x-small'
# mpl.rcParams['ytick.labelsize'] = 'x-small'
#
# def publish(ax: plt.Axes):
#     ax.tick_params(axis = 'both', which = 'both', direction = 'in')
#     ax.tick_params(axis = 'both', which = 'major', length = 7)
#     ax.tick_params(axis = 'both', which = 'minor', length = 4)

plt.style.use("publication")
mpl.rcParams["lines.markersize"] = 1.5
mpl.rcParams["legend.handletextpad"] = 0.08
mpl.rcParams["legend.markerscale"] = 3
mpl.rcParams["scatter.marker"] = "."
mpl.rcParams["font.size"] = 14
mpl.rcParams["axes.labelsize"] = 14
mpl.rcParams["legend.fontsize"] = 12
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12


parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] base [options]",
    description="MSD plots for all temperatures/simulations for one misorientation.",
)
parser.add_argument("base", help="basename of the saved files")
parser.add_argument(
    "-T",
    type=int,
    nargs=2,
    help="Minimum and maximum temperatures to plot MSD of (default: all temperatures)",
    default=[0, 100000],
)
parser.add_argument(
    "--dt",
    type=float,
    help="The timestep used in the simulations",
    default=0.002,
)
parser.add_argument(
    "-n", type=int, nargs="+", help='"Dir" numbers to plot (default: all)'
)
parser.add_argument(
    "-i",
    action="store_true",
    help="Interactive plotting (WARNING: uses a lot of memory!)",
)
parser.add_argument(
    "-2",
    dest="two",
    action="store_true",
    help="Create the 2D MSD plot (xy, xz, and yz)",
)
parser.add_argument(
    "-3", dest="three", action="store_true", help="Create the 3D MSD plot (xyz)"
)

args = parser.parse_args()

rootdir = "."
Ts = []
for f in os.listdir(rootdir):
    d = os.path.join(rootdir, f)
    if os.path.isdir(d):
        val = int(f[1:])
        if val >= args.T[0] and val <= args.T[1]:
            Ts.append(f[1:])
Ts.sort()

if not args.n:
    dirs = list(
        set(os.path.basename(os.path.dirname(i)) for i in glob.glob("*/*/dir_*/"))
    )
    dir_nums = [int("".join(filter(str.isdigit, d))) for d in dirs]
    dir_nums.sort()
    args.n = dir_nums

if args.i:
    figs = dict()
    axs = dict()

for T in Ts:
    for n in args.n:
        if not os.path.isdir(f"T{T}/large_r/dir_final_{n}"):
            continue
        msd = np.recfromcsv(
            f"T{T}/large_r/dir_final_{n}/MSD_U.dat",
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
        xy = msd["msd_x"] + msd["msd_y"]
        xz = msd["msd_x"] + msd["msd_z"]
        yz = msd["msd_y"] + msd["msd_z"]
        time = msd["timestep"] * args.dt
        msd = append_fields(msd, "time", time, usemask=False)
        if args.two:
            msd = append_fields(msd, "msd_xy", xy, usemask=False)
            msd = append_fields(msd, "msd_xz", xz, usemask=False)
            msd = append_fields(msd, "msd_yz", yz, usemask=False)

        if args.i:
            figs[(T, n)], axs[(T, n)] = plt.subplots()
            axs[(T, n)].plot(msd["time"], msd["msd_x"], ".", color="#1b9e77", label="X")
            axs[(T, n)].plot(msd["time"], msd["msd_y"], ".", color="#d95f02", label="Y")
            axs[(T, n)].plot(msd["time"], msd["msd_z"], ".", color="#7570b3", label="Z")
            axs[(T, n)].set_xlabel("Time (ps)")
            axs[(T, n)].set_ylabel(r"MSD ($\AA^2$)")
            lgnd = plt.legend(frameon=False, loc="best")
            for handle in lgnd.legendHandles:
                handle._sizes = [30]
            figs[(T, n)].savefig(
                f"{args.base}_T{T}_dir{n}_msd.png", bbox_inches="tight"
            )

        else:
            fig, ax = plt.subplots()
            ax.plot(msd["time"], msd["msd_x"], ".", color="#1b9e77", label="X")
            ax.plot(msd["time"], msd["msd_y"], ".", color="#d95f02", label="Y")
            ax.plot(msd["time"], msd["msd_z"], ".", color="#7570b3", label="Z")
            ax.set_xlabel("Time (ps)")
            ax.set_ylabel(r"MSD ($\AA^2$)")
            # ax.yaxis.set_ticks([0,5,10,15])
            lgnd = plt.legend(frameon=False, loc="best")
            for handle in lgnd.legendHandles:
                handle._sizes = [30]
            fig.savefig(f"{args.base}_T{T}_dir{n}_msd.png", bbox_inches="tight")
            plt.close()

            if args.two:
                fig, ax = plt.subplots()
                ax.plot(msd["time"], msd["msd_xy"], ".", color="#1b9e77", label="XY")
                ax.plot(msd["time"], msd["msd_xz"], ".", color="#d95f02", label="XZ")
                ax.plot(msd["time"], msd["msd_yz"], ".", color="#7570b3", label="YZ")
                ax.set_xlabel("Time (ps)")
                ax.set_ylabel(r"MSD ($\AA^2$)")
                lgnd = plt.legend(frameon=False, loc="best")
                for handle in lgnd.legendHandles:
                    handle._sizes = [30]
                fig.savefig(f"{args.base}_T{T}_dir{n}_2D_msd.png", bbox_inches="tight")
                plt.close()

            if args.three:
                fig, ax = plt.subplots()
                ax.plot(msd["time"], msd["msd_xyz"], ".", color="#1b9e77", label="XYZ")
                ax.set_xlabel("Time (ps)")
                ax.set_ylabel(r"MSD ($\AA^2$)")
                fig.savefig(f"{args.base}_T{T}_dir{n}_3D_msd.png", bbox_inches="tight")
                plt.close()


if args.i:
    cont = True
    while cont:
        T = int(input(f"Enter the temperature to plot ({list(range(900,1401,50))}): "))
        n = int(input(f"Enter the directory to plot({list(range(1,4))}): "))

        figs[(T, n)].show()

        cont = promptForContinue()
