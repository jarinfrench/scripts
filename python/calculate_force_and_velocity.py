#! /usr/bin/env python3

import argparse
import itertools
import math
import os
from sys import exit

import matplotlib.pyplot as plt
import numpy as np


# Creates an iterator over the previous item (a), current item (b), and the next item
# (c) as a tuple, for each item that has all three.
def threes(iterator):
    "s ->(s0,s1,s2),(s1,s2,s3),(s2,s3,s4), ..."
    a, b, c = itertools.tee(iterator, 3)
    next(b, None)
    next(c, None)
    next(c, None)
    return zip(a, b, c)


def calculateLatticeParam(T, c, potential=0):
    dbFile = "/home/jarinf/projects/scripts/lattice_params.db"
    data = []
    with open(dbFile, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            else:
                data.append(line)
    if potential == 0:
        print(f"There are {len(data)} fits:")
        for i in range(len(data)):
            print(f"  {i + 1} - {data[i].split()[0]}")
        potential = int(input("Please specify the fit to use: "))

    while potential > len(data) or potential < 1:
        print(
            "Invalid number. Additional potentials should be added to the database file."
        )
        potential = int(input("Please specify the fit to use: "))

    name, c_min, c_max, T_min, T_max, c3, c2, c1, T3, T2, T1, cT, cT2, c2T, a0 = [
        i if idx == 0 else float(i) for idx, i in enumerate(data[potential - 1].split())
    ]
    print(f"Using the {name} potential")

    if T < float(T_min) or T > float(T_max) or c < float(c_min) or c > float(c_max):
        print(
            f"Temperature or concentration out of fitted range ({T_min} <= T <= {T_max}, {c_min} <= c <= {c_max})"
        )
        exit(3)
    else:
        return (
            c**3 * c3
            + c**2 * c2
            + c * c1
            + T**3 * T3
            + T**2
            + T2
            + T * T1
            + cT * c * T
            + cT2 * c * T**2
            + c**2 * T * c2T
            + a0
        )


def onClick(event):
    global ix, iy
    if event.inaxes == ax.axes:
        ix, iy = (event.xdata, event.ydata)
        fig.canvas.mpl_disconnect(cid)
        plt.close()
    return


def onButtonPress(event):
    global use1
    if event.key == "t":
        use1 = not use1
        plotN(t, n1, n2, use1, fig)
        plt.show()


def onClose(event):
    global ix, iy
    if ix is None:
        exit(0)


def plotN(t, n1, n2, use1, fig):
    plt.clf()
    global ax
    ax = fig.add_subplot(111)
    if use1:
        ax.plot(t, n1, "ro")
    else:
        ax.plot(t, n2, "ro")
    plt.text(0.45, 0.6, "Press 't' to toggle between n1 and n2", transform=ax.transAxes)
    plt.text(0.45, 0.55, "Click where the growth stops", transform=ax.transAxes)


def calculateR(N, a0, Lz):
    return math.sqrt(N * a0**3 / (4 * math.pi * Lz))


def calculateForce(r, gamma):
    if r == 0:
        return 0
    else:
        return gamma / (r * 0.1)  # The 0.1 converts the value to GPa


def calculateVelocity(ns, ts, a0, Lz):
    # area is given as pi*r^2, and N*a0^3/(4*Lz), and we are calculating r for a set of
    # three points.
    coeff = math.sqrt(a0**3 / (4 * math.pi * Lz))  # fixed for a specific system
    if not len(ns) == 3 or not len(ts) == 3:
        print("Error calculating velocity of grain boundary")
    sqrtN = [math.sqrt(n) for n in ns]  # The number of atoms in the grain
    rs = [coeff * i for i in sqrtN]
    fit = np.polyfit(ts, rs, 1)
    return fit[0] * 100  # The 100 converts the value to m/s


parser = argparse.ArgumentParser(
    description="Calculates the velocity and forces for an assumed cylindrical grain boundary given a data file in the format <timestep> <n grain 1> <n grain 2>"
)
parser.add_argument("t", metavar="T", type=float, help="Temperature of the simulation")
parser.add_argument("c", type=float, help="Concentration of solute in the structure")
parser.add_argument("l", metavar="Lz", type=float, help="Thickness of the grain")
parser.add_argument("a", metavar="a0", type=float, help="Lattice parameter at 0 K")
parser.add_argument(
    "gamma",
    type=float,
    help="Random grain boundary energy value (use 1 for the 'reduced' force)",
)
parser.add_argument(
    "-p",
    "--potential",
    type=int,
    help="Number of the potential to use from the database file",
    default=0,
)
parser.add_argument(
    "-g",
    "--graph",
    action="store_true",
    help="Option to display a graph showing the grain growth (N vs t) for help in determining when growth stops",
)
parser.add_argument(
    "-u",
    "--use",
    type=int,
    choices=[1, 2],
    help="If the data set to be used in calculating the data (between n1 and n2) is known, specify with this option",
)

args = parser.parse_args()

# May change this to be an argument
dataFile = (
    "data.txt"  # data file containing the timestep, and number of atoms in each grain
)
dataOutfile = "force_velocity_data.txt"
a0 = calculateLatticeParam(args.t, args.c, args.potential)

data = []
last_point = -1
current_point = 0
if os.path.exists("slope_calc.txt"):
    # If this file exists, use it's last fitted point to end the dataFile read early
    with open("slope_calc.txt", "r") as f:
        # format: area_data.txt dA/dt = <> error = <> r_sq = <> fit to points <>:<>
        for line in f:
            pass
        last_point = int(line.split(":")[-1]) - 1
with open(dataFile, "r") as f:
    for _ in range(3):
        next(f)
    for line in f:
        if last_point == current_point:
            break
        data.append(line)
        current_point += 1

vel = []
force = []
r = []
t = []
n1 = []
n2 = []
ts = []
n1s = []
n2s = []

# separate out the data in the ways we will use it
for n, i in enumerate(threes(data)):
    ts.append(
        [
            int(i[0].split()[0]) * 0.002,
            int(i[1].split()[0]) * 0.002,
            int(i[2].split()[0]) * 0.002,
        ]
    )
    n1s.append([int(i[0].split()[1]), int(i[1].split()[1]), int(i[2].split()[1])])
    n2s.append([int(i[0].split()[2]), int(i[1].split()[2]), int(i[2].split()[2])])
    t.append(int(i[1].split()[0]) * 0.002)
    n1.append(int(i[1].split()[1]))
    n2.append(int(i[1].split()[2]))

if args.use is None or args.use == 1:
    use1 = True
else:
    use1 = False

if args.graph:
    fig = plt.figure()
    plotN(t, n1, n2, use1, fig)
    cid = fig.canvas.mpl_connect("button_press_event", onClick)
    cid2 = fig.canvas.mpl_connect("key_press_event", onButtonPress)
    fig.canvas.mpl_connect("close_event", onClose)
    plt.show()
    idx = (np.abs(t - ix)).argmin()
    tstop = t[idx]

for i in range(len(n1)):
    if use1:
        r.append(calculateR(n1[i], args.a, args.l))
        vel.append(calculateVelocity(n1s[i], ts[i], args.a, args.l))
    else:
        r.append(calculateR(n2[i], args.a, args.l))
        vel.append(calculateVelocity(n2s[i], ts[i], args.a, args.l))

    force.append(calculateForce(r[i], args.gamma))

written = False
with open(dataOutfile, "w") as f:
    f.write("# Time, n_grain, r, force, and instantaneous velocity of the boundary.\n")
    f.write("# Note that a negative velocity indicates grain shrinking.\n")
    for i in range(len(force)):
        if args.graph and t[i] > tstop and not written:
            f.write("\n\n# Note that after this point, growth appears to stop\n")
            written = True
        f.write(
            "{time} {num} {rad} {force} {vel}\n".format(
                time=t[i],
                num=(n1[i] if use1 else n2[i]),
                rad=r[i],
                force=force[i],
                vel=vel[i],
            )
        )
