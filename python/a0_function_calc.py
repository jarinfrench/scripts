#! /usr/bin/env python3

import argparse
import csv
import statistics as stats
from inspect import cleandoc
from sys import exit

import matplotlib.pyplot as plt
import myModules as mine
import numpy as np
import scipy.optimize as optimize

# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def func3D(data, z0, x1, y1, x2, y2, x3, y3, xy, xy2, x2y):
    x = data[0]  # fraction
    y = data[1]  # temperature
    return (
        x3 * x**3
        + x2 * x**2
        + x1 * x
        + z0
        + y3 * y**3
        + y2 * y**2
        + y1 * y
        + xy * x * y
        + xy2 * x * y * y
        + x2y * x * x * y
    )


def func2D(data, y0, m1, m2, m3):
    x = np.array(data)  # temperature
    return m3 * x**3 + m2 * x**2 + m1 * x + y0


def SurfacePlot(func, data, fittedParameters):
    def fracVsTemp(data, y0, x1, x2):
        x = np.array(data)
        return x2 * x**2 + x1 * x + y0

    plt.figure()
    axes = plt.axes(projection="3d")

    x_data = data[:, 0]  # fraction
    y_data = data[:, 1]  # temperature
    z_data = data[:, 2]  # lattice parameter

    # determine the bounds of the surface
    # xmin = min(x_data)
    # xmax = max(x_data)
    # ymin = min(y_data)
    # ymax = max(y_data)

    # guess = (np.mean(y_data), 0, 1)  # intercept, linear slope, parabolic slope
    # I need to figure out a way to generate the boundary conditions of the surface...
    xModel = np.linspace(min(x_data), max(x_data), 50)
    yModel = np.linspace(min(y_data), max(y_data), 50)
    X, Y = np.meshgrid(xModel, yModel)

    Z = func(np.array([X, Y]), *fittedParameters)
    # R = Z
    # R = np.where(
    #     Y + 250000 * X * X + 13886.3636363636 * X < 1872.5, Z, np.nan
    # )  # for U_Mo_Xe_ternary_eam - UXe
    # R = np.where(
    #     Y - 5994.006 * X * X + 2261.072 * X <= 2011.888, Z, np.nan
    # )  # for U_Mo_Xe_ternary_eam - UMo
    # R = np.where(Y - 4000.0 / 3.0 * X <= 1600.0, Z, np.nan) # for U_Mo_adp
    # R = np.where((X > 0.03) & (Y < 1600), Z, np.nan)
    # R = np.where((X > 0.12) & (Y < 1700), R, np.nan)
    # R = np.where((X > 0.15) & (Y < 1800), R, np.nan)
    # R = np.where((X > 0.24) & (Y < 1900), R, np.nan)
    # R = np.where((X > 0.27) & (Y < 2000), R, np.nan)

    # NOTE: set vmin, vmax according to the min/max of the lattice parameter, giving a
    # bit of space
    # Then, set_zlim is just 0.01 lower/higher (respectively)
    axes.plot_surface(
        X,
        Y,
        Z,
        alpha=0.75,
        cmap=cm.coolwarm,
        linewidth=0,
        antialiased=False,
        vmin=min(z_data) * 0.99,
        vmax=max(z_data) * 1.01,
        zorder=2,
    )
    axes.set_zlim(min(z_data) * 0.99 - 0.01, max(z_data) * 1.01 + 0.01)
    axes.scatter(x_data, y_data, z_data, c="black", alpha=1, zorder=1)
    axes.set_xlabel("Impurity atomic fraction")  # change as needed.
    axes.set_ylabel("Temperature")
    axes.set_zlabel("Lattice Parameter")

    plt.show()
    plt.close("all")


parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] file a0 runs(3)",
    description="Perform a fit to the lattice parameter as a function of temperature"
    + "and concentration",
)
parser.add_argument(
    "file", help="The file containing the temperature and lattice parameters"
)
parser.add_argument("a0", type=float, help="The estimated lattice_parameter")
parser.add_argument(
    "runs",
    type=int,
    nargs="?",
    default=3,
    help="The number of runs per temperature. Default: 3",
)
parser.add_argument(
    "-T",
    "--temperature-only",
    action="store_true",
    help="Calculate a fit using temperature only",
)
parser.add_argument(
    "-d",
    "--delimiter",
    default=" ",
    help="Specify the delimiter for reading the data file. Default: ' '",
)

args = parser.parse_args()

data = []

with open(args.file) as f:
    reader = csv.reader(
        f, delimiter=args.delimiter, quotechar='"'
    )  # change the delimiter as needed
    for row in reader:
        if len(row) < 5:
            print(f"Only found {len(row)} elements: did you indicate the delimiter?")
            exit(1)
        data.append(tuple(row))

if args.temperature_only:
    # data format: Temp, ax, ay, az, solid(1) or liquid(0)
    data = np.array(
        [
            (float(i[0]), float(i[1]), float(i[2]), float(i[3]))
            for i in data
            if int(i[4]) == 1
        ]
    )
else:
    # data format: fraction, Temp, ax, ay, az, solid(1) or liquid(0)
    data = np.array(
        [
            (float(i[0]), float(i[1]), float(i[2]), float(i[3]), float(i[4]))
            for i in data
            if int(i[5]) == 1
        ]
    )

max_T = np.max(np.max(data, 1), 0)
data = data[np.argsort(data[:, 0])]  # sort by the first column
avg_data = []

for chunk in mine.chunker(data, args.runs):
    d = np.mean(chunk, axis=0)
    if args.temperature_only:
        avg_data.append(
            np.append(
                d[0:4],
                [
                    np.mean(chunk[:, 1:4]),
                    stats.stdev(mine.flatten(chunk[:, 1:4])) / np.sqrt(args.runs * 3),
                ],
            )
        )
    else:
        avg_data.append(
            np.append(
                d[0:5],
                [
                    np.mean(chunk[:, 2:5]),
                    stats.stdev(mine.flatten(chunk[:, 2:5])) / np.sqrt(args.runs * 3),
                ],
            )
        )

print(f"Maximum temperature is {max_T}")

if args.temperature_only:
    x_data = np.array([i[0] for i in avg_data])  # temperature
    y_data = np.array([i[4] for i in avg_data])  # average lattice parameter
    err = np.array(
        [i[5] if not i[5] == 0 else 1 for i in avg_data]
    )  # standard deviation / sqrt(n)
    guess = (
        args.a0,
        0,
        1,
        1,
    )  # lattice_parameter, linear slope, parabolic slope, cubic slope
    params, pcov = optimize.curve_fit(func2D, x_data, y_data, guess, err)
    # gives in order: 'z' intercept, linear coefficient,
    # parabolic, and cubic coefficient
    print("z0 = {}\ny1 = {}\ny2 = {}\ny3 = {}".format(*params))
    x_plt = [i for i in np.arange(min(x_data), max(x_data))]
    y_plt = func2D(x_plt, *params)
    plt.plot(x_data, y_data, "r.", x_plt, y_plt, "k-", markersize=8)
    plt.show()
else:
    x_data = np.array([i[0] for i in avg_data])  # fraction
    y_data = np.array([i[1] for i in avg_data])  # temperature
    z_data = np.array([i[5] for i in avg_data])  # average lattice parameter
    err = np.array(
        [i[6] if not i[6] == 0 else 1 for i in avg_data]
    )  # standard deviation / sqrt(n)
    # lattice_parameter, linear fraction slope, linear temperature slope,
    # parabolic fraction slope, parabolic temperature slope
    # cubic fraction slope, cubic temperature slope, fraction*temperature slope
    # fraction*temperature^2 slope, fraction^2*temperature slope
    guess = (args.a0, 0, 1, 1, 1, 1, 1, 1, 1, 1)
    params, pcov = optimize.curve_fit(func3D, [x_data, y_data], z_data, guess, err)
    print(
        cleandoc(
            """z0 = {}
        x1 = {}
        y1 = {}
        x2 = {}
        y2 = {}
        x3 = {}
        y3 = {}
        xy = {}
        xy2 = {}
        x2y = {}""".format(
                *params
            )
        )
    )
    SurfacePlot(
        func3D, np.array([[i, j, k] for i, j, k in zip(x_data, y_data, z_data)]), params
    )
