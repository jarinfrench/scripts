#! /usr/bin/env python3

import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np
import pwlf  # PieceWise Linear Fitting

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file', description = "Plot the MSD in 1, 2, and 3D from a LAMMPS simulation")
parser.add_argument('file', help = "The data file containing the MSD data.")
parser.add_argument('-d','--dt', type = float, help = "The timestep used in the simulation. Default = 0.002", default = 0.002)

args = parser.parse_args()

data = np.loadtxt(args.file)
oneDMax = np.max(data[:,1:4])
twoDMax = np.max(np.array([np.add(data[:,1], data[:,2]), np.add(data[:,1], data[:,3]), np.add(data[:,2], data[:,3])]))

fig = plt.figure(figsize=(12,6))
grid = plt.GridSpec(2, 6, wspace = 0.5, hspace = 0.4)
x_plt = fig.add_subplot(grid[0,0])
y_plt = fig.add_subplot(grid[0,1])
z_plt = fig.add_subplot(grid[0,2])
xy_plt = fig.add_subplot(grid[1,0])
xz_plt = fig.add_subplot(grid[1,1])
yz_plt = fig.add_subplot(grid[1,2])
xyz_plt = fig.add_subplot(grid[0:, 3:])

time = [i * args.dt for i in data[:,0]]

x_model = pwlf.PiecewiseLinFit(time, data[:,1]) # slopes are angstroms^2/ps
y_model = pwlf.PiecewiseLinFit(time, data[:,2])
z_model = pwlf.PiecewiseLinFit(time, data[:,3])
xy_model = pwlf.PiecewiseLinFit(time, np.add(data[:,1], data[:,2]))
xz_model = pwlf.PiecewiseLinFit(time, np.add(data[:,1], data[:,3]))
yz_model = pwlf.PiecewiseLinFit(time, np.add(data[:,2], data[:,3]))
xyz_model = pwlf.PiecewiseLinFit(time, data[:,4])

x_model.fit(2)
y_model.fit(2)
z_model.fit(2)
xy_model.fit(2)
xz_model.fit(2)
yz_model.fit(2)
xyz_model.fit(2)
min_break = min([x_model.fit_breaks[1], y_model.fit_breaks[1], z_model.fit_breaks[1],
                 xy_model.fit_breaks[1], xz_model.fit_breaks[1], yz_model.fit_breaks[1],
                 xyz_model.fit_breaks[1]])

# WARNING: This fitting find the ideal breakpoint for _each_ dimension. Thus, there
# is a high chance that the different lines will not have the same break points
# It is not likely, however, that the slope will drastically change if a smaller
# number of points are included.
x_D = x_model.slopes[0] / 2
y_D = y_model.slopes[0] / 2
z_D = z_model.slopes[0] / 2
xy_D = xy_model.slopes[0] / 4
xz_D = xz_model.slopes[0] / 4
yz_D = yz_model.slopes[0] / 4
xyz_D = xyz_model.slopes[0] / 6

fitted_time = np.linspace(0, min_break, 2)
x_plt.scatter(time, data[:,1], 1)
x_plt.plot(fitted_time, x_model.predict(fitted_time), 'r--')

y_plt.scatter(time, data[:,2], 1)
y_plt.plot(fitted_time, y_model.predict(fitted_time), 'r--')

z_plt.scatter(time, data[:,3], 1)
z_plt.plot(fitted_time, z_model.predict(fitted_time), 'r--')

xy_plt.scatter(time, np.add(data[:,1], data[:,2]), 1)
xy_plt.plot(fitted_time, xy_model.predict(fitted_time), 'r--')

xz_plt.scatter(time, np.add(data[:,1], data[:,3]), 1)
xz_plt.plot(fitted_time, xz_model.predict(fitted_time), 'r--')

yz_plt.scatter(time, np.add(data[:,2], data[:,3]), 1)
yz_plt.plot(fitted_time, yz_model.predict(fitted_time), 'r--')

xyz_plt.scatter(time, data[:,4], 1)
xyz_plt.plot(fitted_time, xyz_model.predict(fitted_time), 'r--')
xyz_plt.text(0.05, 0.65, f"$D_x$ = {x_D:.2e} $\AA^2/ps$\n$D_y$ = {y_D:.2e}\n$D_z$ = {z_D:.2e}\n"
             f"$D_{{xy}}$ = {xy_D:.2e}\n$D_{{xz}}$ = {xz_D:.2e}\n$D_{{yz}}$ = {yz_D:.2e}\n"
             f"$D_{{xyz}}$ = {xyz_D:.2e}", ha = 'left', transform = xyz_plt.transAxes)

xyz_plt.text(0.50, 0.05, f"Fit ends at time t = {min_break:.0f} ps", ha = 'left', transform = xyz_plt.transAxes)

with open("diffusion_data.csv", 'w') as f:
    csvfile = csv.writer(f, delimiter=' ', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
    csvfile.writerow(["#Time_Stop", "Dx", "Dy", "Dz", "Dxy", "Dxz", "Dyz", "Dxyz"])
    csvfile.writerow([min_break, x_D, y_D, z_D, xy_D, xz_D, yz_D, xyz_D])

x_plt.set_title("X")
y_plt.set_title("Y")
z_plt.set_title("Z")
xy_plt.set_title("XY")
xz_plt.set_title("XZ")
yz_plt.set_title("YZ")
xyz_plt.set_title("XYZ")

x_plt.set_ylim([0,oneDMax])
y_plt.set_ylim([0,oneDMax])
z_plt.set_ylim([0,oneDMax])
xy_plt.set_ylim([0,twoDMax])
xz_plt.set_ylim([0,twoDMax])
yz_plt.set_ylim([0,twoDMax])
xyz_plt.set_ylim([0,np.max(data[:,4])])

x_plt.set_xlim([0,max(time)])
y_plt.set_xlim([0,max(time)])
z_plt.set_xlim([0,max(time)])
xy_plt.set_xlim([0,max(time)])
xz_plt.set_xlim([0,max(time)])
yz_plt.set_xlim([0,max(time)])
xyz_plt.set_xlim([0,max(time)])

x_plt.xaxis.set_major_locator(plt.MaxNLocator(4))
y_plt.xaxis.set_major_locator(plt.MaxNLocator(4))
z_plt.xaxis.set_major_locator(plt.MaxNLocator(4))
xy_plt.xaxis.set_major_locator(plt.MaxNLocator(4))
xz_plt.xaxis.set_major_locator(plt.MaxNLocator(4))
yz_plt.xaxis.set_major_locator(plt.MaxNLocator(4))
xyz_plt.xaxis.set_major_locator(plt.MaxNLocator(4))

fig.text(0.5, 0.04, 'Time (ps)', ha = 'center', fontsize = 16)
fig.text(0.04, 0.5, r"MSD ($\AA^2$)", va = 'center', rotation = 'vertical', fontsize = 16)

plt.savefig("MSD_diffusion_plot.png", bbox_inches='tight')
# plt.show()
