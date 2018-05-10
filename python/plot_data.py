#! /usr/bin/python

from __future__ import division, print_function
import matplotlib.pyplot as plt
from sys import argv
from os.path import basename, splitext
import csv

# Simple python script to create a plot of two-dimensional data.
filename = []
if len(argv) < 2:
    filename = input("Please enter the filename containing the data: ")

else:
    for i in range(len(argv) - 1):
        filename.append(argv[i+1])

x_label = raw_input("Please enter the x-axis label: ")
y_label = raw_input("Please enter the y-axis label: ")
plot_title = raw_input("Please enter the title of the plot: ")

plot_style = ['bo-', 'go-', 'ro-', 'co-', 'mo-', 'yo-', 'ko-',
              'bv-', 'gv-', 'rv-', 'cv-', 'mv-', 'yv-', 'kv-',
              'b^-', 'g^-', 'r^-', 'c^-', 'm^-', 'y^-', 'k^-',
              'b<-', 'g<-', 'r<-', 'c<-', 'm<-', 'y<-', 'k<-',
              'b>-', 'g>-', 'r>-', 'c>-', 'm>-', 'y>-', 'k>-',
              'b1-', 'g1-', 'r1-', 'c1-', 'm1-', 'y1-', 'k1-',
              'b2-', 'g2-', 'r2-', 'c2-', 'm2-', 'y2-', 'k2-',
              'b3-', 'g3-', 'r3-', 'c3-', 'm3-', 'y3-', 'k3-',
              'b4-', 'g4-', 'r4-', 'c4-', 'm4-', 'y4-', 'k4-',
              'b8-', 'g8-', 'r8-', 'c8-', 'm8-', 'y8-', 'k8-',
              'bs-', 'gs-', 'rs-', 'cs-', 'ms-', 'ys-', 'ks-',
              'bp-', 'gp-', 'rp-', 'cp-', 'mp-', 'yp-', 'kp-',
              'bP-', 'gP-', 'rP-', 'cP-', 'mP-', 'yP-', 'kP-',
              'b*-', 'g*-', 'r*-', 'c*-', 'm*-', 'y*-', 'k*-',
              'bh-', 'gh-', 'rh-', 'ch-', 'mh-', 'yh-', 'kh-',
              'bH-', 'gH-', 'rH-', 'cH-', 'mH-', 'yH-', 'kH-',
              'b+-', 'g+-', 'r+-', 'c+-', 'm+-', 'y+-', 'k+-',
              'bx-', 'gx-', 'rx-', 'cx-', 'mx-', 'yx-', 'kx-',
              'bX-', 'gX-', 'rX-', 'cX-', 'mX-', 'yX-', 'kX-',
              'bD-', 'gD-', 'rD-', 'cD-', 'mD-', 'yD-', 'kD-',
              'bd-', 'gd-', 'rd-', 'cd-', 'md-', 'yd-', 'kd-',
              'b|-', 'g|-', 'r|-', 'c|-', 'm|-', 'y|-', 'k|-',
              'b_-', 'g_-', 'r_-', 'c_-', 'm_-', 'y_-', 'k_-']
x_max = 0.0
x_min = 0.0
for i in range(len(filename)):
    x_data = [0.0]
    y_data = [0.0]
    fin = open(filename[i])
    reader = csv.reader(fin)

    for data in reader:
        if not (data):
            break
        elif len(data) > 2:
            print("Too many values in line")
            continue
        else:
            x_data.append(float(data[0]))
            y_data.append(float(data[1]))

    # Sort the data
    xy1 = zip(x_data, y_data)
    xy1.sort()
    x_data = [x for x,y in xy1]
    y_data = [y for x,y in xy1]

    x_max = max(max(x_data), x_max)
    x_min = min(min(x_data), x_min)
    # Plot the results
    plt.plot(x_data, y_data, plot_style[i], label=splitext(basename(filename[i]))[0])

plt.xlabel(x_label)
plt.ylabel(y_label)
#plt.ylabel(r"Energy (J/m$^2$)")
plt.xlim(x_min, x_max)
plt.title(plot_title)
plt.legend(loc='best')
plt.show()
