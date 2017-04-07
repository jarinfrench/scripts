#! /usr/bin/python

from __future__ import division, print_function
import matplotlib.pyplot as plt
from sys import argv

# Simple python script to create a plot of two-dimensional data.

if len(argv) != 2:
    filename = input("Please enter the filename containing the data: ")
    #x_label = input("Please enter the x-axis label: ")
    #y_label = input("Please enter the y-axis label: ")
else:
    filename = argv[1]
    #x_label = argv[2]
    #y_label = argv[3]

x_data = [0.0]
y_data = [0.0]

fin = open(filename)

while True:
    data = fin.readline().split()
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

# Plot the results
plt.plot(x_data, y_data, 'bo-')
plt.xlabel("Angle (degrees)")
plt.ylabel(r"Energy (J/$m^2$)")
plt.xlim(0, max(x_data)*1.1)
plt.title("Grain Boundary Energy")
plt.show()
