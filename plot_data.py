#! /usr/bin/python

from __future__ import division, print_function
import matplotlib.pyplot as plt
from sys import argv

# Simple python script to create a plot of two-dimensional data.

#if len(argv) != 4:
#    filename = input("Please enter the filename containing the data: ")
#    x_label = input("Please enter the x-axis label: ")
#    y_label = input("Please enter the y-axis label: ")
#else:
#    filename = argv[1]
#    x_label = argv[2]
#    y_label = argv[3]

fin = open("100Tilt_individual_energy_rcut2.txt")
fin2 = open("100Tilt_individual_energy_rcut2.5.txt")
fin3 = open("100Tilt_individual_energy_rcut3.txt")
x_data = [0.0]
x2 = [0.0]
x3 = [0.0]
y_data = [0.0]
y2 = [0.0]
y3 = [0.0]

while True:
    data = fin.readline().split()
    data2 = fin2.readline().split()
    data3 = fin3.readline().split()
    if not (data or data2 or data3):
        break
    elif len(data) > 2:
        print("Too many values in line")
        continue
    else:
        x_data.append(float(data[0]))
        y_data.append(float(data[1]))
        x2.append(float(data2[0]))
        y2.append(float(data2[1]))
        x3.append(float(data3[0]))
        y3.append(float(data3[1]))

# Sort the data
xy1 = zip(x_data, y_data)
xy2 = zip(x2, y2)
xy3 = zip(x3, y3)
xy1.sort()
xy2.sort()
xy3.sort()
x_data = [x for x,y in xy1]
y_data = [y for x,y in xy1]
x2 = [x for x,y in xy2]
y2 = [y for x,y in xy2]
x3 = [x for x,y in xy3]
y3 = [y for x,y in xy3]

# calculate the residuals
y_data_res = [0.0 for i in range(len(y_data))]
y2_res = [0.0 for i in range(len(y2))]
y3_res = [0.0 for i in range(len(y3))]
for i in range(len(y_data)):
    y_avg = (y_data[i] + y2[i] + y3[i]) / 3
    y_data_res[i] = abs(y_data[i] - y_avg)
    y2_res[i] = abs(y2[i] - y_avg)
    y3_res[i] = abs(y3[i] - y_avg)

print("Average residual for r_cut = 2.0:", sum(y_data_res)/len(y_data_res))
print("Average residual for r_cut = 2.5:", sum(y2_res)/len(y2_res))
print("Average residual for r_cut = 3.0:", sum(y3_res)/len(y3_res))

# Plot the results
plt.plot(x_data, y_data, 'bo-', label="r_cut = 2.0")
plt.plot(x2, y2, 'r^-', label="r_cut = 2.5")
plt.plot(x3, y3, 'gs-', label="r_cut = 3.0")
plt.xlabel("Angle (degrees)")
plt.ylabel(r"Energy (J/m^2)")
plt.xlim(0, max(x_data)*1.1)
plt.legend(loc='best')
plt.title("Grain Boundary Energy")
plt.figure(2)
plt.plot(x_data, y_data_res, 'bo-', label="r_cut = 2.0 residuals")
plt.plot(x_data, y2_res, 'r^-', label="r_cut = 2.5 residuals")
plt.plot(x_data, y3_res, 'gs-', label="r_cut = 3.0 residuals")
plt.plot(x_data, [0 for i in range(len(x_data))], 'k-')
#plt.xlabel(x_label)
#plt.ylabel(y_label)
plt.xlabel("Angle (degrees)")
plt.ylabel(r"Energy (J/m^2)")
plt.xlim(0, max(x_data)*1.1)
#plt.axis((0, max(x_data), 0, 1.1*max(y_data)))
plt.legend(loc='best')
plt.title("Grain Boundary Energy Residuals")

# Plot the individual data sets
plt.figure(3)
plt.plot(x_data, y_data, 'bo-')
plt.xlabel("Angle (degrees)")
plt.ylabel(r"Energy (J/m^2)")
plt.xlim(0, max(x_data)*1.1)

plt.figure(4)
plt.plot(x2, y2, 'r^-')
plt.xlabel("Angle (degrees)")
plt.ylabel(r"Energy (J/m^2)")
plt.xlim(0, max(x2)*1.1)

plt.figure(5)
plt.plot(x3, y3, 'gs-')
plt.xlabel("Angle (degrees)")
plt.ylabel(r"Energy (J/m^2)")
plt.xlim(0, max(x3)*1.1)
plt.show()
