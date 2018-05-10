from __future__ import division, print_function
from sys import argv, exit
import numpy as np
from myModules import *

if not len(argv) > 1:
    filename = []
    filename.append(raw_input("Please enter the file to be processed.  Note that this should be a file generated from the \"calculate_displacement.cpp\" script: "))
else:
    filename = []
    for i in range(len(argv) - 1):
        filename.append(argv[i + 1])

print("Please enter the coordinates for the crystallographic axes (e.g. [1, 0, 0] for x):")
isScalar = True
x_vec = []
y_vec = []
z_vec = []
while isScalar or (not len(x_vec) == 3):
    x_vec = input("X coordinates: ")
    isScalar = checkScalarInput(x_vec)
while isScalar or (not len(y_vec) == 3):
    y_vec = input("Y coordinates: ")
    isScalar = checkScalarInput(y_vec)
while isScalar or (not len(z_vec) == 3):
    z_vec = input("Z coordinates: ")
    isScalar = checkScalarInput(z_vec)

# verify that the coordinates are orthogonal
if not checkOrthogonality(x_vec,y_vec,z_vec):
    print("Crystal axes are not orthogonal!")
    exit(1)

x_norm = np.linalg.norm(x_vec)
y_norm = np.linalg.norm(y_vec)
z_norm = np.linalg.norm(z_vec)

rotation_matrix = np.array([x_vec/x_norm, y_vec/y_norm, z_vec/z_norm])
rotation_matrix_inv = np.transpose(rotation_matrix)
converted_jump_vector = []

for i in range(len(filename)):
    print("Processing file \"%s\""%filename[i])
    with open(filename[i], 'r') as f:
        # the first line contains the variable names
        f.readline();
        for line in f:
            data = line.split()
            jump_vector = np.array([float(data[7]), float(data[8]), float(data[9])])
            converted_jump_vector.append(rotation_matrix_inv.dot(normalize(jump_vector)))
