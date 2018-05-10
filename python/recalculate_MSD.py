#! /usr/bin/env python
from __future__ import division, print_function
from itertools import izip

def calculateCom(filename):
    com_U = [0,0,0]
    com_O = [0,0,0]
    com_total = [0,0,0]
    with open(filename, "r") as fopen:
        junk = fopen.readline()
        junk = fopen.readline()
        junk = fopen.readline()
        n_atoms = int(fopen.readline())
        junk = fopen.readline()
        junk = fopen.readline()
        junk = fopen.readline()
        junk = fopen.readline()
        junk = fopen.readline()
        n_U = n_atoms / 3
        n_O = n_atoms - n_U
        mass_U = 238.0289
        mass_O = 15.9994
        total_mass_U = mass_U * n_U
        total_mass_O = mass_O * n_O
        total_mass = total_mass_U + total_mass_O
        for line in fopen:
            atom_type = int(line.split()[1])
            x_data = float(line.split()[6])
            y_data = float(line.split()[7])
            z_data = float(line.split()[8])

            if atom_type == 1:
                com_U =     [i + mass_U * j for i,j in izip(com_U,[x_data,y_data,z_data])]
                com_total = [i + mass_U * j for i,j in izip(com_U,[x_data,y_data,z_data])]
            elif atom_type == 2:
                com_O = [i + mass_O * j for i,j in izip(com_O,[x_data,y_data,z_data])]
                com_total = [i + mass_O * j for i,j in izip(com_O,[x_data,y_data,z_data])]
        for i in range(len(com_U)):
            com_U[i] = com_U[i] / total_mass_U
            com_O[i] = com_O[i] / total_mass_O
            com_total[i] = com_total[i] / total_mass
    return com_U, com_O, com_total




file_O = "A_MSD_O.dat"
file_total = "A_MSD.dat"
file_U = "A_MSD_U_approx.dat"
f_U = open(file_U, "w")
#com_U_original, com_O_original, com_total_original = calculateCom("0.dump")
f_U.write("# Time-averaged data for fix my_ave_o\n")
f_U.write("# TimeStep c_msd_O[1] c_msd_O[2] c_msd_O[3] c_msd_O[4]\n")
f_U.write("# This is an approximate calculation based on the values from the O and total MSD\n")
f_U.write("VARIABLES=\"Timestep\" \"MSD_x\" \"MSD_y\" \"MSD_z\" \"MSD_total\"\n")

with open(file_O) as f_O, open(file_total) as f_total:
    for line_O, line_total in izip(f_O, f_total):

        if line_O.startswith("#") and line_total.startswith("#"):
            continue
        if int(line_O.split()[0]) != int(line_total.split()[0]):
            print("Data does not match!")
            print(repr(line_O.split()[0])+" != "+repr(line_total.split()[0]))
        else:
            data_O = line_O.split()
            data_total = line_total.split()
            data_U = []
            for i,j in izip(data_O, data_total):
                data_U.append((float(j)*12 - float(i)*8)/4)
            f_U.write("%d %8.6f %8.6f %8.6f %8.6f\n"%(data_U[0], data_U[1], data_U[2], data_U[3],data_U[4]))
