#! /usr/bin/env python3

from sys import argv, exit
import argparse, os
from myModules import *

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file [file file file ...] -s scale (-s scale2 -s scale3) [options]', description = "Expand a lattice by a constant.")
parser.add_argument('file', nargs = '+', help = "The file(s) to process")
parser.add_argument('-s','--scale', action = 'append', type = float, required = True, help = "Scaling factor(s)")
parser.add_argument('-o','--output', default = "new_structure", help = "Output filename base")

args = parser.parse_args()

uniform = False # one scaling factor
xy_only = False # two scaling factors, scaling the xy-direction with the first one, and the z direction with the second
scale_all = False # scale all directions independently

if len(args.scale) == 1:
    uniform = True
    if args.scale[0] == 1.0:
        print("No lattice scaling needed. Exiting...")
        exit(0)
elif len(args.scale) == 2:
    xy_only = True
elif len(args.scale) == 3:
    scale_all = True
else:
    print("Too many scaling factors passed in. Only using the first three.")
    scale_all = True
    args.scale = args.scale[0:3]

for file in args.file:
    if args.output.find(".dat") != -1: # if the .dat extension was found
        output_filename = args.output # the exact filename will be used
    else: # otherwise, we use the basename and check for already existing files.
        output_filename = "{base}.dat".format(base = args.output)
        output_filename = verify_new_file(output_filename)

    if not os.path.splitext(file)[1] == ".dat":
        print("Incorrect file type {}".format(os.path.splitext(file)[1]))
        exit(1)

    fout = open(output_filename, 'w')
    with open(file) as f:
        # read the header data
        lower_bounds = {"x":0.0, "y":0.0, "z":0.0, "xu":0.0, "yu":0.0, "zu":0.0}
        comment = f.readline()
        data_format = comment.split(":")[1] # we end up with a string that looks like: " [ID type x y z]\n"
        if uniform:
            fout.write("The {file} structure has been expanded by a factor of {scale}:{data_str}\n".format(file = file, scale = args.scale[0], data_str = data_format))
        elif xy_only:
            fout.write("The {file} structure has been expanded by a factor of x,y * {scale1}, z * {scale2}:{data_str}\n".format(file = file, scale1 = args.scale[0], scale2 = args.scale[1], data_str = data_format))
        elif scale_all:
            fout.write("The {file} structure has been expanded by a factor of x * {scale1}, y * {scale2}, z * {scale3}:{data_str}\n".format(file = file, scale1 = args.scale[0], scale2 = args.scale[1], scale3 = args.scale[2], data_str = data_format))
        else:
            print("Error: unknown scaling relationship.")
            exit(1)
        data_format = data_format.strip(" []\n").split()
        tmp = f.readline() # blank line
        n_atoms = f.readline() # number of atoms
        fout.write("{}".format(n_atoms))
        n_types = f.readline() # number of atom types
        fout.write("{}".format(n_types))
        tmp = f.readline() # x bounds
        xlo = float(tmp.split()[0]) # lower bound
        xhi = float(tmp.split()[1]) # upper bound
        lower_bounds["x"] = xlo
        lower_bounds["xu"] = xlo
        fout.write("{:8.6f} {:8.6f} xlo xhi\n".format(xlo, (xhi-xlo) * args.scale[0] + xlo)) # the x direction is always scaled by args.scale[0]
        tmp = f.readline() # y bounds
        ylo = float(tmp.split()[0]) # lower bound
        yhi = float(tmp.split()[1]) # upper bound
        lower_bounds["y"] = ylo
        lower_bounds["yu"] = ylo
        if uniform or xy_only: # the y direction is scaled by args.scale[0] in the uniform case, and in the x-y case
            fout.write("{:8.6f} {:8.6f} ylo yhi\n".format(ylo, (yhi-ylo) * args.scale[0] + ylo))
        elif scale_all: # otherwise it is scaled by args.scale[1]
            fout.write("{:8.6f} {:8.6f} ylo yhi\n".format(ylo, (yhi-ylo) * args.scale[1] + ylo))
        else:
            print("Error: unknown scaling relationship.")
            os.exit(1)
        tmp = f.readline() # z bounds
        zlo = float(tmp.split()[0]) # lower bound
        zhi = float(tmp.split()[1]) # upper bound
        lower_bounds["z"] = zlo
        lower_bounds["zu"] = zlo
        if uniform: # z scales with args.scale[0] if one value is given
            fout.write("{:8.6f} {:8.6f} zlo zhi\n\n".format(zlo, (zhi-zlo) * args.scale[0] + zlo))
        elif xy_only: # z scales with args.scale[1] if two values are given
            fout.write("{:8.6f} {:8.6f} zlo zhi\n\n".format(zlo, (zhi-zlo) * args.scale[1] + zlo))
        elif scale_all: # z scales with args.scale[2] if three (or more) values are given
            fout.write("{:8.6f} {:8.6f} zlo zhi\n\n".format(zlo, (zhi-zlo) * args.scale[2] + zlo))
        tmp = f.readline() # blank line
        tmp = f.readline() # "Atoms"
        fout.write("{}\n".format(tmp))
        tmp = f.readline() # blank line

        for line in f:
            data = line.split()
            for i in range(len(data)):
                # This strips off the [,], and ' ' at the ends, and divides it into it's component parts
                if data_format[i] in ["x", "xu"]:
                    data[i] = "{:8.6f}".format((float(data[i]) - lower_bounds[data_format[i]]) * args.scale[0] + lower_bounds[data_format[i]])
                elif data_format[i] in ["y", "yu"]:
                    if uniform or xy_only:
                        data[i] = "{:8.6f}".format((float(data[i]) - lower_bounds[data_format[i]]) * args.scale[0] + lower_bounds[data_format[i]])
                    elif scale_all:
                        data[i] = "{:8.6f}".format((float(data[i]) - lower_bounds[data_format[i]]) * args.scale[1] + lower_bounds[data_format[i]])
                    else:
                        print("Error: unknown scaling relationship.")
                        exit(1)
                elif data_format[i] in ["z", "zu"]:
                    if uniform:
                        data[i] = "{:8.6f}".format((float(data[i]) - lower_bounds[data_format[i]]) * args.scale[0] + lower_bounds[data_format[i]])
                    elif xy_only:
                        data[i] = "{:8.6f}".format((float(data[i]) - lower_bounds[data_format[i]]) * args.scale[1] + lower_bounds[data_format[i]])
                    elif scale_all:
                        data[i] = "{:8.6f}".format((float(data[i]) - lower_bounds[data_format[i]]) * args.scale[2] + lower_bounds[data_format[i]])
                    else:
                        print("Error: unknown scaling relationship.")
                        exit(1)
            update = " ".join(data)
            fout.write("{}\n".format(update))
