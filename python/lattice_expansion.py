#! /usr/bin/env python3

from sys import argv, exit
import argparse, os
from myModules import *

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file [file file file ...] -s scale [options]', description = "Expand a lattice by a constant.")
parser.add_argument('file', nargs = '+', help = "The file(s) to process")
parser.add_argument('-s','--scale', type = float, required = True, help = "Scaling factor")
parser.add_argument('-o','--output', default = "new_structure", help = "Output filename base")

args = parser.parse_args()

for file in args.file:
    output_filename = "{base}.dat".format(base = args.output)
    output_filename = verify_new_file(output_filename)

    if not os.path.splitext(file)[1] == ".dat":
        print("Incorrect file type {}".format(os.path.splitext(file)[1]))
        os.exit(1)

    fout = open(output_filename, 'w')
    with open(file) as f:
        # read the header data
        lower_bounds = {"x":0.0, "y":0.0, "z":0.0, "xu":0.0, "yu":0.0, "zu":0.0}
        comment = f.readline()
        data_format = comment.split(":")[1] # we end up with a string that looks like: " [ID type x y z]\n"
        fout.write("The {file} structure has been expanded by a factor of {scale}:{data_str}\n".format(file = file, scale = args.scale, data_str = data_format))
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
        fout.write("{} {} xlo xhi\n".format(xlo, (xhi-xlo) * args.scale + xlo))
        tmp = f.readline() # y bounds
        ylo = float(tmp.split()[0]) # lower bound
        yhi = float(tmp.split()[1]) # upper bound
        lower_bounds["y"] = ylo
        lower_bounds["yu"] = ylo
        fout.write("{} {} ylo yhi\n".format(ylo, (yhi-ylo) * args.scale + ylo))
        tmp = f.readline() # z bounds
        zlo = float(tmp.split()[0]) # lower bound
        zhi = float(tmp.split()[1]) # upper bound
        lower_bounds["z"] = zlo
        lower_bounds["zu"] = zlo
        fout.write("{} {} zlo zhi\n\n".format(zlo, (zhi-zlo) * args.scale + zlo))
        tmp = f.readline() # blank line
        tmp = f.readline() # "Atoms"
        fout.write("{}\n".format(tmp))
        tmp = f.readline() # blank line

        for line in f:
            data = line.split()
            for i in range(len(data)):
                # This strips off the [,], and ' ' at the ends, and divides it into it's component parts
                if data_format[i] in ["x", "y", "z", "xu", "yu", "zu"]:
                    data[i] = str((float(data[i]) - lower_bounds[data_format[i]]) * args.scale + lower_bounds[data_format[i]])
            update = " ".join(data)
            fout.write("{}\n".format(update))
