#! /usr/bin/env python3

from sys import argv, exit
import argparse, os
from myModules import *

file_format = "LAMMPS_<element(s)>_N<number_of_atoms>[_<extra_info>].dat"

def determineElements(filename):
    # Note that the filename must be in the format LAMMPS_<element(s)>_N<number_of_atoms>[_<extra_info>].<extension>
    substrings = filename.split('_')
    if not substrings[0] == "LAMMPS" or not substrings[2].startswith("N"):
        print("File not formatted correctly.  File must be in the format {}".format(file_format))
        exit(1)
    else:
        return substrings[1]

def getFittedPotentials(db_file):
    fits = dict()
    with open(db_file) as f:
        for line in nonblank_lines(f):
            if line.startswith("#"):
                continue
            else:
                data = line.split()
                fits[data[0]] = [float(i) for i in data[1:]]
                if not len(fits[data[0]]) == 8:
                    print("Database entry for {} not formatted correctly.".format(data[0]))
                    del fits[data[0]]
    return fits

def calculateScaleFactor(fit, T, a0):
    if T >= fit[0] and T <= fit[1]:
        A = fit[2]
        B = fit[3]
        C = 0.0
    elif T > fit[1] and T <= fit[4]:
        A = fit[5]
        B = fit[6]
        C = fit[7]
    else:
        print("Temperature out of fitted range ({Tmin} K - {Tmax} K)".format(Tmin = fit[0], Tmax = fit[4]))
        exit(2)

    return (A + B * T + C * T**2)/a0


parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file T potential a0 [file file file ...] [options]', description = "Expands a lattice by a constant (a, b, and c can be specified separately).")
parser.add_argument('file', nargs = '+', help = "The file(s) to process")
parser.add_argument('T', type = int, help = "The temperature to expand the lattice at.")
parser.add_argument('potential', help = "The name of the potential (as specified in the database file)")
parser.add_argument('a0', type = float, help = "lattice constant at 0K (or at lowest stable T)")
#parser.add_argument('-n', '--no-element', action = 'store_true', help = "Flag specifying that the file name does not contain the element, and thus should not be looked for")
parser.add_argument('-a', type = float, help = "New lattice constant 1 (x direction)")
parser.add_argument('-b', type = float, help = "New lattice constant 2 (y direction)")
parser.add_argument('-c', type = float, help = "New lattice constant 3 (z direction)")
parser.add_argument('-d','--database', help = "location of database file to use", default = os.environ['HOME'] + "/projects/scripts/lattice_params.db")

args = parser.parse_args()

fits = getFittedPotentials(args.database)
fit  = fits[args.potential]

for file in args.file:
    output_filename = "{base}_T{temp}_expanded.dat".format(base = os.path.splitext(file)[0], temp = args.T)
    scale_factor = calculateScaleFactor(fit, args.T, args.a0)
    fout = open(output_filename, 'w')
    element = determineElements(file)
    with open(file) as f:
        # Read the data.  Header formats change based on file extension, so:
        if not os.path.splitext(file)[1] == ".dat":
            print("File type {} not handled by this script.  Please convert to .dat (LAMMPS) format.".format(os.path.splitext(file)[1]))
            os.exit(1)

        # read the header data
        lower_bounds = {"x":0.0, "y":0.0, "z":0.0, "xu":0.0, "yu":0.0, "zu":0.0}
        comment = f.readline()
        data_format = comment.split(":")[1] # we end up with a string that looks like: " [ID type x y z]"
        fout.write("This {el} structure has been expanded at T = {temp} K:{data_str}\n".format(el = element, temp = args.T, data_str = data_format))
        data_format = data_format.strip("[] ").split()
        tmp = f.readline() # blank line
        n_atoms = f.readline() # number of atoms
        fout.write("{}\n".format(n_atoms))
        n_types = f.readline() # number of atom types
        fout.write("{}\n".format(n_types))
        tmp = f.readline() # x bounds
        xlo = float(tmp.split()[0]) # lower bound
        xhi = float(tmp.split()[1]) # upper bound
        lower_bounds["x"] = xlo
        lower_bounds["xu"] = xlo
        fout.write("{} {} xlo xhi\n".format(xlo, (xhi-xlo) * scale_factor + xlo))
        tmp = f.readline() # y bounds
        ylo = float(tmp.split()[0]) # lower bound
        yhi = float(tmp.split()[1]) # upper bound
        lower_bounds["y"] = ylo
        lower_bounds["yu"] = ylo
        fout.write("{} {} ylo yhi\n".format(ylo, (yhi-ylo) * scale_factor + ylo))
        tmp = f.readline() # z bounds
        zlo = float(tmp.split()[0]) # lower bound
        zhi = float(tmp.split()[1]) # upper bound
        lower_bounds["z"] = zlo
        lower_bounds["zu"] = zlo
        fout.write("{} {} zlo zhi\n".format(zlo, (zhi-zlo) * scale_factor + zlo))
        tmp = f.readline() # blank line
        tmp = f.readline() # "Atoms"
        fout.write("{}\n".format(tmp))
        tmp = f.readline() # blank line

        for line in f:
            data = line.split()
            for i in range(len(data)):
                # This strips off the [,], and ' ' at the ends, and divides it into it's component parts
                if data_format[i] in ["x", "y", "z", "xu", "yu", "zu"]:
                    data[i] = str((float(data[i]) - lower_bounds[data_format[i]]) * scale_factor + lower_bounds[data_format[i]])
            update = " ".join(data)
            fout.write("{}\n".format(update))
    fout.close()
