#! /usr/bin/env python3

from sys import argv, exit
import argparse, os
from myModules import *

def is_num(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

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
    if args.output.find(".dat") != -1 or args.output.find(".lmp") != 1: # if the .dat or .lmp extension is found
        output_filename = args.output # the exact filename will be used
    else: # otherwise, we use the basename and check for already existing files.
        output_filename = "{base}.dat".format(base = args.output)
        output_filename = verify_new_file(output_filename)

    if not (os.path.splitext(file)[1] == ".dat" or os.path.splitext(file)[1] == ".lmp"):
        print("Incorrect file type {}".format(os.path.splitext(file)[1]))
        exit(1)

    fout = open(output_filename, 'w')
    scaling_atoms = False # flag to indicate what part of the file we are parsing
    with open(file) as f:
        # edit the comment line to give additional info about file creation
        # read the header data
        lower_bounds = {"x":0.0, "y":0.0, "z":0.0, "xu":0.0, "yu":0.0, "zu":0.0}
        comment = f.readline()
        data_format = comment.split(":")[-1] # we end up with a string that looks like: " [ID type x y z]\n"
        if uniform:
            fout.write(f"The {file} structure has been expanded by a factor of {args.scale[0]}:{data_format}")
        elif xy_only:
            fout.write(f"The {file} structure has been expanded by a factor of x,y * {args.scale[0]}, z * {args.scale[1]}:{data_format}")
        elif scale_all:
            fout.write(f"The {file} structure has been expanded by a factor of x * {args.scale[0]}, y * {args.scale[1]}, z * {args.scale[2]}:{data_format}")
        else:
            print("Error: unknown scaling relationship.")
            exit(1)
        data_format = data_format.strip(" []\n").split()

        # Parse the rest of the header normally: look for specific keywords to apply the scaling
        for line in f:
            if not line or line is "\n":
                fout.write("\n")
                continue
            if scaling_atoms:
                data = line.split()
                for i in range(len(data)):
                    if data_format[i] in ["x", "xu"]:
                        data[i] = f"{(float(data[i]) - lower_bounds[data_format[i]]) * args.scale[0] + lower_bounds[data_format[i]]:.6f}"
                    elif data_format[i] in ["y", "yu"]:
                        if uniform or xy_only:
                            data[i] = f"{(float(data[i]) - lower_bounds[data_format[i]]) * args.scale[0] + lower_bounds[data_format[i]]:.6f}"
                        else: # The only other option is "scale_all"
                            data[i] = f"{(float(data[i]) - lower_bounds[data_format[i]]) * args.scale[1] + lower_bounds[data_format[i]]:.6f}"
                    elif data_format[i] in ["z", "zu"]:
                        if uniform:
                            data[i] = f"{(float(data[i]) - lower_bounds[data_format[i]]) * args.scale[0] + lower_bounds[data_format[i]]:.6f}"
                        elif xy_only:
                            data[i] = f"{(float(data[i]) - lower_bounds[data_format[i]]) * args.scale[1] + lower_bounds[data_format[i]]:.6f}"
                        else:
                            data[i] = f"{(float(data[i]) - lower_bounds[data_format[i]]) * args.scale[2] + lower_bounds[data_format[i]]:.6f}"
                update = " ".join(data)
                fout.write(f"{update}\n")
                continue
            else:
                if scaling_atoms and not all(is_num(i) for i in line.split()): # done with the Atoms section
                    scaling_atoms = False
                    fout.write(line) # we don't change this line
                    continue
                if line.strip().startswith("Atoms"):
                    scaling_atoms = True
                    fout.write(line) # we don't change this line
                    continue

                line_results = [scale_keyword in line for scale_keyword in ["xlo", "ylo", "zlo"]]
                if any(line_results):
                    lo = float(line.split()[0]) # lower bound
                    hi = float(line.split()[1]) # upper bound
                    if line_results[0]:
                        lower_bounds["x"] = lower_bounds["xu"] = lo
                        fout.write(f"{lo:.6f} {(hi - lo) * args.scale[0] + lo:.6f} xlo xhi\n")
                    elif line_results[1]:
                        lower_bounds["y"] = lower_bounds["yu"] = lo
                        if uniform or xy_only: # y direction is scaled by the same factor as x in the uniform and x-y case
                            fout.write(f"{lo:.6f} {(hi - lo) * args.scale[0] + lo:.6f} ylo yhi\n")
                        elif scale_all: # the y direction is scaled by a different factor than x
                            fout.write(f"{lo:.6f} {(hi - lo) * args.scale[1] + lo:.6f} ylo yhi\n")
                    else: # the only other option is line_results[2] == True
                        lower_bounds["z"] = lower_bounds["zu"] = lo
                        if uniform:
                            fout.write(f"{lo:.6f} {(hi - lo) * args.scale[0] + lo:.6f} zlo zhi\n")
                        elif xy_only:
                            fout.write(f"{lo:.6f} {(hi - lo) * args.scale[1] + lo:.6f} zlo zhi\n")
                        else: # otherwise each dimension is scaled separately, so
                            fout.write(f"{lo:.6f} {(hi - lo) * args.scale[2] + lo:.6f} zlo zhi\n")
                else:
                    fout.write(line)
