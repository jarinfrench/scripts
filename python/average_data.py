#! /usr/bin/env python3

import numpy as np
import scipy as sp
from sys import argv, exit
import argparse, itertools, contextlib, natsort

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file file [file ...]', description = "Averages values across multiple files into one file.")
parser.add_argument('file1', metavar = 'file', nargs = 1, help = "The files to average together")
parser.add_argument('file2', metavar = 'file', nargs = '+', help = argparse.SUPPRESS)
parser.add_argument('-n','--not-constant', action = 'store_true', help = "Flag specifying whether the first value is constant across all files or not")
parser.add_argument('-s', '--style', choices = ['a','h','g'], default = 'a', help = "Averaging style - (a)rithmetic (default), (h)armonic, or (g)eometric")
args = parser.parse_args()
args.files = natsort.natsorted(args.file1 + args.file2)

diff = []
for j in range(len(args.files) - 1):
    diff.append([i for i in range(len(args.files[0])) if args.files[0][i] != args.files[j+1][i]])
if min(min(diff)) == 0:
    average_data_file = "average.txt"
else:
    average_data_file = args.files[0][:min(min(diff))] + "average.txt"

f2 = open(average_data_file,'w')
warn_string_base = "\n# Empty line(s) found in file "
warn = False

with contextlib.ExitStack() as stack:
    fs = [stack.enter_context(open(i, 'r')) for i in args.files]
    for rows in itertools.zip_longest(*fs):
        try:
            cont = all(i.startswith("#") for i in rows)
        except: # If we have reached the end of any file, we stop averaging
            break
        if cont:
            continue
        try:
            for n,i in enumerate(rows):
                while not i.split() or i.startswith("#"):
                    if not warn:
                        warn_string = warn_string_base + args.files[n] + "\n"
                        warn = True
                    tmp = list(rows)
                    tmp[n] = fs[n].readline()
                    rows = tuple(tmp)
                    i = rows[n]
            data = [i.split() for i in rows]
        except:
            break
        if not args.not_constant:
            if not all(v[0] == data[0][0] for v in data):
                print("Error: unable to average values because the first values are different.  If this is expected, include the -n or --not-constant flag.")
                exit(10)
        if warn:
            f2.write(warn_string)
            warn = False
        string = ""
        for i in range(len(data[0])):
            if args.style == 'a':
                string += "{val} ".format(val=np.mean([float(v[i]) for v in data]))
            elif args.style == 'h':
                string += "{val} ".format(val=sp.stats.hmean([float(v[i]) for v in data]))
            elif args.style == 'g':
                string += "{val} ".format(val=sp.stats.mstats.gmean([float(v[i]) for v in data]))
            else:
                print("Incorrect averaging parameter passed.  Something is wrong!  You entered {}, but we need 'a', 'h' or 'g'".format(args.style))
        f2.write(string + "\n")
f2.close()
