#! /usr/bin/env python3

import numpy as np
from sys import argv, exit
import argparse, itertools, contextlib, natsort

parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file file [file ...]', description = "Averages values across multiple files into one file.")
parser.add_argument('file1', metavar = 'file', nargs = 1, help = "The files to average together")
parser.add_argument('file2', metavar = 'file', nargs = '+', help = argparse.SUPPRESS)
parser.add_argument('-n','--not-constant', action = 'store_true', help = "Flag specifying whether the first value is constant across all files or not")
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

with contextlib.ExitStack() as stack:
    fs = [stack.enter_context(open(i, 'r')) for i in args.files]
    for rows in itertools.zip_longest(*fs):
        if all(i.startswith("#") for i in rows):
            continue
        try:
            data = [i.split() for i in rows]
        except:
            break
        if not args.not_constant:
            if not all(v[0] == data[0][0] for v in data):
                print("Error: unable to average values because the first values are different.  If this is expected, include the -n or --not-constant flag.")
                exit(10)
        string = ""
        for i in range(len(data[0])):
            string += "{val} ".format(val=np.mean([float(v[i]) for v in data]))
        f2.write(string + "\n")
f2.close()
