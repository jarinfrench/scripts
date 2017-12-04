#! /usr/bin/python

from __future__ import division, print_function
from itertools import izip
import numpy as np
from sys import argv, exit

# get our filenames and number of files
if len(argv) == 1:
    num_files = int(input("Please specify the number of files to average together: "))
    files = [None] * num_files
    for i in range(num_files):
        files[i] = raw_input("Please specify file " + repr(i+1) + ": ")
else:
    num_files = len(argv) - 1
    files = [None] * num_files
    for i in range(num_files):
        files[i] = argv[i+1]

if num_files == 1:
    print("Cannot average only one file")
    exit(10)

# open up a file stream for each file
f = [open(i, 'r') for i in files]
diff = []
for j in range(num_files - 1):
    diff.append([i for i in xrange(len(files[0])) if files[0][i] != files[j+1][i]])
if min(min(diff)) == 0:
    average_data_file = "average.txt"
else:
    average_data_file = files[0][:min(min(diff))] + "average.txt"

f2 = open(average_data_file,'w')

# Now let's read the first line of each file and discard it (it's a comment line)
# for i in range(num_files):
    # comment = f[i].readline()

# Now let's read the actual data
data = [None] * num_files
for rows in izip(*f): # for each row in each file
    if not all(a.startswith("#") for a in rows):
        for i in range(num_files):
            temp = rows[i].split() # Split the data by spaces
            data_app = []
            for j in range(len(temp)): # for however many elements there are in temp
                data_app.append(float(temp[j])) # combine them into one list as a float
                data[i] = data_app # and put that list into data
        if all(v[0] == data[0][0] for v in data):
            string = "%d"%data[0][0]
            for i in range(len(temp)-1):
                string += " %7.5f"%(np.mean([v[i+1] for v in data]))
            f2.write(string+"\n")
        else:
            print("Error: unable to average values because timesteps are different!")
            exit(10)
for i in range(num_files):
    f[i].close()
f2.close()
