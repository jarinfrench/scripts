#! /usr/bin/env python3

# Simple script to parse the orientation matrices database and grab the Euler
# angles.

from sys import argv
from myModules import *
import argparse

parser = argparse.ArgumentParser(usage = '%(prog)s file', description = "Parses the orientation matrix database to extract the Euler angles")
parser.add_argument('file', help = "The database file")

args = parser.parse_args()

f1 = open(args.file,'r')

count = 0

while True:
    data = f1.readline().split()
    if not data:
        break
    elif not len(data) == 6:
        continue
    else:
        count = count + 1
        axis = data[0][1:4]
        try:
            degree = int(data[0][9:12]) - 1
        except:
            try:
                degree = int(data[0][9:11]) - 1
            except:
                degree = int(data[0][9]) - 1
        filename2 = 'angle_files/' + str(degree) + '_about_' + axis + '.tex'
        if count == 2:
            count = 0
        if count == 1:
            f2 = open(filename2, 'w')
            f2.write("Euler angles\n")
            f2.write("Angles for a %d degree misorientation\n"%(degree))
            f2.write("Axis = " + axis + '\n')
            f2.write("B100\n")
        else:
            f2 = open(filename2, 'a')
        z1 = str(rad2deg(float(data[3][1:])))
        x = str(rad2deg(float(data[4])))
        z2 = str(rad2deg(float(data[5])))

        f2.write(z1 + ' ' + x + ' ' + z2 + ' 01.00\n')
        f2.close()

f1.close()
