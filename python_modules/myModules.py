# This file contains all of the functions that I have found useful over time.
# A short description of each function precedes the definition.
from __future__ import division,print_function # makes division and printing easier
from math import sin, cos, pi
from itertools import takewhile,repeat

# This function converts degrees to radians.  Argument is assumed to be in degrees
def deg2rad(x):
    return x * pi / 180.0

# This function converts radians to degrees.  Argument is assumed to be in radians
def rad2deg(x):
    return x * 180.0 / pi

# This function displays a matrix in an easy-to-see format.  This will display
# a matrix of any size and dimension
def displayMat(m):
    for i in range(0, len(m)):
        for j in range(0,len(m[i])):
            print("%2.6f\t"%(m[i][j]),end="")
        print("\n")
    print("\n")
    return

# This function checks an argument list for a -s or a --save command to save
# something (generally used for command-line arguments)
def check4Save(args):
    if "-s" in args or "--save" in args: # We want to save the file
        try:
            index = args.index("-s")
        except:
            index = args.index("--save")
        del args[index]
        return True, args
    else:
        return False, args

# Similar to the check4Save function, this function checks an argument list for
# a -q or a --quiet command to suppress terminal output (generally used for
# command line arguments).
def check4Quiet(args):
    if "-q" in args or "--quiet" in args: # We don't want terminal output
        try:
            index = args.index("-q")
        except:
            index = args.index("--quiet")
        del args[index]
        return True, args
    else:
        return False, args

# This function counts the number of lines in a file
def countFileLines(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen if buf )
