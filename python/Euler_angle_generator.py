#! /usr/bin/env python3
#
# This script will generate a series of Euler Angle files for any arbitrary set
# inputs:
# angle - the maximum angle (starting from 0) that the file will generate
# axis - the axis of interest.  Generally, this will be 100, 110, or 111
# type - Twist or tilt boundary
#      - NOTE: This could possibly be made to handle ANY arbitrary axis
# Author: Jarin French
# This script is based off of work done by John-Michael Bradley

from sys import argv, exit
from math import atan
from numbers import Number
import argparse

# Helper functions
def getTiltType():
    try:
        _tilt_type = input("Is this a symmetric tilt (symm) or an asymmetric tilt (asymm)? ")
    except:
        _tilt_type = -1 # error code to run the while loop again
    return _tilt_type

def writeHeader(i, _axis, tex_file_base):
    # This write the four-line header ignored by MOOSE objects that contains
    # information about the angles.
    deg = float(i)
    tex_filename = tex_file_base + "%d.tex" % i
    tex_file = open(tex_filename, 'w')
    tex_file.write("Angles for a " + str(deg) + " degree symmetric tilt\n")
    tex_file.write("Axis = [%d]\n\n" % _axis)
    tex_file.write("B 2\n")
    return tex_file

# See the argparse docs for a better way to do this.
def check_range(val):
    try:
        value = int(val)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err))

    # A max angle of 0 won't print anything helpful.
    if value < 1 or value >= 360:
        message = "Expected 1 <= angle < 360, got value = {}".format(value)
        raise argparse.ArgumentTypeError(message)

    return value

def three_digit_int(val):
    try:
        value = int(val)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err))

    if len(str(value)) > 3:
        message = "Expected three digit integer, got {} digit integer".format(len(str(value)))
        raise argparse.ArgumentTypeError(message)

    return value


parser = argparse.ArgumentParser(usage = "%(prog)s [-h] angle axis type",
    description = "Generates the Euler angles for a GB given the axis, angle, and GB type")
parser.add_argument("angle", type = check_range, help = "The maximum angle")
parser.add_argument("axis", type = three_digit_int, help = "The rotation axis")
parser.add_argument("type", choices = ["twist", "tilt"], help = "The type of grain boundary")

args = parser.parse_args()

axis = [] # This takes the individual components of args.axis as part of a vector

# Determine what axis we're looking at by reading the string version of args.axis char by char
for i in range(3 - len(str(args.axis))):
    axis.append(0.0)
axis = axis + [float(i) for i in str(args.axis)]

# TODO: This can be taken out, but it requires some changes in the non-basic set analysis
# Sort the axis from greatest to least.  x-axis first, y-axis second, z-axis third
for i in range(0,len(axis)):
    for j in range(i + 1, len(axis)):
        if axis[j] > axis[i]:
            axis[i], axis[j] = axis[j], axis[i]

# This abstract line says: if 1 is equal to the length of the resulting list
# when you take out all the zeros, you have a x00 axis.
is_x00 = (1 == len([value for value in axis if value != 0]))
is_xx0 = (2 == len([value for value in axis if value != 0])) # same here, but checking for xx0
is_xxx = (3 == len([value for value in axis if value != 0])) # ditto, but xxx

# Normalize the axes <-- TODO: Needs work for anything besides the basic 100, 110, and 111 sets
if is_x00: # No matter what number is in the x position, it can always be normalized to be 100
    axis[0] = 1
    is100 = True
    is110 = False
    is111 = False
else:
    is100 = False

if is_xx0:
    # check to see if the x's are the same number
    if axis[0] == axis[1]:
        axis[0],axis[1] = 1,1
        is110 = True
        is111 = False
    else:
        is110 = False

if is_xxx:
    # Check to see if the x's are the same number
    # This line says "Count the number of times that the first element equals
    # the other elements, and check if that's equal to the length"
    if axis.count(axis[0]) == len(axis):
        axis = [1] * 3 # create a list of 1's
        is111 = True
    else:
        is111 = False

if not True in {is100, is110, is111}: # Not a basic axis set
    if not 0 in axis:
        # This checks to see if the axis can be scaled down at all (i.e.
        # 642 can become 321).  This is specific to an axis without 0's
        if axis[0] % axis[2] == 0 and axis[1] % axis[2] == 0:
            axis[0] = axis[0] / axis[2]
            axis[1] = axis[1] / axis[2]
            axis[2] = 1
    else: # Same thing, but for the xx0 axis
        if axis[0] % axis[1] == 0:
            axis[0] = axis[0] / axis[1]
            axis[1] = 1

# Determine the base file name.  Extra word if its a tilt boundary
# NOTE: Asymmetric tilt may not be possible with this code.
if args.type == 'tilt':
    _tilt_type = None
    while not _tilt_type in {"symm", "asymm"}:
        if _tilt_type == -1:
            print("Please type either \"symm\" or \"asymm\"")
        _tilt_type = getTiltType()
    if _tilt_type == 'asymm':
        asymm_weight = float(input("How strong is the asymmetry? "))
        while asymm_weight > 1 or asymm_weight < 0:
            print("Please input a number between 0 and 1 inclusive.  Note that an asymmetry of 1 will be treated as a twist boundary.")
            asymm_weight = float(input("How strong is the asymmetry? "))
        tex_file_base = "{axis}_{tilt_type}_{type}_weight_{:.2f}".format(str(args.axis).zfill(3), _tilt_type, args.type, asymm_weight)
    else:
        tex_file_base = "{}_{}_{}".format(str(args.axis).zfill(3), _tilt_type, args.type)
else:
    tex_file_base = "{}_{}" %(str(args.axis).zfill(3), args.type) # Base for the file name

# The Euler Angles will now be calculated

# First, check to see if twist or tilt
if args.type == 'tilt':
    # If we're only rotating about the x-axis (only a non-zero value in the first spot),
    # the rotation is really easy.
    # Simple cases first
    # NOTE: The 'asymm' tilt is actually used for the twist rotations! Asymmetric
    # tilt is a 3D plot, and (so far) cannot be represented here
    if is100:
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)
            if _tilt_type == 'symm':
                tex_file.write("00.00 " + str(i / 2.0) + " 00.00 1.00\n")
                tex_file.write("00.00 " + str(-i / 2.0) + " 00.00 1.00\n")
            elif _tilt_type == 'asymm':
                phi1 = i * asymm_weight
                phi2 = i - phi1
                tex_file.write("00.00 " + str(phi1) + " 00.00 1.00\n")
                tex_file.write("00.00 " + str(phi2) + " 00.00 1.00\n")
        tex_file.close()

    elif is110:
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)
            if _tilt_type == 'symm':
                tex_file.write("45.00 " + str(i/2.00) + " 00.00 1.00\n")
                tex_file.write("45.00 " + str(-i/2.00) + " 00.00 1.00\n")
            elif _tilt_type == 'asymm':
                phi1 = i * asymm_weight
                phi2 = i - phi1
                tex_file.write("45.00 " + str(phi1) + " 00.00 1.00\n")
                tex_file.write("45.00 " + str(phi2) + " 00.00 1.00\n")
        tex_file.close()

    elif is111:
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)
            if _tilt_type == 'symm':
                phi_plus = 37.9381 + i / 2.00
                phi_minus = 37.9381 - i / 2.00
                # Map everything back to the 0-360 equivalent space
                if phi_plus > 360:
                    phi_plus = phi_plus - 360.00
                if phi_minus < 0:
                    phi_minus = phi_minus + 360.00
                tex_file.write("290.104 " + str(phi_plus) + " 110.104 1.00\n")
                tex_file.write("290.104 " + str(phi_minus) + " 110.104 1.00\n")

            elif _tilt_type == 'asymm':
                i1 = i * asymm_weight
                i2 = i - i1
                phi1 = 37.9381 + i1
                phi2 = 37.9381 - i2
                if phi1 > 360:
                    phi1 = phi1 - 360.00
                if phi2 > 360:
                    phi2 = phi2 - 360.00
                if phi1 < 0:
                    phi1 = phi1 + 360.00
                if phi2 < 0:
                    phi2 = phi2 + 360.00
                tex_file.write("290.104 " + str(phi1) + " 110.104 1.00\n")
                tex_file.write("290.104 " + str(phi2) + " 110.104 1.00\n")
        tex_file.close()

    # Now the more general cases
    elif is_xx0:
        # Calculate how much around the Z axis we need to rotate
        z_rotation = atan(axis[1]/axis[0]);
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)
            if _tilt_type == 'symm':
                tex_file.write(str(z_rotation) + " " + str(i/2.00) + " 00.00 1.00\n")
                tex_file.write(str(z_rotation) + " " + str(-i/2.00) + " 00.00 1.00\n")
            elif _tilt_type == 'asymm':
                tex_file.write(str(z_rotation) + " " + str(i) + " 00.00 1.00\n")
                tex_file.write(str(z_rotation) + " 00.00 00.00 1.00\n")
            tex_file.close()

    elif is_xxx:
        print("Sorry, this case has not yet been implemented")
        exit()

    else: # Just in case something weird happens
        print("ERROR: Unrecognized axis")
        exit()

elif args.type == 'twist':
    if is100:
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)
            tex_file.write("00.00 " + str(i) + " 00.00 1.00\n")
            tex_file.write("00.00 00.00 00.00 1.00\n")
        tex_file.close()

    elif is110:
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)
            tex_file.write("45.00 " + str(i) + " 00.00 1.00\n")
            tex_file.write("45.00 00.00 00.00 1.00\n")
        tex_file.close()

    elif is111:
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)

    elif is_xx0:
        for i in range(0, args.angle + 1):
            tex_file = writeHeader(i, args.axis, tex_file_base)

    elif is_xxx:
        print("Sorry, this case has not yet been implemented")
        exit()

    else:
        print("ERROR: Unrecognized axis")
        exit()

else:
    print("ERROR: Unrecognized boundary type")
