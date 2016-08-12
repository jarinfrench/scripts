#! /opt/moose/miniconda/bin/python
#
# This script will generate a series of Euler Angle files for any arbitrary set
# inputs:
# _max_angle - the maximum angle (starting from 0) that the file will generate
# _axis - the axis of interest.  Generally, this will be 100, 110, or 111
# _type - Twist or tilt boundary
#      - NOTE: This could possibly be made to handle ANY arbitrary axis
# Author: Jarin French
# This script is based off of work done by John-Michael Bradley

from sys import argv
from math import atan
from numbers import Number

# Helper functions
def getMaxAngle():
    try:
        max_angle = int(input("Please enter the max angle: "))
    except:
        print("Please make sure to enter only numbers")
        max_angle = -1
    return max_angle

def getMisorientationType():
    try:
        _type = input("Please enter the type of misorientation (\"twist\" or \"tilt\"): ")
    except:
        printQuoteError()
        _type = -1
    return _type

def getAxis():
    try:
        _axis = int(input("Please enter the axis of rotation: "))
    except:
        print("Please make sure to enter only numbers")
        _axis = -1
    return _axis

def printQuoteError():
    print("ERROR: You probably forgot to put your response in quotes")
    return

def largeAngleWarning():
    print("WARNING: Maximum angle is over the symmetry limit.")
    return

def promptChangeAngle():
    try:
        change_angle = input("Would you like to lower the angle to the symmetry limit? ")
    except:
        printQuoteError()
        change_angle = -1
    return change_angle

def getTiltType():
    try:
        _tilt_type = input("Is this a symmetric tilt (symm) or an asymmetric tilt (asymm)? ")
    except:
        printQuoteError()
        _tilt_type = -1
    return _tilt_type

def writeHeader(i, _axis, tex_file_base):
    deg = float(i)
    tex_filename = tex_file_base + "%d.tex" % i
    tex_file = open(tex_filename, 'w')
    tex_file.write("Angles for a " + str(deg) + " degree symmetric tilt\n")
    tex_file.write("Axis = [%d]\n\n" % _axis)
    tex_file.write("B 2\n")
    return tex_file

# Input error checking - Look into using Python's built in methods?
if len(argv) != 4:
    print("WARNING: Not enough command line arguments.")
    _max_angle = 0
    _axis = 0
    _type = None
    while _max_angle < 1:
        if _max_angle == -1:
            print("Please enter an integer between 1 and 360")
        _max_angle = getMaxAngle()

    while _axis < 1:
        if _axis == -1:
            print("Please enter an integer of length 3 (i.e. 100 or 111)")
        _axis = getAxis()

    while not _type in {'twist','tilt'}:
        if _type == -1 or isinstance(_type, Number):
            print("Please type \"twist\" or \"tilt\"")
        _type = getMisorientationType()
else:
    script, _max_angle, _axis, _type = argv

axis = [None]*3 # This takes the individual components of _axis as part of a vector

if len(str(_axis)) > 3: # axis length greater than 3
    print("ERROR: Argument 2 must by a 3 digit number like \'100\'")
    exit()

if not isinstance(_max_angle, Number): # max angle not a number
    print("ERROR: Argument 1 must be a number")
    exit()

if type(_max_angle) != int: # max angle not an integer <-- may remove this... I will probably use floats eventually
    print("ERROR: Argument 1 must be an integer")
    exit()

if _max_angle < 0: # max angle is negative
    print("ERROR: Argument 1 must be positive")
    exit()

if not _type in ['twist', 'tilt']: # type is not twist or tilt
    print("ERROR: _type must be either \'twist\' or \'tilt\'") # TODO: Is it possible to figure out a mixed boundary euler angle?
    exit()

# Warnings
if len(str(_axis)) < 3: # axis length less than 3
    print("WARNING: Assumed trailing zeros")
    for i in range(0, len(str(_axis))):
        axis[i] = int(str(_axis)[i])
    trailing_zeros = 3 - len(str(_axis))
    for i in range(0, int(trailing_zeros)):
        axis[i + 1] = 0

if _max_angle > 360:
    largeAngleWarning()
    change_angle = None
    while not change_angle.lower().startswith('y') or change_angle.lower().startswith('n'):
        if change_angle == -1:
            print("Please type \'y\' or \'n\'")
        change_angle = promptChangeAngle()
    if change_angle.lower().startswith('y'):
        _max_angle = 360

# Determine what axis we're looking at by reading the string version of _axis char by char
for i in range(0,len(str(_axis))):
    axis[i] = float(str(_axis)[i])

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
if is_x00:
    axis[0] = 1
    is100 = True
else:
    is100 = False

if is_xx0:
    # check to see if the x's are the same number
    if axis[0] == axis[1]:
        axis[0],axis[1] = 1,1
        is110 = True
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
        if axis[0] % axis[2] == 0 and axis[1] % axis[2] == 0:
            axis[0] = axis[0] / axis[2]
            axis[1] = axis[1] / axis[2]
            axis[2] = 1
    else:
        if axis[0] % axis[1] == 0:
            axis[0] = axis[0] / axis[1]
            axis[1] = 1

# Change _axis for the filename
_axis = int(str(int(axis[0])) + str(int(axis[1])) + str(int(axis[2])))

# Determine the base file name.  Extra word if its a tilt boundary
if _type == 'tilt':
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
        tex_file_base = "%d_%s_%s_weight_%2.2f" %(_axis, _tilt_type, _type, asymm_weight)
    else:
        tex_file_base = "%d_%s_%s_" %(_axis, _tilt_type, _type)
else:
    tex_file_base = "%d_%s_" %(_axis, _type) # Base for the file name

# The Euler Angles will now be calculated

# First, check to see if twist or tilt
if _type == 'tilt':
    # If we're only rotating about the x-axis (only a non-zero value in the first spot),
    # the rotation is really easy.
    # Simple cases first
    # NOTE: The 'asymm' tilt is actually used for the twist rotations! Asymmetric
    # tilt is a 3D plot, and (so far) cannot be represented here
    # NOTE: Perhaps the 'asymm' tilt can be represented by a factor between 0 and
    # 1? This could weight the differences heavier to one side, instead of just
    # having everything depend on the one grain.
    if is100:
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)
            if _tilt_type == 'symm':
                tex_file.write("00.00 " + str(i/2) + " 00.00 1.00\n")
                tex_file.write("00.00 " + str(-i/2) + " 00.00 1.00\n")
            elif _tilt_type == 'asymm':
                phi1 = i * asymm_weight
                phi2 = i - phi1
                tex_file.write("00.00 " + str(phi1) + " 00.00 1.00\n")
                tex_file.write("00.00 " + str(phi2) + " 00.00 1.00\n")
        tex_file.close()

    elif is110:
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)
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
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)
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
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)
            if _tilt_type == 'symm':
                tex_file.write(str(z_rotation) + " " + str(i/2.00) + " 00.00 1.00\n")
                tex_file.write(str(z_rotation) + " " + str(-i/2.00) + " 00.00 1.00\n")
            elif _tilt_type == 'asymm':
                tex_file.write(str(z_rotation) + " " + str(i) + " 00.00 1.00\n")
                tex_file.write(str(z_rotation) + " 00.00 00.00 1.00\n")
            tex_file.close()

    elif is_xxx:
        print("Sorry, I haven't figured this out yet")
        exit()

    else:
        print("ERROR: Unrecognized axis")
        exit()

elif _type == 'twist':
    if is100:
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)
            tex_file.write("00.00 " + str(i) + " 00.00 1.00\n")
            tex_file.write("00.00 00.00 00.00 1.00\n")
        tex_file.close()

    elif is110:
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)
            tex_file.write("45.00 " + str(i) + " 00.00 1.00\n")
            tex_file.write("45.00 00.00 00.00 1.00\n")
        tex_file.close()

    elif is111:
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)

    elif is_xx0:
        for i in range(0, _max_angle + 1):
            tex_file = writeHeader(i, _axis, tex_file_base)

    elif is_xxx:
        print("Sorry, I haven't figured this out yet")
        exit()

    else:
        print("ERROR: Unrecognized axis")
        exit()

else:
    print("ERROR: Unrecognized boundary type")
