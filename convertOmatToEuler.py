#! /opt/moose/miniconda/bin/python
#
# This script accepts as an argument a list of 9 numbers. These 9 numbers are
# assumed to represent, *in order*, the values in the tensor:
# _o_matrix =
#    arg1   arg2   arg3
#    arg4   arg5   arg6
#    arg7   arg8   arg9
#
# These values make up the orientation matrix of a grain.  This script will
# then calculate the Euler angles.
#
# Options:
#   -q --quiet                    Suppresses the output of this script
#   --help                        Displays this information
#
# Output:
#  The output is an Euler angle set in the form phi_1, Phi, phi_2, representing
#  the Bunge convention of crystallographic angles (ZXZ).

from __future__ import print_function, division
from math import sin, cos, acos, atan2, pi
from numpy import array, linalg
from sys import argv
from myModules import *

# Helper functions
def displayHelp():
    print('''
    This script accepts as an argument a list of 9 numbers. These 9 numbers are
    assumed to represent, *in order*, the values in the tensor:
    _o_matrix =
        arg1   arg2   arg3
        arg4   arg5   arg6
        arg7   arg8   arg9

    These values make up the orientation matrix of a grain.  This script will
    then calculate the Euler angles.

    Options:
       -q --quiet                    Suppresses the output of this script
       --help                        Displays this information

    Output:
      The output is an Euler angle set in the form phi_1, Phi, phi_2, representing
      the Bunge convention of crystallographic angles (ZXZ).
    ''')
    return

quiet, argv = check4Quiet(argv) # Checks for suppressing output

if "--help" in argv: # Help info
    displayHelp()
    exit()

if len(argv) < 10:
    print("ERROR: Not enough arguments.")
    displayHelp()
    exit()
else:
    script, a11, a12, a13, a21, a22, a23, a31, a32, a33 = argv

# make sure 0 is positive, otherwise we get problems with atan2
a = array([[float(a11), float(a12), float(a13)],[float(a21), float(a22), float(a23)], [float(a31), float(a32), float(a33)]])
# a = a / linalg.norm(a) # NOTE: normalizing produces INCORRECT results
if a[2][2] == 0:
    a[2][2] = abs(a[2][2])

if a[0][2] == 0:
    a[0][2] = abs(a[0][2])

if a[1][2] == 0:
    a[1][2] = abs(a[1][2])

if a[2][0] == 0:
    a[2][0] = abs(a[2][0])

if a[2][1] == 0:
    a[2][1] = abs(a[2][1])

# This next part assumes that the rotation matrix is a passive rotation.
Phi = acos(a[2][2])
phi_2 = atan2(a[0][2],a[1][2])
if a[2][1] == 0:
    phi_1 = atan2(a[2][0], abs(a[2][1])) # -0 causes an issue with atan2
else:
    phi_1 = atan2(a[2][0], -a[2][1])

print("Euler angles are %2.4f, %2.4f, and %2.4f"%(rad2deg(phi_1), rad2deg(Phi), rad2deg(phi_2)))
# Double check that this creates the same matrix.
c1 = cos(phi_1)
c2 = cos(Phi)
c3 = cos(phi_2)
s1 = sin(phi_1)
s2 = sin(Phi)
s3 = sin(phi_2)

b11 = round(c1*c3 - c2*s1*s3, 9)
b12 = round(c3*s1 + c1*c2*s3, 9)
b13 = round(s2*s3, 9)
b21 = round(-c1*s3 - c2*c3*s1, 9)
b22 = round(c1*c2*c3 - s1*s3, 9)
b23 = round(c3*s2, 9)
b31 = round(s1*s2, 9)
b32 = round(-c1*s2, 9)
b33 = round(c2, 9)

b = array([[b11, b12, b13],[b21, b22, b23], [b31, b32, b33]])
#b = b / linalg.norm(b)
result = [[b11 == a[0][0], b12 == a[0][1], b13 == 10], [b21 == a[1][0], b22 == a[1][1], b23 == a[1][2]], [b31 == a[2][0], b32 == a[2][1], b33 == a[2][2]]]

# TODO: determine how to display which results do not match.
# Still some work to do....
