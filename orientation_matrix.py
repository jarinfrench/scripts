#! /opt/moose/miniconda/bin/python

# This script will calculate the orientation matrices for any given misorientation
# for any of the high-symmetry axes.
# Arguments:
#
#   _axis: The axis of orientation (type: int)
#   _misorientation: The angle of misorientation (type: float)
#   _gbnorm: The grain boundary normal of the two grains (type: int)
#            --------OR--------
#    (with option -e or --euler)
#   _z1: The first rotation angle  (Z ) (type: float)
#   _x:  The second rotation angle (X') (type: float)
#   _z2: The third rotation angle  (Z") (type: float)
#
# If the option -e or --euler-angles is entered, the calculation skips to simply
# output the orientation matrices.  Otherwise, the Euler angles are calculated from
# the axis, orientation, and grain boundary normal, and then the orientation matrix is
# created through the use of the Rodrigues Rotation Formula, which is:
# R = I + sin(theta) * K + (1 - cos(theta))*K^2
# where I is the identity matrix, theta is the misorientation angle, and K is
# the skew-symmetric matrix formed by the axis of rotation:
# K = 0  -kz  ky
#     kz  0  -kx
#    -ky  kx  0
# where the vector k is the unit vector defining the axis of rotation, or using
# a set of predefined rotations for each axis (default is the predefined rotations).
# The Euler angles are calculated in this case simply for the file to be written
# to.  If the user does not specify to save, then the angles are not used for
# anything.
#
# Options:
# -e --euler <_z1> <_x> <_z2>           Returns the  Bunge orientation matrix
#                                       based on the euler angles provided.
#
# -f --file <filename>                  Reads the file filename and uses the
#                                       Euler angles from them to calculate the
#                                       orientation matrix.
#
# --rrf                                 Calculates the matrices using the Rodrigues
#                                       Rotation Formula
#
# -a --angles                           Displays the Euler angles.  Can be used
#                                       in conjunctions with -q or --quiet to
#                                       display only the Euler angles.
#
# -s --save                             Saves the resultant orientation matrix to
#                                       a database (orientation_matrix_database.m)
#                                       with the accompanying Euler angles.
#
# -q --quiet                            Suppresses output of the orientation matrices
#
# --help                                Display this help info
#
# Output:
# For an Euler angle set, the ouput is simply its orientation matrix.
# For the misorientations, the first matrix is the 'P' orientation matrix, and
# the second matrix is the 'Q' orientation matrix (see Bulatov et al., Acta Mater
# 65 (2014) 161-175).

from __future__ import division,print_function # To avoid numerical problems with division, and for ease of printing
from sys import argv # for CLI arguments
from math import cos,sin,pi, acos, asin, atan2, sqrt # Trig functions
from numbers import Number # For checking input parameters
from os.path import exists # For checking existence of a file
from numpy import array, linalg
from myModules import * # imports my functions from the file myModules.py

orientation_matrix = []

# Helper functions
def displayHelp():
    print('''
    This script will calculate the orientation matrices for any given misorientation
    for any of the high-symmetry axes.
    Arguments:

      _axis: The axis of orientation (type: int or list)
      _misorientation: The angle of misorientation (type: float)
      _gbnorm: The grain boundary normal of the two grains (type: int or list)
                --------OR--------
      (with option -e or --euler)
      _z1: The first rotation angle  (Z ) (type: float)
      _x:  The second rotation angle (X') (type: float)
      _z2: The third rotation angle  (Z") (type: float)

    If the option -e or --euler-angles is entered, the calculation skips to simply
    output the orientation matrices.  Otherwise, the Euler angles are calculated from
    the axis, orientation, and misorientation type, and then the orientation matrix is
    created through the use of the Rodrigues Rotation Formula, which is:
    R = I + sin(theta) * K + (1 - cos(theta))*K^2
    where I is the identity matrix, theta is the misorientation angle, and K is
    the skew-symmetric matrix formed by the axis of rotation:
    K = 0  -kz  ky
        kz  0  -kx
       -ky  kx  0
    Where the vector k is the unit vector defining the axis of rotation, or using
    a set of predefined rotations for each axis (default is the predefined rotations).
    The Euler angles are calculated in this case simply for the file to be written
    to.  If the user does not specify to save, then the angles are not used for
    anything.

    Options:
    -e --euler <_z1> <_x> <_z2>           Returns the  Bunge orientation matrix
                                          based on the euler angles provided.

    -f --file <filename>                  Reads the file filename and uses the
                                          Euler angles from them to calculate the
                                          orientation matrix.

    --rrf                                 Calculates the matrices using the Rodrigues
                                          Rotation Formula

    -a --angles                           Displays the Euler angles.  Can be used
                                          in conjunctions with -q or --quiet to
                                          display only the Euler angles.

    -s --save                             Saves the resultant orientation matrix to
                                          a database (orientation_matrix_database.m)
                                          with the accompanying Euler angles.

    -q --quiet                            Suppresses output of the orientation matrices

    --help                                Display this help info

    Output:
    For an Euler angle set, the ouput is simply its orientation matrix.
    For the misorientations, the first matrix is the 'P' orientation matrix, and
    the second matrix is the 'Q' orientation matrix (see Bulatov et al., Acta Mater
    65 (2014) 161-175).
    ''')
    return

def displayAngles(z1, x, z2): # Displays an Euler angle set (Bunge convention)
    print("Euler angles:")
    print("Z\t\tX\t\tZ")
    print("----------------------------------------")
    print("%2.4f\t%2.4f\t\t%2.4f\n\n"%(rad2deg(z1),rad2deg(x),rad2deg(z2)))
    return

def check4RRF(args): # Check the args for the rrf command
    if "--rrf" in args:
        index = args.index("--rrf")
        del args[index]
        return True, args
    else:
        return False, args

def check4Euler(args): # Check the args for the -a or --angles command
    if "-a" in args or "--angles" in args:
        try:
            index = args.index("-a")
        except:
            index = args.index("--angles")
        del args[index]
        return True, args
    else:
        return False, args

# This function determines the type of misorientation: twist, tilt, or mixed
def defineMisorientation(axis,gbnorm):
    # Make sure axis and gbnorm are arrays:
    axis = array(axis)
    gbnorm = array(gbnorm)

    # Make sure the vectors are normalized
    axis_norm = axis / linalg.norm(axis)
    gbnorm_norm = gbnorm / linalg.norm(gbnorm)

    # Now take the dot product
    dotp = axis_norm.dot(gbnorm_norm.T)

    if dotp == 0.0:
        return 'tilt'

    elif dotp == 1.0 or dotp == -1.0:
        return 'twist'

    else:
        return 'mixed'


# Calculates the Bunge orientation matrix.  Arguments are first z rotation, x rotation, and second z rotation
def calcRotMat(_z1,_x,_z2):
    c1 = cos(_z1)
    c2 = cos(_x)
    c3 = cos(_z2)
    s1 = sin(_z1)
    s2 = sin(_x)
    s3 = sin(_z2)

    rot_mat = array([[c1*c3 - c2*s1*s3, -c1*s3 - c2*c3*s1,  s1*s2],
                     [c3*s1 + c1*c2*s3,  c1*c2*c3 - s1*s3, -c1*s2],
                     [s2*s3           ,  c3*s2           ,  c2  ]])
    return rot_mat

def matMult(m1,m2): # Multiplies two matrices together
    result = m1.dot(m2)
    return result

# Write the matrix and angles to a file
def writeMat(m, z1, x, z2, grain, axis):
    # This is to avoid issues with duplicates
    if _z1 == 0:
        _z1 = abs(_z1)
    if _x == 0:
        _x = abs(_x)
    if _z2 == 0:
        _z2 = abs(_z2)

    lastVal = 1
    #tex_filename = "orientation_matrix_database_ints.m"
    #tex_filename = "orientation_matrix_database_ints_rrf.m"
    tex_filename = "orientation_matrix_database.m"
    var_name = "%s%d%s"%(grain, axis, _type)
    #var_name = "%s%d%srrf"%(grain, axis, _type)
    if not exists(tex_filename):
        tex_file = open(tex_filename, "a")
        tex_file.write("%Database for orientation matrices for specified Euler Angles\n")
        tex_file.write("%-------------------------------------------------------------------\n")
        tex_file.write("%Orientation Matrix                                     Euler Angles\n")
        tex_file.write("%s(:,:,%d)=[%2.6f  %2.6f  %2.6f           %%%2.4f  %2.4f  %2.4f\n"%(var_name, lastVal, m[0][0], m[0][1], m[0][2], _z1, _x, _z2))
        tex_file.write("%2.6f  %2.6f  %2.6f\n"%(m[1][0], m[1][1], m[1][2]))
        tex_file.write("%2.6f  %2.6f  %2.6f];\n"%(m[2][0], m[2][1], m[2][2]))
        tex_file.write("%-------------------------------------------------------------------\n")
        tex_file.close()
    else:
        f = open(tex_filename,"r")
        while True:
            data = f.readline().split()
            if not data:
                break
            elif len(data) != 6:
                continue
            else:
                assert data[0][0] in {'P', 'Q'}, "Unknown orientation matrix type (should be \'P\' or \'Q\')"
                if not "%d%s"%(axis,_type) in {data[0][1:8], data[0][1:9]}:
                    lastVal = 0
                elif data[0][5] == 'w': # We're looking at either P or Q1xxtwist(i)
                    try:
                        try:
                            if data[0][0] == 'P': # Handles anything 3 digits long
                                lastVal = int(data[0][14:17]) - 1
                                #lastVal = int(data[0][17:20]) - 1
                            else:
                                lastVal = int(data[0][14:17])
                                #lastVal = int(data[0][17:20])
                        except:
                            if data[0][0] == 'P': # Handles anything 2 digits long
                                lastVal = int(data[0][14:16]) - 1
                                #lastVal = int(data[0][17:19]) - 1
                            else: # data[0][0] == 'Q'
                                lastVal = int(data[0][14:16])
                                #lastVal = int(data[0][17:19])
                    except:
                        if data[0][0] == 'P': # One digit case
                            lastVal = int(data[0][14]) - 1
                            #lastVal = int(data[0][17]) - 1
                        else: # data[0][0] == 'Q'
                            lastVal = int(data[0][14])
                            #lastVal = int(data[0][17])
                elif data[0][5] == 'i': # P or Q1xxtilt(i)
                    try:
                        try:
                            if data[0][0] == 'P':
                                lastVal = int(data[0][13:16]) - 1
                                #lastVal = int(data[0][16:19]) - 1
                            else:
                                lastVal = int(data[0][13:16])
                                #lastVal = int(data[0][16:19])
                        except:
                            if data[0][0] == 'P':
                                lastVal = int(data[0][13:15]) - 1
                                #lastVal = int(data[0][16:18]) - 1
                            else: # data[0][0] == 'Q'
                                lastVal = int(data[0][13:15])
                                #lastVal = int(data[0][16:18])
                    except:
                        if data[0][0] == 'P':
                            lastVal = int(data[0][13]) - 1
                            #lastVal = int(data[0][16]) - 1
                        else: # data[0][0] == 'Q'
                            lastVal = int(data[0][13])
                            #lastVal = int(data[0][16])
                else:
                    print("Error: Unknown last index. Line 286")
                    exit()

                if grain == 'Q' and _type == 'twist': # We run into problems if we're doing twist matrices for the Q grain - As is now, this will cause the 'Q' twist matrices to ALWAYS be written
                    unique = True
                elif data[0][0] == grain and (data[0][4:8] == _type or data[0][4:9] == _type) and data[3] == ('%' + "%2.4f"%_z1) and data[4] == "%2.4f"%_x and data[5] == "%2.4f"%_z2:
                    unique = False
                    break
                else:
                    unique = True
        if unique:
            tex_file = open(tex_filename, "a")
            tex_file.write("%s(:,:,%d)=[%2.6f  %2.6f  %2.6f           %%%2.4f  %2.4f  %2.4f\n"%(var_name, lastVal + 1, m[0][0], m[0][1], m[0][2], _z1, _x, _z2))
            tex_file.write("%2.6f  %2.6f  %2.6f\n"%(m[1][0], m[1][1], m[1][2]))
            tex_file.write("%2.6f  %2.6f  %2.6f];\n"%(m[2][0], m[2][1], m[2][2]))
            tex_file.write("%-------------------------------------------------------------------\n")
            tex_file.close()
    return

if "--help" in argv: # Help info
    displayHelp()
    exit()

save, argv = check4Save(argv) # Save the file?  Delete the save argument
quiet, argv = check4Quiet(argv) # Checks for suppressing output. Delete the quiet argument.
useRRF, argv = check4RRF(argv) # Checks for using the RRF method. Delete the rrf argument.
dispEuler, argv = check4Euler(argv) # Checks for displaying the Euler angles.  Delete the angle argument

# If the arguments come from a file...
if "-f" in argv or "--file" in argv: #input arguments come from file
    try:
        try:
            index = argv.index("-f")
        except:
            index = argv.index("--file")
    except:
        print("ERROR: Unable to find filename. Line 322")
        exit()
    filename = argv[index + 1]
    try:
        f1 = open(filename, 'r')
    except:
        print("ERROR: Unable to read file. Line 328", filename)

    while True: # Read the file line by line.
        line = f1.readline()
        if not line: # break if we don't read anything
            break;
        data = line.split()
        if len(data) != 4: # If there are less than 4 parts to the data, move along (format of file MUST be _z1 _x _z2 1.00)
            continue
        else:
            # Convert the data to stuff we can use
            _z1 = float(data[0])
            _x  = float(data[1])
            _z2 = float(data[2])

            _z1 = deg2rad(_z1)
            _x  = deg2rad(_x)
            _z2 = deg2rad(_z2)
            orientation_matrix = calcRotMat(_z1, _x, _z2)
            if not quiet:
                displayMat(orientation_matrix)
            if save:
                writeMat(orientation_matrix, _z1, _x, _z2,'P', _axis, _type)
# Input is a set of euler angles
elif "-e" in argv or "--euler-angles" in argv:
    try:
        try:
            index = argv.index("-e")
        except:
            index = argv.index("--euler-angles")
    except:
        print("ERROR: Unable to read Euler angles. Line 359")
        exit()
    _z1 = float(argv[index + 1])
    _x =  float(argv[index + 2])
    _z2 = float(argv[index + 3])
    _z1 = deg2rad(_z1)
    _x  = deg2rad(_x)
    _z2 = deg2rad(_z2)

    orientation_matrix = calcRotMat(_z1, _x, _z2)

    if not quiet:
        displayMat(orientation_matrix)
    if save:
        writeMat(orientation_matrix, _z1, _x, _z2, 'P', _axis, _type)

else:
    if len(argv) < 4:
        print("ERROR: Not enough command line arguments. Line 377")
        print("Input either an axis, misorientation, and misorientation type, or a ZXZ Euler angle set with the option -e or --euler-angles.")
        displayHelp()
        exit()
    try:
        _axis = int(argv[1])
        _misorientation = float(argv[2])
        _gbnorm = argv[3]
    except:
        print("ERROR: Command line argument(s) is (are) not of correct type.  Please enter an int for argument 1, a float for argument 2, and an int for argument 3. Line 386")
        exit()

    if not len(str(_axis)) == 3: # axis length greater than 3
        print("ERROR: Argument 1 must by a 3 digit number like \'100\'.  Line 390")
        exit()

    axis = [None]*3
    gbnorm = [None]*3
    j = 0
    for i in range(0, len(str(_axis))):
        axis[i] = int(str(_axis)[i])
        try:
            gbnorm[i] = int(_gbnorm[j])
            j = j+1
        except:
            #print(_gbnormal[i:i+2])
            gbnorm[i] = int(_gbnorm[j:j + 2])
            j = j + 2
    _type = defineMisorientation(axis, gbnorm) # Determine the type of misorientation from the axis and gb normal.

    print("These grains have a %s boundary"%(_type))
#------------------------------------------------------------------------------#
#-------------------------------The Actual Calculations------------------------#
#------------------------------------------------------------------------------#
    # First convert to a quaternion

    axis = axis / linalg.norm(axis)
    q = [None]*2
    _z1 = [None]*2
    _x = [None]*2
    _z2 = [None]*2
    at1 = [None]*2
    at2 = [None]*2
    q[0] = [cos((_misorientation/2)/2), sin((_misorientation/2)/2)*axis[0], sin((_misorientation/2)/2)*axis[1], sin((_misorientation/2)/2)*axis[2]]
    q[1] = [cos((-_misorientation/2)/2), sin((-_misorientation/2)/2)*axis[0], sin((-_misorientation/2)/2)*axis[1], sin((-_misorientation/2)/2)*axis[2]]

    # Conversion to Euler angles using the method in the MTEX MATLAB code
    for i in range(0, len(q)):
        at1[i] = atan2(q[i][3],q[i][0])
        at2[i] = atan2(q[i][1],q[i][2])

        _z1[i] = at1[i] - at2[i] + pi / 2.0
        _x[i] = 2*atan2(sqrt(q[i][1]**2 + q[i][2]**2), sqrt(q[i][0]**2 + q[i][3]**2))
        _z2[i] = at1[i] + at2[i] + 3 * pi / 2.0

#---------------------------------------------------------------------------------------------------#
    # Using the Rodrigues Rotation Formula, defined as R = I + sin(theta) * K + (1 - cos(theta))*K^2
    # with K = [0 -k_z, k_y; k_z, 0, -k_x; -k_y, k_x, 0], and the components of
    # k coming from the vector being rotated about.  Theta is specified by the misorientation.
    if useRRF:
        K = array([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        theta = [deg2rad(-_misorientation / 2), deg2rad(_misorientation / 2)]

        for i in range(0, len(theta)):
            R = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]) + sin(theta[i]) * K + (1 - cos(theta[i])) * matMult(K,K)
            orientation_matrix = R

            if not quiet:
                displayMat(orientation_matrix)

            if dispEuler:
                displayAngles(_z1[i], _x[i], _z2[i])

            if save:
                assert i < 2, "ERROR: Too many values for the second Euler angle. Line 442"
                if i == 0:
                    writeMat(orientation_matrix, _z1[i], _x[i], _z2[i], 'P', _axis, _type)
                else:
                    writeMat(orientation_matrix, _z1[i], _x[i], _z2[i], 'Q', _axis, _type)
#----------------------------------------------------------------------------------------------------#
    else:
        for i in range(0,len(_x)):
            orientation_matrix = calcRotMat(_z1[i], _x[i], _z2[i])

            if not quiet:
                displayMat(orientation_matrix)

            if dispEuler:
                displayAngles(_z1[i], _x[i], _z2[i])

            if save:
                assert i < 2, "ERROR: Too many values for the second Euler angle. Line 459"
                if i == 0:
                    writeMat(orientation_matrix, _z1[i], _x[i], _z2[i], 'P', _axis, _type)
                else:
                    writeMat(orientation_matrix, _z1[i], _x[i], _z2[i], 'Q', _axis, _type)
