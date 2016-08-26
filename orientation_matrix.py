#! /opt/moose/miniconda/bin/python

# This script will calculate the orientation matrices for any given misorientation
# for any of the high-symmetry axes.
# Arguments:
#
#   _axis: axis of orientation (needs to be either 100, 110, or 111) (type: int or list)
#   _misorientation: The angle of misorientation (type: float)
#   _type: Type of misorientation (twist or tilt) (type: string)
#            --------OR--------
#    (with option -e or --euler-angles)
#   _z1: The first rotation angle  (Z ) (type: float)
#   _x:  The second rotation angle (X') (type: float)
#   _z2: The third rotation angle  (Z") (type: float)
#
# If the option -e or --euler-angles is entered, the calculation skips to simply
# output the orientation matrices.  Otherwise, the euler angles are calculated from
# the axis, orientation, and misorientation type, and then the orientation matrix is
# created through the use of the Rodrigues Rotation Formula, which is:
# R = I + sin(theta) * K + (1 - cos(theta))*K^2
# where I is the identity matrix, and K is the skew-symmetric matrix formed by
# the axis of rotation:
# K = 0  -kz  ky
#     kz  0  -kx
#    -ky  kx  0
# Where the vector k is the unit vector defining the axis of rotation.
#
# Options:
# -e --euler-angles <_z1> <_x> <_z2>    Returns the  Bunge orientation matrix
#                                       based on the euler angles provided.
# -m --misorientation <_axis> <_misorientation> <_type>
#                                       Returns the Bunge orientation matrix after
#                                       calculating the Euler angles from the
#                                       information provided (default behavior)
#
# -s --save                             Saves the resultant orientation matrix to
#                                       a database (orientation_matrix_database.tex)
#                                       with the accompanying Euler angles.
#
# -q --quiet                            Suppresses output of the orientation matrices
#
# --help                                Display this help info
#
# Output:
# For an Euler angle set, the ouput is simply its orientation matrix.
# For the misorientations, the first matrix is the 'P' orientation matrix, and
# the second matrix is the 'Q' orientation matrix (see Bulatov et al. Acta Mater
# 65 (2014) 161-175).

from __future__ import division,print_function # To avoid numerical problems with division, and for ease of printing
from sys import argv # for CLI arguments
from math import cos,sin,pi # Trig functions
from numbers import Number # For checking input parameters
from os.path import exists # For checking existence of a file
from numpy import array, linalg

orientation_matrix = []

# Helper functions
def displayHelp():
    print('''
    This script will calculate the orientation matrices for any given misorientation
    for any of the high-symmetry axes.
    Arguments:

       _axis: axis of orientation (needs to be either 100, 110, or 111) (type: int or list)
       _misorientation: The angle of misorientation (type: float)
       _type: Type of misorientation (twist or tilt) (type: string)
                  --------OR--------
        (with option -e or --euler-angles)
       _z1: The first rotation angle  (Z ) (type: float)
       _x:  The second rotation angle (X') (type: float)
       _z2: The third rotation angle  (Z") (type: float)

    If the option -e or --euler-angles is entered, the calculation skips to simply
    output the orientation matrices.  Otherwise, the euler angles are calculated from
    the axis, orientation, and misorientation type, and then the orientation matrix is
    created through the use of the Rodrigues Rotation Formula, which is:
    R = I + sin(theta) * K + (1 - cos(theta))*K^2
    where I is the identity matrix, and K is the skew-symmetric matrix formed by
    the axis of rotation:
    K = 0  -kz  ky
        kz  0  -kx
       -ky  kx  0
    Where the vector k is the unit vector defining the axis of rotation.

    Options:
    -e --euler-angles <_z1> <_x> <_z2>  Returns the  Bunge orientation matrix
                                        based on the euler angles provided.
    -m --misorientation <_axis> <_misorientation> <_type>
                                        Returns the Bunge orientation matrix after
                                        calculating the Euler angles from the
                                        information provided (default behavior)

    -s --save                           Saves the resultant orientation matrix to
                                        a database (orientation_matrix_database.tex)
                                        with the accompanying Euler angles.

    -q --quiet                          Suppresses output of the orientation matrices

    --help                              Display this help info

    Output:
    For an Euler angle set, the ouput is simply its orientation matrix.
    For the misorientations, the first matrix is the 'P' orientation matrix, and
    the second matrix is the 'Q' orientation matrix (see Bulatov et al., Acta Mater
    65 (2014) 161-175).
    ''')
    return

def calcRotMat(_z1,_x,_z2): # Calculates the Bunge orientation matrix.  Arguments are first z rotation, x rotation, and second z rotation
    c1 = cos(_z1)
    c2 = cos(_x)
    c3 = cos(_z2)
    s1 = sin(_z1)
    s2 = sin(_x)
    s3 = sin(_z2)

    rot_mat = array([[ c1*c3 - c2*s1*s3,  c3*s1 + c1*c2*s3,  s2*s3],
                     [-c1*s3 - c2*c3*s1,  c1*c2*c3 - s1*s3,  c3*s2],
                     [ s1*s2           , -c1*s2           ,  c2  ]])
    return rot_mat

def deg2rad(x): # Converts degrees to radians.  Argument is in degrees
    return x * pi / 180.0

def displayMat(m): # Displays the matrix
    for i in range(0, len(m)):
        for j in range(0,len(m[i])):
            print("%2.6f\t"%(m[i][j]),end="")
        print("\n")
    print("\n")
    return

def matMult(m1,m2): # Multiplies two matrices together
    result = m1.dot(m2)
    return result

def check4Save(args): # Check the args for a save command
    if "-s" in args or "--save" in args: # We want to save the file
        try:
            index = args.index("-s")
        except:
            index = args.index("--save")
        del args[index]
        return True, args
    else:
        return False, args

def check4Quiet(args): # Check the args for a quiet command
    if "-q" in args or "--quiet" in args: # We don't want terminal output
        try:
            index = args.index("-q")
        except:
            index = args.index("--quiet")
        del args[index]
        return True, args
    else:
        return False, args

def writeMat(m, _z1, _x, _z2, grain, axis, _type): # Write the matrix and angles to a file
    # This is to avoid issues with duplicates
    if _z1 == 0:
        _z1 = abs(_z1)
    if _x == 0:
        _x = abs(_x)
    if _z2 == 0:
        _z2 = abs(_z2)

    lastVal = 1
    #tex_filename = "orientation_matrix_database_ints.tex"
    tex_filename = "orientation_matrix_database_ints_rrf.tex"
    #var_name = "%s%d%s"%(grain, axis, _type)
    var_name = "%s%d%srrf"%(grain, axis, _type)
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
                                #lastVal = int(data[0][14:17]) - 1
                                lastVal = int(data[0][17:20]) - 1
                            else:
                                #lastVal = int(data[0][14:17])
                                lastVal = int(data[0][17:20])
                        except:
                            if data[0][0] == 'P': # Handles anything 2 digits long
                                #lastVal = int(data[0][14:16]) - 1
                                lastVal = int(data[0][17:19]) - 1
                            else: # data[0][0] == 'Q'
                                #lastVal = int(data[0][14:16])
                                lastVal = int(data[0][17:19])
                    except:
                        if data[0][0] == 'P': # One digit case
                            #lastVal = int(data[0][14]) - 1
                            lastVal = int(data[0][17]) - 1
                        else: # data[0][0] == 'Q'
                            #lastVal = int(data[0][14])
                            lastVal = int(data[0][17])
                elif data[0][5] == 'i': # P or Q1xxtilt(i)
                    try:
                        try:
                            if data[0][0] == 'P':
                                #lastVal = int(data[0][13:16]) - 1
                                lastVal = int(data[0][16:19]) - 1
                            else:
                                #lastVal = int(data[0][13:16])
                                lastVal = int(data[0][16:19])
                        except:
                            if data[0][0] == 'P':
                                #lastVal = int(data[0][13:15]) - 1
                                lastVal = int(data[0][16:18]) - 1
                            else: # data[0][0] == 'Q'
                                #lastVal = int(data[0][13:15])
                                lastVal = int(data[0][16:18])
                    except:
                        if data[0][0] == 'P':
                            #lastVal = int(data[0][13]) - 1
                            lastVal = int(data[0][16]) - 1
                        else: # data[0][0] == 'Q'
                            #lastVal = int(data[0][13])
                            lastVal = int(data[0][16])
                else:
                    print("Error: Unknown last index")
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

save, argv = check4Save(argv) # Save the file?  Delete the save argument
quiet, argv = check4Quiet(argv) # Checks for suppressing output
# Error checking for input arguments
if "--help" in argv: # Help info
    displayHelp()
    exit()
if "-f" in argv or "--file" in argv: #input arguments come from file
    try:
        try:
            index = argv.index("-f")
        except:
            index = argv.index("--file")
    except:
        print("ERROR: Unable to find filename")
        exit()
    filename = argv[index + 1]
    try:
        f1 = open(filename, 'r')
    except:
        print("ERROR: Unable to read file", filename)

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
        print("ERROR: Unable to read Euler angles")
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

    exit()
else:
    if len(argv) < 4:
        print("ERROR: Not enough command line arguments.")
        print("Input either an axis, misorientation, and misorientation type, or a ZXZ Euler angle set with the option -e or --euler-angles.")
        displayHelp()
        exit()
    _axis = int(argv[1])
    _misorientation = float(argv[2])
    _type = argv[3]

    # Additional input error checking
    if len(str(_axis)) > 3: # axis length greater than 3
        print("ERROR: Argument 1 must by a 3 digit number like \'100\'")
        exit()

    if not _axis in {100, 110, 111}: # TODO: determine from ANY axis
        print("ERROR: Unrecognized axis.  Please enter \'100\', \'110\', or \'111\'")
        exit()

    try:
        _misorientation = float(_misorientation)
    except:
        print("ERROR: Argument 2 must be a number")
        exit()

    _type = _type.lower() # makes the string all lower case
    if not _type in ['twist', 'tilt']: # type is not twist or tilt
        print("ERROR: Argument 3 must be either \'twist\' or \'tilt\'")
        exit()

#------------------------------------------------------------------------------#
#-------------------------------The Actual Calculations------------------------#
#------------------------------------------------------------------------------#
    _x = [None]*2
    # First, check to see if twist or tilt
    assert _axis in {100, 110, 111}, "ERROR: Unrecognized axis."
    if _type == 'tilt': # Assuming symmetric tilt ONLY
        if _axis == 100:
            _z1   = 0.00
            _x[0] = _misorientation / 2.0
            _x[1] = - _misorientation / 2.0
            _z2   = 0.00

        elif _axis == 110:
            _z1   = 45.00
            _x[0] = _misorientation / 2.0
            _x[1] = - _misorientation / 2.0
            _z2   = 0.00

        elif _axis == 111:
            _z1   = 290.104
            _x[0] = 37.9381 + _misorientation / 2.0
            _x[1] = 37.9381 - _misorientation / 2.0
            _z2   = 110.104

    elif _type == 'twist': # NOTE: for twist misorientations, the second grain _does_not_move_!!!
        if _axis == 100:
            _z1   = 0.00
            _x[0] = _misorientation
            _x[1] = 0.00
            _z2   = 0.00

        elif _axis == 110:
            _z1   = 45.0
            _x[0] = _misorientation
            _x[1] = 0.00
            _z2   = 0.0

        elif _axis == 111: # Except in the case of the <111> twist.
            # For this special case we assume the misorientation is the same as
            # for the tilt.  Once the orientation matrix is calculated, we then
            # rotate THAT matrix by the rotation necessary to get from the <111>
            # axis to the <100> axis.
            # EDIT 23 Aug 2016: It now follows the same pattern, but uses a different formula to calculate the orientation matrix
            _z1   = 290.104
            _x[0] = 37.9381 + _misorientation
            _x[1] = 0.0
            _z2   = 110.104

    else:
        print("ERROR: Unrecognized boundary type")
        exit()

    _z1   = deg2rad(_z1)
    _x[0] = deg2rad(_x[0])
    _x[1] = deg2rad(_x[1])
    _z2   = deg2rad(_z2)
#    for i in range(0,len(_x)):
#        orientation_matrix = calcRotMat(_z1, _x[i], _z2)
#---------------------------------------------------------------------------------------------------#
    # Using the Rodrigues Rotation Formula, defined as R = I + sin(theta) * K + (1 - cos(theta))*K^2
    # with K = [0 -k_z, k_y; k_z, 0, -k_x; -k_y, k_x, 0], and the components of
    # k coming from the vector being rotated about.  Theta is specified by the misorientation.
    k = [None]*3
    for i in range(0,len(str(_axis))):
        k[i] = float(str(_axis)[i])
    k = k / linalg.norm(k)
    K = array([[0, -k[2], k[1]],[k[2], 0, -k[0]], [-k[1], k[0], 0]])
    if _type == 'twist':
        theta = [deg2rad(_misorientation), 0]
    else:
        theta = [deg2rad(_misorientation / 2), deg2rad(-_misorientation / 2)]

    for i in range(0, len(theta)):
        R = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]) + sin(theta[i]) * K + (1 - cos(theta[i])) * matMult(K,K)
        orientation_matrix = R
#----------------------------------------------------------------------------------------------------#

        if not quiet:
            displayMat(orientation_matrix)

        if save:
            assert i < 2, "ERROR: Too many values for the second Euler angle."
            if i == 0:
                writeMat(orientation_matrix, _z1, _x[i], _z2, 'P', _axis, _type)
            else:
                writeMat(orientation_matrix, _z1, _x[i], _z2, 'Q', _axis, _type)
