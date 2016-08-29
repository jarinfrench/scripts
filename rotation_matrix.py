#! /opt/moose/miniconda/bin/python

# This script will generate the rotation matrix for the given misorientation axis
# Input:
#    _rotation_axis     This specifies the axis around which the grains are rotated
#                       rotated. (type: int (100), list ([1, 0, 0]) or string ('100'))
#
#                       boundary normal to be determined. (type: string) TODO: determine if a mixed boundary is possible
#    _gbnormal          This specifies the boundary plane normal.  If no option
#                       for this is given, an assumed normal of (010) is used
#                       for tilt boundaries, and an assumed normal of (-100) is
#                       used for twist boundaries.
# Options:
#    -s --save                           Saves the resultant rotation matrix to
#                                        a database (rotation_matrix_database.m)
#                                        with the accompanying rotation axis
#                                        and misorientation type.
#
#    -q --quiet                          Suppresses output of the rotation matrix
#
#    --help                              Display this help info
#
# Output:
# The output displayed will be the resultant rotation matrix for the given
# misorientation.

from __future__ import division,print_function # Automatically divides as floats, and considers print() a function
from sys import argv # for CLI arguments
from numpy import array,linalg # for matrix operations
from os.path import exists # For checking existence of a file
from myModules import * # For using my user-defined functions

# Helper functions
def displayHelp():
    print('''
    This script will generate the rotation matrix for the given misorientation axis
    Input:
        _rotation_axis      This specifies the axis around which the grains are
                            rotated. (type: int (100), list ([1, 0, 0]) or string ('100'))

        _gbnormal           This specifies the boundary plane normal.  If no option
                            for this is given, an assumed normal of (010) is used
                            for tilt boundaries, and an assumed normal of (-100) is
                            used for twist boundaries. (type: int (100), list
                            ([1, 0, 0]) or string ('100'))
    Options:
        -s --save           Saves the resultant rotation matrix to a database
                            (rotation_matrix_database.m) with the accompanying
                            rotation axis and misorientation type.

        -q --quiet          Suppresses output of the rotation matrix

        --help              Display this help info.
    Output:
    The output displayed will be the resultant rotation matrix for the given
    misorientation.
    ''')
def rotVecToZ(vec): # Creates the rotation matrix to rotate vec to the z direction
    # REALLY make sure vec is normalized
    vec = vec / linalg.norm(vec)

    # Initialize our vectors
    v1 = array([[0.,0.,0.]])
    v0 = array([[0.,0.,0.]])

    # Temp vector that gives a prototype of v1 by looking at the smallest component of vec
    w = array([[abs(vec[0][0]), abs(vec[0][1]), abs(vec[0][2])]])
    if ( (w[0][2] >= w[0][1] and w[0][1] >= w[0][0]) or (w[0][1] >= w[0][2] and w[0][2] >= w[0][0]) ):
        v1[0][0] = 1.0
    elif ( (w[0][2] >= w[0][0] and w[0][0] >= w[0][1]) or (w[0][0] >= w[0][2] and w[0][2] >= w[0][1]) ):
        v1[0][1] = 1.0
    else:
        v1[0][2] = 1.0

    # Gram-Schmidt method to find v1
    v1 = v1 - ((v1.dot(vec.T))*vec)
    v1 = v1 / linalg.norm(v1)

    # v0 = v1 x vec
    v0[0][0] = v1[0][1]*vec[0][2] - v1[0][2]*vec[0][1]
    v0[0][1] = v1[0][2]*vec[0][0] - v1[0][0]*vec[0][2]
    v0[0][2] = v1[0][0]*vec[0][1] - v1[0][1]*vec[0][0]

    # Rotation matrix is just:
    rot = array([[v0[0][0], v0[0][1], v0[0][2]],
                 [v1[0][0], v1[0][1], v1[0][2]],
                 [vec[0][0],vec[0][1],vec[0][2]]])
    return rot

def rotVec1ToVec2(vec1, vec2):
    rot1_to_z = rotVecToZ(vec1)
    rot2_to_z = rotVecToZ(vec2)
    return (rot2_to_z.T).dot(rot1_to_z)

def writeMat(m, _axis, gbnormal): # Write the matrix and angles to a file
    tex_filename = "rotation_matrix_database.m"
    normName = _axis + '_' + gbnormal
    var_name = "rot%snorm"%(normName)
    if not exists(tex_filename):
        tex_file = open(tex_filename, "a")
        tex_file.write("%Database for rotation matrices for specified misorientation axes\n")
        tex_file.write("%----------------------------------------------------------------\n")
        tex_file.write("%Rotation Matrix\n")
        tex_file.write("%s=[%2.4f  %2.4f  %2.4f\n"%(var_name, m[0][0], m[0][1], m[0][2]))
        tex_file.write("%2.4f  %2.4f  %2.4f\n"%(m[1][0], m[1][1], m[1][2]))
        tex_file.write("%2.4f  %2.4f  %2.4f];\n"%(m[2][0], m[2][1], m[2][2]))
        tex_file.write("%----------------------------------------------------------------\n")
        tex_file.close()
    else:
        # Check for already written
        numlines = countFileLines(tex_filename)
        if numlines <= 4:
            tex_file = open(tex_filename, "a")
            tex_file.write("%s=[%2.4f  %2.4f  %2.4f\n"%(var_name, m[0][0], m[0][1], m[0][2]))
            tex_file.write("%2.4f  %2.4f  %2.4f\n"%(m[1][0], m[1][1], m[1][2]))
            tex_file.write("%2.4f  %2.4f  %2.4f];\n"%(m[2][0], m[2][1], m[2][2]))
            tex_file.write("%----------------------------------------------------------------\n")
            tex_file.close()
        else:
            f = open(tex_filename,"r")
            while True:
                data = f.readline().split()
                if not data:
                    break
                elif len(data[0]) > 14:
                    if not data[0][14] == '=':
                        continue
                    else:
                        if data[0][0:10] == var_name:
                            unique = False
                        else:
                            unique = True
            if unique:
                tex_file = open(tex_filename, "a")
                tex_file.write("%s=[%2.4f  %2.4f  %2.4f\n"%(var_name, m[0][0], m[0][1], m[0][2]))
                tex_file.write("%2.4f  %2.4f  %2.4f\n"%(m[1][0], m[1][1], m[1][2]))
                tex_file.write("%2.4f  %2.4f  %2.4f];\n"%(m[2][0], m[2][1], m[2][2]))
                tex_file.write("%----------------------------------------------------------------\n")
                tex_file.close()

# Check what we were given...
if "--help" in argv: # Help info
    displayHelp()
    exit()

save, argv = check4Save(argv)
quiet, argv = check4Quiet(argv) # Checks for suppressing output

if len(argv) != 3: # if not both values given
    print("ERROR: Incorrect number of command line arguments. Line 151")
    displayHelp()
    exit()
else: # len(argv) == 3
    script, _rotation_axis, _gbnormal = argv

if not type(_gbnormal) in {int, str}:
    print("ERROR: Grain boundary normal type is incorrect.  Please enter an int or a string. Line 157")
    print("You entered %s with type %s"%(str(_gbnormal), type(_gbnormal)))
    exit()
else:
    if type(_gbnormal) == int:
        _gbnormal = '0' + '0' + '0' + str(_gbnormal)
        _gbnormal =_gbnormal[-3:] # gets the last three characters
        assert _gbnormal != '000', "ERROR: invalid boundary normal. Line 173"

    assert(len(_rotation_axis) == 3), "ERROR: Something went wrong converting _gbnormal into a string. Line 177"

if not type(_rotation_axis) in {int, str}:
    print("ERROR: Grain boundary rotation axis type is incorrect.  Please enter an int or a string. Line 180")
    print("You entered %s with type %s"%(str(_rotation_axis), type(_rotation_axis)))
    exit()
else: # Convert anything besides a string into a string
    if type(_rotation_axis) == int:
        _rotation_axis = '0' + '0' + '0' + str(_rotation_axis)
        _rotation_axis = _rotation_axis[-3:] # Get the last three characters
        assert _rotation_axis != '000', "ERROR: invalid rotation axis. Line 189"

    assert(len(_rotation_axis) == 3), "ERROR: Something went wrong converting _rotation_axis into a string. Line 193"

# Now that the input is taken care of, do the work
axis = array([[1, 0, 0]]) # This is the axis that we rotate the grain boundary normal to

# This part converts _gbnormal to an array for use in the rotation functions
gbnormal = array([[None]*3])
j = 0
indices = []
for i in range(0,len(gbnormal[0])):
    try:
        gbnormal[0][i] = int(_gbnormal[j])
        j = j+1
    except:
        #print(_gbnormal[i:i+2])
        gbnormal[0][i] = int(_gbnormal[j:j + 2])
        j = j + 2
        indices.append(i)


# So much work...
gbnorm = ''
for i in range(0,len(gbnormal[0])):
    if gbnormal[0][i] < 0:
        gbnorm = gbnorm + str(abs(gbnormal[0][i])) + 'bar'
    else:
        gbnorm = gbnorm + str(gbnormal[0][i])
gbnormal = gbnormal / linalg.norm(gbnormal) # Normalize the gbnormal vector.
axis = axis / linalg.norm(axis) # Just to be sure...
rotation_matrix = rotVec1ToVec2(gbnormal, axis)



if not quiet:
    displayMat(rotation_matrix)
if save:
    writeMat(rotation_matrix, _rotation_axis, gbnorm)
