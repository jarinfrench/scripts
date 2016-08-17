#! /opt/moose/miniconda/bin/python

# This script will generate the rotation matrix for the given misorientation axis
# Input:
#    _rotation_axis     This specifies the axis around which the grains are rotated
#                       rotated. (type: int (100), list ([1, 0, 0]) or string ('100'))
#    _mis_type          This is either 'twist' or 'tilt.' This allows the grain
#                       boundary normal to be determined. (type: string) TODO: determine if a mixed boundary is possible
#    _gbnormal   This specifies the boundary plane normal.  If no option
#                       for this is given, an assumed normal of (010) is assumed
#                       for tilt boundaries, and an assumed normal of (-100) is
#                       assumed for twist boundaries.
# Options:
#    -s --save                           Saves the resultant rotation matrix to
#                                        a database (rotation_matrix_database.tex)
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

# TODO: This needs a major overhaul.  This assumes that the boundary plane is
# perpendicular to the rotation axis (or that the boundary plane normal is
# parallel to the rotation axis), and this is often not a valid assumption.

from __future__ import division,print_function # Automatically divides as floats, and considers print() a function
from sys import argv # for CLI arguments
from numpy import array,linalg # for matrix operations
from os.path import exists # For checking existence of a file
from itertools import takewhile,repeat

# Helper functions
def displayHelp():
    print('''
    This script will generate the rotation matrix for the given misorientation axis
    Input:
        _rotation_axis      This specifies the axis around which the grains are rotated
                            rotated. (type: int (100), list ([1, 0, 0]) or string ('100'))
        _mis_type           This is either 'twist' or 'tilt.' This allows the grain
                            boundary normal to be determined. (type: string)
    Options:
        -s --save           Saves the resultant rotation matrix to a database
                            (rotation_matrix_database.tex) with the accompanying
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

def displayMat(m): # Displays the matrix
    for i in range(0, len(m)):
        for j in range(0,len(m[i])):
            print("%2.6f\t"%(m[i][j]),end="")
        print("\n")
    print("\n")
    return

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

def writeMat(m, axis, _type ): # Write the matrix and angles to a file
    tex_filename = "rotation_matrix_database.tex"
    var_name = "rot%s%s"%(axis,_type)
    if not exists(tex_filename):
        tex_file = open(tex_filename, "a")
        tex_file.write("%Database for rotation matrices for specified misorientation axes\n")
        tex_file.write("%----------------------------------------------------------------\n")
        tex_file.write("%Rotation Matrix                                    Axis     Type\n")
        tex_file.write("%s=[%2.4f  %2.4f  %2.4f                 %%%s      %s\n"%(var_name, m[0][0], m[0][1], m[0][2], axis, _type.title()))
        tex_file.write("%2.4f  %2.4f  %2.4f\n"%(m[1][0], m[1][1], m[1][2]))
        tex_file.write("%2.4f  %2.4f  %2.4f];\n"%(m[2][0], m[2][1], m[2][2]))
        tex_file.write("%----------------------------------------------------------------\n")
        tex_file.close()
    else:
        # Check for already written
        numlines = rawbigcount(tex_filename)
        if numlines <= 4:
            tex_file = open(tex_filename, "a")
            tex_file.write("%s=[%2.4f  %2.4f  %2.4f                 %%%s      %s\n"%(var_name, m[0][0], m[0][1], m[0][2], axis, _type.title()))
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
                elif len(data) != 5:
                    continue
                else:
                    if data[3] == (str('%') + axis) and data[4] == _type.title():
                        unique = False
                    else:
                        unique = True
            if unique:
                tex_file = open(tex_filename, "a")
                tex_file.write("%s=[%2.4f  %2.4f  %2.4f                 %%%s      %s\n"%(var_name, m[0][0], m[0][1], m[0][2], axis, _type.title()))
                tex_file.write("%2.4f  %2.4f  %2.4f\n"%(m[1][0], m[1][1], m[1][2]))
                tex_file.write("%2.4f  %2.4f  %2.4f];\n"%(m[2][0], m[2][1], m[2][2]))
                tex_file.write("%----------------------------------------------------------------\n")
                tex_file.close()


def rawbigcount(filename): # Counts the number of lines in a file
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen if buf )

# Check what we were given...
if "--help" in argv: # Help info
    displayHelp()
    exit()

save, argv = check4Save(argv)
quiet, argv = check4Quiet(argv) # Checks for suppressing output

if len(argv) != 4: # if not all three values given
    if len(argv) != 3: # Check for _rotation_axis and _mis_type, otherwise
        print("ERROR: Incorrect number of command line arguments.")
        displayHelp()
        exit()
    else: # len(argv) == 3
        script, _rotation_axis, _mis_type = argv
        assert _mis_type in {'twist', 'tilt'}, "ERROR: Misorientation type not recognized"
        if _mis_type == 'twist':
            _gbnormal = '010'
        else:
            _gbnormal = '-100'
else: # len(argv) == 4
    script, _rotation_axis, _mis_type, _gbnormal = argv
    if not type(_gbnormal) in {int, list, str}:
        print("ERROR: Grain boundary normal type is incorrect.  Please enter an int, a list, or a string")
        print("You entered %s with type %s"%(str(_gbnormal), type(_gbnormal)))
        exit()
    else:
        if type(_gbnormal) == str:
            pass # Do nothing
        elif type(_gbnormal) == int:
            _gbnormal = '0' + '0' + '0' + str(_gbnormal)
            _gbnormal =_gbnormal[-3:] # gets the last three characters
            assert _gbnormal != '000', "ERROR: invalid boundary normal."
        else type(_gbnormal) == list:
            _gbnormal = ''.join(str(i) for i in _gbnormal)

        assert(len(_rotation_axis) == 3), "ERROR: Something went wrong converting _gbnormal into a string."

if not type(_rotation_axis) in {int, list, str}:
    print("ERROR: Grain boundary rotation axis type is incorrect.  Please enter an int, a list, or a string")
    print("You entered %s with type %s"%(str(_rotation_axis), type(_rotation_axis)))
    exit()
else: # Convert anything besides a string into a string
    if type(_rotation_axis) == str:
        pass # Do nothing
    elif type(_rotation_axis) == int:
        _rotation_axis = '0' + '0' + '0' + str(_rotation_axis)
        _rotation_axis = _rotation_axis[-3:] # Get the last three characters
        assert _rotation_axis != '000', "ERROR: invalid rotation axis."
    else type(_rotation_axis) == list:
        _rotation_axis = ''.join(str(i) for i in _rotation_axis)

    assert(len(_rotation_axis) == 3), "ERROR: Something went wrong converting _rotation_axis into a string."

# Now that the input is taken care of, do the work
axis = array([[1, 0, 0]]) # This is the axis that we rotate the grain boundary normal to
#if _mis_type == 'twist':
#    gbnormal = array([[None]*3])
#    for i in range(0, len(_rotation_axis)):
#        gbnormal[0][i] = int(_rotation_axis[i])
#elif _mis_type == 'tilt':
#    if _rotation_axis == '100':
#        gbnormal = array([[0, 1, 0]])
#    elif _rotation_axis == '110':
#        gbnormal = array([[1, -1, 0]])
#    elif _rotation_axis == '111':
#        gbnormal = array([[1, -1, 0]])
#    else:
#        print("ERROR: Non-standard axis used.")
#        exit()
#else:
#    print("ERROR: Misorientation type not recognized")
#    exit()
gbnormal == array([[None]*3])
for i in len(_gbnormal):
    gbnormal[0][i] = int(_gbnormal[i])

print(gbnormal)

# So much work...
gbnormal = gbnormal / linalg.norm(gbnormal) # Normalize the gbnormal vector.
axis = axis / linalg.norm(axis) # Just to be sure...
rotation_matrix = rotVec1ToVec2(gbnormal, axis)

if not quiet:
    displayMat(rotation_matrix)
if save:
    writeMat(rotation_matrix, _rotation_axis, _mis_type)
