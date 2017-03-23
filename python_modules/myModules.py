# This file contains all of the functions that I have found useful over time.
# A short description of each function precedes the definition.
from __future__ import division,print_function # makes division and printing easier
from math import sin, cos, pi, sqrt, atan2
from itertools import takewhile,repeat
from numpy import array, linalg

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

# This function converts an axis-angle representation to a quaternion representation
# Assumes that the misorientation is in radians.  Axis must be a 3D vector (list)
def axis2quat(axis, angle):
    axis = axis / linalg.norm(axis)
    return [cos(angle / 2.0), axis[0] * sin(angle / 2.0), axis[1] * sin(angle / 2.0), axis[2] * sin(angle / 2.0)]

# This function converts a quaternion to a matrix representation.
def quat2mat(q):
    e0 = q[0];
    e1 = q[1];
    e2 = q[2];
    e3 = q[3];

    m = array([[e0^2+e1^2-e2^2-e3^2, 2*(e1*e2-e0*e3)     , 2*(e1*e3+e0*e2)],
               [2*(e1*e2+e0*e3)    , e0^2-e1^2+e2^2-e3^2 , 2*(e2*e3-e0*e1)],
               [2*(e1*e3-e0*e2)    , 2*(e2*e3+e0*e1)     , e0^2-e1^2-e2^2+e3^2]])


    return m/(e0^2+e1^2+e2^2+e3^2)

# This function converts a quaternion to Euler angles
def quat2euler(q):
    chi = sqrt((q[0]**2 + q[3]**2)*(q[1]**2 + q[2]**2))

    if chi == 0:
        if q[1] == 0 and q[2] == 0:
            Phi = 0
            phi1 = atan2(2*q[0]*q[3], q[0]**2 - q[3]**2)
            phi2 = 0
        elif q[0] == 0 and q[3] == 0:
            Phi = pi
            phi1 = atan2(s*q[1]*q[2], q[1]**2 - q[2]**2)
            phi2 = 0
        else:
            print("ERROR: Uncalculable Euler angles.")
            exit()
    else:
        Phi = atan2(2 * chi, q[0]**2 + q[3]**2 - q[1]**2 - q[2]**2)
        phi1 = atan2((q[0]*q[2] + q[1]*q[3]) / (2*chi), (q[0]*q[1] - q[2]*q[3]) / (2*chi))
        phi2 = atan2((q[1]*q[3] - q[0]*q[2]) / (2*chi), (q[0]*q[1] + q[2]*q[3]) / (2*chi))
    return phi1, Phi, phi2

# Calculates the Bunge orientation matrix.  Arguments are first z rotation, x rotation, and second z rotation
# Arguments must be in radians!
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

# Calculates the Bunge Orientation matrix using the Rodrigues Rotation Formula
# axis must be a 3D vector (of type list), and the misorientation must in radians!
def calcRotMatRRF(axis, misorientation):
    K = array([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    theta = [-_misorientation / 2, _misorientation / 2]
    R1 = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]) + sin(theta[0]) * K + (1 - cos(theta[0])) * K.dot(K)
    R2 = array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]) + sin(theta[1]) * K + (1 - cos(theta[1])) * K.dot(K)

    return R1, R2
