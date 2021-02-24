# This file contains all of the functions that I have found useful over time.
# A short description of each function precedes the definition.
# from __future__ import division,print_function # makes division and printing easier
from math import sin, cos, pi, sqrt, atan2
from itertools import takewhile,repeat
import numpy as np
import sys, os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

plot_style = ['bo', 'go', 'ro', 'co', 'mo', 'yo', 'ko',
              'bv', 'gv', 'rv', 'cv', 'mv', 'yv', 'kv',
              'b^', 'g^', 'r^', 'c^', 'm^', 'y^', 'k^',
              'b<', 'g<', 'r<', 'c<', 'm<', 'y<', 'k<',
              'b>', 'g>', 'r>', 'c>', 'm>', 'y>', 'k>',
              'b1', 'g1', 'r1', 'c1', 'm1', 'y1', 'k1',
              'b2', 'g2', 'r2', 'c2', 'm2', 'y2', 'k2',
              'b3', 'g3', 'r3', 'c3', 'm3', 'y3', 'k3',
              'b4', 'g4', 'r4', 'c4', 'm4', 'y4', 'k4',
              'b8', 'g8', 'r8', 'c8', 'm8', 'y8', 'k8',
              'bs', 'gs', 'rs', 'cs', 'ms', 'ys', 'ks',
              'bp', 'gp', 'rp', 'cp', 'mp', 'yp', 'kp',
              'bP', 'gP', 'rP', 'cP', 'mP', 'yP', 'kP',
              'b*', 'g*', 'r*', 'c*', 'm*', 'y*', 'k*',
              'bh', 'gh', 'rh', 'ch', 'mh', 'yh', 'kh',
              'bH', 'gH', 'rH', 'cH', 'mH', 'yH', 'kH',
              'b+', 'g+', 'r+', 'c+', 'm+', 'y+', 'k+',
              'bx', 'gx', 'rx', 'cx', 'mx', 'yx', 'kx',
              'bX', 'gX', 'rX', 'cX', 'mX', 'yX', 'kX',
              'bD', 'gD', 'rD', 'cD', 'mD', 'yD', 'kD',
              'bd', 'gd', 'rd', 'cd', 'md', 'yd', 'kd',
              'b|', 'g|', 'r|', 'c|', 'm|', 'y|', 'k|',
              'b_', 'g_', 'r_', 'c_', 'm_', 'y_', 'k_']

# This function converts degrees to radians.  Argument is assumed to be in degrees
def deg2rad(x):
    return x * pi / 180.0

# This function converts radians to degrees.  Argument is assumed to be in radians
def rad2deg(x):
    return x * 180.0 / pi

# This function displays a matrix in an easy-to-see format.  This will display
# a matrix of any size and dimension
def displayMat(m):
    for i in range(len(m)):
        for j in range(len(m[i])):
            print(f"{(m[i][j]):.6f}",end="")
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
    axis = axis / np.linalg.norm(axis)
    return [cos(angle / 2.0), axis[0] * sin(angle / 2.0), axis[1] * sin(angle / 2.0), axis[2] * sin(angle / 2.0)]

# This function converts a quaternion to a matrix representation.
def quat2mat(q):
    e0 = q[0];
    e1 = q[1];
    e2 = q[2];
    e3 = q[3];

    m = np.array([[e0^2+e1^2-e2^2-e3^2, 2*(e1*e2-e0*e3)     , 2*(e1*e3+e0*e2)],
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
            phi1 = atan2(2*q[1]*q[2], q[1]**2 - q[2]**2)
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
    K = np.array([[0, -axis[2], axis[1]],[axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    theta = [-_misorientation / 2, _misorientation / 2]
    R1 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]) + sin(theta[0]) * K + (1 - cos(theta[0])) * K.dot(K)
    R2 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],[0.0, 0.0, 1.0]]) + sin(theta[1]) * K + (1 - cos(theta[1])) * K.dot(K)

    return R1, R2

# Checks that the three vectors passed in are mutually orthogonal
def checkOrthogonality(x, y, z):
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    if (x.dot(y) < 1e-15 and x.dot(z) < 1e-15 and y.dot(x) < 1e-15):
        return True
    else:
        return False

# Checks an input value to determine if it is a scalar or not
def checkScalarInput(input):
    try:
        x=len(input)
        return False
    except:
        # if we are here, len(input) did not work
        return True # Therefore input is a scalar

def normalize(vector):
    return vector / np.linalg.norm(vector)

def promptForContinue():
    while True:
        if sys.version_info.major == 2:
            cont = raw_input("Would you like to continue? (y|n): ").lower()
        elif sys.version_info.major >= 3:
            cont = input("Would you like to continue? (y|n): ").lower()
        if len(cont) == 1 and cont in ['y', 'n']:
            break
        print("Please enter either Y or N")
    if cont == 'y':
        return True
    else:
        return False

def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

def verify_new_file(filebase):
    n = 1
    file = filebase
    while os.path.isfile(file):
        file = "{base}_{num}{ext}".format(base = os.path.splitext(filebase)[0], num = n, ext = os.path.splitext(filebase)[1])
        n += 1

    return file

# flattens a nested list of lists (of lists of lists...)
def flatten(data_list):
    if data_list == []:
        return data_list
    if isinstance(data_list[0], list):
        return flatten(data_list[0]) + flatten(data_list[1:])
    elif isinstance(data_list, np.ndarray):
        return np.ndarray.flatten(data_list)
    return data_list[:1] + flatten(data_list[1:])

# simple way to return either the length of a list, or the "length" of a single number
def depth(data):
    try:
        return len(data)
    except:
        return 1

# counts the number of times each item appears in a list
def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs) > 1)

# searches a list for a element, and returns the index of the element
def find(lst, elem):
    for i, x in enumerate(lst):
        if x == elem:
            return i
    return None

# splits the passed in list 'seq' into chunks of size 'size'
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))
