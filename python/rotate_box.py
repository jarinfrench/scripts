from visual import *
import argparse

parser = argparse.ArgumentParser(usage = 'python %(prog)s x y z', description = "A simple python script to visualize rotation of a cube about a specific axis.")
parser.add_argument('axis', type = int, nargs = 3, help = "The desired rotation axis")
args = parser.parse_args()

# axis = input("Please enter the axis to be rotated about: ")

redbox = box(pos=vector(0,0,0), size=(2,2,2), color=color.red, opacity=0.2)
pointer = arrow(pos = (0,0,0), axis=(args.axis[0],args.axis[1],args.axis[2]), length=2, color=color.white)
for i in range(360):
    rate(30)
    redbox.rotate(angle=radians(1), axis=(args.axis[0],args.axis[1],args.axis[2]), origin=(0,0,0))
    pointer.rotate(angle=radians(1), axis=(args.axis[0],args.axis[1],args.axis[2]), origin=(0,0,0))
