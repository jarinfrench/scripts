from visual import *

axis = input("Please enter the axis to be rotated about: ")

redbox = box(pos=vector(0,0,0), size=(2,2,2), color=color.red, opacity=0.2)
pointer = arrow(pos = (0,0,0), axis=(int(axis[0]),int(axis[1]),int(axis[2])), length=2, color=color.white)
for i in range(360):
    rate(30)
    redbox.rotate(angle=radians(1), axis=(int(axis[0]),int(axis[1]),int(axis[2])), origin=(0,0,0))
    pointer.rotate(angle=radians(1), axis=(int(axis[0]),int(axis[1]),int(axis[2])), origin=(0,0,0))
