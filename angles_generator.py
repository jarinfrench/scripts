
def genAngleFile(f, maxAngle): #assumes f is a filename, not an open file
    f = open(f,"w")
    for i in range(0, maxAngle+1):
        f.write("%d\n"%(i))


file_100_tilt = "Tilt100Angles.tex"
file_110_tilt = "Tilt110Angles.tex"
file_111_tilt = "Tilt111Angles.tex"
file_100_twist = "Twist100Angles.tex"
file_110_twist = "Twist110Angles.tex"
file_111_twist = "Twist111Angles.tex"

max100tilt = 90 # Also max 110 twist
max110tilt = 180
max111tilt = 60 # Also max 111 twist
max100twist = 45

genAngleFile(file_100_tilt, max100tilt)
genAngleFile(file_110_tilt, max110tilt)
genAngleFile(file_111_tilt, max111tilt)
genAngleFile(file_100_twist, max100twist)
genAngleFile(file_110_twist, max100tilt)
genAngleFile(file_111_twist, max111tilt)
