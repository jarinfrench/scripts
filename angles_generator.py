def genAngleFile(f, maxAngle): #assumes f is a filename, not an open file
    f = open(f,"w")
    for i in range(0, maxAngle+1):
        f.write("%d\n"%(i))


file_100_angles = "AngleFile100.tex"
file_110_angles = "AngleFile110.tex"
file_111_angles = "AngleFile111.tex"

maxAngle = 360

genAngleFile(file_100_angles, maxAngle)
genAngleFile(file_110_angles, maxAngle)
genAngleFile(file_111_angles, maxAngle)
