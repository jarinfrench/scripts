def genAngleFile(f, maxAngle): #assumes f is a filename, not an open file
    f = open(f,"w")
    for i in range(0, maxAngle+1):
        f.write("%d\n"%(i)) # Writes each i to a new line in file f


# A bunch of file names.  Uses the convention required in genOrientationMatrices.sh
file_100_angles = "AngleFile100.tex"
file_110_angles = "AngleFile110.tex"
file_111_angles = "AngleFile111.tex"

# This can be chosen specifically for each set based on symmetries
# i.e. for 100 tilt, max angle could be 90, 100 twist could be 45,
# 110 tilt could be 180, 110 twist could be 90, and both 111 sets could be 60.
# This just creates all of the angles from 0 to 360, giving a full period of
# sorts to describe each boundary.
maxAngle = 360

genAngleFile(file_100_angles, maxAngle)
genAngleFile(file_110_angles, maxAngle)
genAngleFile(file_111_angles, maxAngle)
