# Author: Axel Seoane, 2018
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import time

#turns a 1D array into a unit vector
def makeUnitVector(rot_axis):
    length = np.sqrt(np.dot(rot_axis,rot_axis))
    for i in range(len(rot_axis)):
        rot_axis[i] = float(rot_axis[i])/float(length)
    return rot_axis

#checks if the rows of a matrix are orthogonal
def isOrthogonal(v):
    orthogonality = 0
    for i in range(len(v[0].A1) - 1):
        for j in range(len(v[0].A1)):
            for k in range(len(v[0].A1)):
                if (i != j):
                    orthogonality += v[i].A1[k] * v[j].A1[k]
    return orthogonality

#Gets the angle between 2 1D vectors
def getAngle(v1,v2):
    return np.arccos(np.dot(v1,v2) /
                    (np.sqrt(np.dot(v1,v1)) * np.sqrt(np.dot(v2,v2))))

#rotates a matrix by any desired axis by any desired angle
def rotate(v, rot_axis, angle):
    rot_axis = makeUnitVector(np.array(rot_axis,dtype=np.float))
    R = np.matrix([[np.cos(angle)+rot_axis[0]**2 * (1 - np.cos(angle)),
                    rot_axis[0]*rot_axis[1] * (1 - np.cos(angle)) - rot_axis[2]*np.sin(angle),
                    rot_axis[0]*rot_axis[2] * (1 - np.cos(angle)) + rot_axis[1]*np.sin(angle)],

                    [rot_axis[1]*rot_axis[0] * (1 - np.cos(angle)) + rot_axis[2]*np.sin(angle),
                    np.cos(angle) + rot_axis[1]**2 * (1 - np.cos(angle)),
                    rot_axis[1]*rot_axis[2] * (1 - np.cos(angle)) - rot_axis[0]*np.sin(angle)],

                    [rot_axis[2]*rot_axis[0] * (1 - np.cos(angle)) - rot_axis[1]*np.sin(angle),
                    rot_axis[2]*rot_axis[1] * (1 - np.cos(angle)) + rot_axis[0]*np.sin(angle),
                    np.cos(angle) + rot_axis[2]**2 * (1 - np.cos(angle))] ])
    #print np.linalg.det(R)
    #print isOrthogonal(R)
    f = v * R
    return f[0].A1[0],f[0].A1[1],f[0].A1[2]

#Rotates a point and should be faster than my rotate function
#All credit for this function goes to Paul Bourke (I found it online)
def PointRotate3D(p1, p2, p0, theta):

    # Translate so axis is at origin
    p = [p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]]
    # Initialize point q
    q = [0.0,0.0,0.0]
    N = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    Nm = np.sqrt(N[0]**2 + N[1]**2 + N[2]**2)

    # Rotation axis unit vector
    n = [N[0]/Nm, N[1]/Nm, N[2]/Nm]

    # Matrix common factors
    c = np.cos(theta)
    t = (1 - np.cos(theta))
    s = np.sin(theta)
    X = n[0]
    Y = n[1]
    Z = n[2]

    # Matrix 'M'
    d11 = t*X**2 + c
    d12 = t*X*Y - s*Z
    d13 = t*X*Z + s*Y
    d21 = t*X*Y + s*Z
    d22 = t*Y**2 + c
    d23 = t*Y*Z - s*X
    d31 = t*X*Z - s*Y
    d32 = t*Y*Z + s*X
    d33 = t*Z**2 + c

    #This is just for testing if it is a pure rotation matrix: det = 1 and orthogonal
    #R = np.matrix([[d11,d12,d13],[d21,d22,d23],[d31,d32,d33]])
    #if ( abs(np.linalg.det(R) - 1) > 1e-8 or abs(isOrthogonal(R)) > 1e-8):
    #    print np.linalg.det(R)
    #    print isOrthogonal(R)

    #            |p.x|
    # Matrix 'M'*|p.y|
    #            |p.z|
    q[0] = d11*p[0] + d12*p[1] + d13*p[2]
    q[1] = d21*p[0] + d22*p[1] + d23*p[2]
    q[2] = d31*p[0] + d32*p[1] + d33*p[2]

    # Translate axis and rotated point back to original location
    return q[0] + p1[0],q[1] + p1[1],q[2] + p1[2]

def data(a,GB, angle):

    #I turn my GB into two points to use PointRotate3D
    p1 = [0.0,0.0,0.0]
    p2 = GB


    xx, yy = np.meshgrid(range(-a ,a +1,1),range(-a ,a +1,1))
    zz = []
    p0 = [0.0,0.0,0.0]

    for i in range(len(yy[0])):
        zz.append(np.copy(yy))

    holder = np.array(np.copy(zz),dtype=np.float)
    xx2 = np.array(np.copy(zz),dtype=np.float)
    yy2 = np.array(np.copy(zz),dtype=np.float)
    zz2 = np.array(np.copy(zz),dtype=np.float)

    for i in range(len(yy[0])):
        for j in range(len(yy[0])):
            for k in range(len(yy[0])):
                zz[i][j][k] = yy[i][j]
                p0 = [xx[j][k], yy[j][k],zz[i][j][k]]
                xx2[i][j][k], yy2[i][j][k],zz2[i][j][k] = PointRotate3D(p1, p2, p0, angle)
                holder[i][j][k] = yy[i][j]

    return xx,yy,zz,xx2,yy2,zz2,holder

def plot3D(a,GB, angle,nPLane,w,tolerance):
    xx,yy,zz,xx2,yy2,zz2,h = data(a,GB, angle)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_xlim(-a , a )
    ax.set_ylim(-a , a )
    ax.set_zlim(-a , a )
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    l = len(xx2[0]) - 1
#    #print ((l+1)/2+1)
#    cluster = ['o','^','*','h','D']
#    c = ['blue','green','red', 'yellow', 'black']
#    for i in range(l + 1):
#        for j in range(l + 1):
#            for k in range(l + 1):
#                if ( j + k == 2*a):
#                    #print xx2[i][j][k],yy2[i][j][k],zz2[i][j][k],(i,j,k), 20 * (5 - i), cluster[j], c[k]
#                    print xx[j][k],yy[j][k],zz[i][j][k],(i,j,k), 20 * (5 - i), cluster[j], c[k]
#                    ax.scatter(xx[j][k],yy[j][k],zz[i][j][k], s= 20 * (5 - i),marker = cluster[j], color = c[k])
#                    #ax.scatter(xx2[i][j][k],yy2[i][j][k],zz2[i][j][k], s= 20 * (3 - i),marker = cluster[j], color = c[k])
#                #if ((i!=j and j!=k and i!=k) or i==j==k==((l+1)/2+1)):
#                #    ax.scatter(xx2[i][j][k],yy2[i][j][k],zz2[i][j][k]), s= 20 * (3 - i),marker = cluster[j], color = c[k])
#
#    #GB = [1,-1,0]
#    #ax.plot([-a*GB[0],a*GB[0]], [-a*GB[1],a*GB[1]], [-a*GB[2],a*GB[2]], 'black')
#    #GB = [1,1,-2]
#    #ax.plot([-a*GB[0],a*GB[0]], [-a*GB[1],a*GB[1]], [-a*GB[2],a*GB[2]], 'black')


    for i in range(l + 1):
        ax.scatter(xx,yy,zz[i], marker = '.',color='blue')
        ax.scatter(xx2[i],yy2[i],zz2[i], marker = '^',color='green')

    for i in range(l + 1):
        for j in range(l + 1):
            ax.scatter(xx2[j][i][l-i],yy2[j][i][l-i],zz2[j][i][l-i],s=5*w,marker='^',color='green')
            ax.scatter(xx[i][i],yy[l-i][l-i],zz[j][i][i],s=5*w,marker='s',color='darkblue')

    for i in range(l+1):
        for j in range(l+1):
            ax.plot(yy[j],xx[j],zz[i][j],color='darkblue',linewidth = w)
            ax.plot(xx[j],yy[j],zz[i][j],color='darkblue',linewidth = w)
            ax.plot(yy[j],zz[i][j],xx[j],color='darkblue',linewidth = w)

    for i in range(l+1):
        for j in range(l+1):
            for k in range(l+1):
                ax.plot([xx2[l-i][l-j][l-k],xx2[l-i][l-j-1][l-k]],[yy2[l-i][l-j][l-k],yy2[l-i][l-j-1][l-k]],[zz2[l-i][l-j][l-k]
                        ,zz2[l-i][l-j-1][l-k]],color='green',linewidth = w)
                ax.plot([xx2[l-i][l-j][l-k],xx2[l-i-1][l-j][l-k]],[yy2[l-i][l-j][l-k],yy2[l-i-1][l-j][l-k]],[zz2[l-i][l-j][l-k]
                        ,zz2[l-i-1][l-j][l-k]],color='green',linewidth = w)
            ax.plot(xx2[i][j],yy2[i][j],zz2[i][j],color='green',linewidth = w)

    plt.show()

def plot001(a,GB, angle,nPlane,w,tolerance):
    GB = [0,0,1]
    xx,yy,h,xx2,yy2,zz2,zz = data(a,GB, angle)

    plt.xlim(-a * 2,a * 2)
    plt.ylim(-a * 2,a * 2)
    #plt.gca().set_linewidth(.2)
    plt.gca().set_aspect('equal', adjustable='box')


    ##Creates the new top and bot grain
    newTop = []
    newBot = []

    l = len(yy[0])
    for i in range(l):
        for j in range(l):
            newTop.append([xx[i][j],yy[i][j]])
            newBot.append([xx2[nPlane][i][j],yy2[nPlane][i][j]])

    ##PLots the atoms
    ln = len(newTop)
    for i in range(ln):
        plt.scatter(newTop[i][0],newTop[i][1],s=20*w,marker='s',color='darkblue')
        plt.scatter(newBot[i][0],newBot[i][1],s=20*w,marker='^',color='green')

    ##Joins the atoms along the 100 direction
    for i in range(0,ln,l):
        h = i + l - 1
        plt.plot([newTop[i][0],newTop[h][0]],[newTop[i][1],newTop[h][1]], linewidth = w, color = 'darkblue')
        plt.plot([newBot[i][0],newBot[h][0]],[newBot[i][1],newBot[h][1]],linewidth = w, color = 'green')


    ##Joins the atoms along the 010 direction
    for i in range(l):
        h = i + ln - l
        plt.plot([newTop[i][0],newTop[h][0]],[newTop[i][1],newTop[h][1]], linewidth = w, color = 'darkblue')
        plt.plot([newBot[i][0],newBot[h][0]],[newBot[i][1],newBot[h][1]],linewidth = w, color = 'green')


    ##Finds csl
    csl = []
    for i in range(ln):
        for j in range(ln):
            #print abs(newTop[i][0] - newBot[i][0])
            if (abs(newTop[i][0] - newBot[j][0]) < tolerance and abs(newTop[i][1] - newBot[j][1]) < tolerance):
                csl.append([newTop[i][0],newTop[i][1]])

    print len(csl), '# of csl'

    ##Plots csl
    for i in range(len(csl)):
        plt.scatter(csl[i][0],csl[i][1],s=100 * w, facecolors='none', edgecolors='red')

    plt.show()

def plot110(a,GB, angle, nPlane,w,tolerance):
    GB = [1,1,0]
    xy = [0,0,1]
    xx,yy,h,xx2,yy2,zz2,zz = data(a,GB, angle)



    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #
    #ax.set_xlim(-a , a )
    #ax.set_ylim(-a , a )
    #ax.set_zlim(-a , a )
    #ax.set_xlabel('X axis')
    #ax.set_ylabel('Y axis')
    #ax.set_zlabel('Z axis')

    ll = len(yy[0])

    newTop = []
    newBot = []

    rotatedGB = np.cross(GB,xy)
    phi = getAngle(GB,xy)
    p1 = [0.0,0.0,0.0]

    for i in range(ll):
        for j in range(ll):
            for k in range(ll):
                if ( j + k == nPlane ):
                    p0 = [xx2[i][j][k],yy2[i][j][k],zz2[i][j][k]]
                    newBot.append(PointRotate3D(p1, rotatedGB, p0, phi))
                    p0 = [xx[j][k],yy[j][k],zz[i][j][k]]
                    newTop.append(PointRotate3D(p1, rotatedGB, p0, phi))
                    #ax.scatter(xx2[j][i][ll-i],yy2[j][i][ll-i],zz2[j][i][ll-i],s=5*w,marker='^',color='green')
                    #ax.scatter(xx[i][i],yy[ll-i][ll-i],zz[j][i][i],s=5*w,marker='s',color='darkblue')
                    #plt.plot(xx[i],yy[i],color='darkblue',linewidth = w)
                    #plt.plot(yy[i],xx[i],color='darkblue',linewidth = w)

    ##This part will highlight the csls
    ln = len(newTop)
    csl = []
    for i in range(ln):
        for j in range(ln):
            #print abs(newTop[i][0] - newBot[i][0])
            if (abs(newTop[i][0] - newBot[j][0]) < tolerance and abs(newTop[i][1] - newBot[j][1]) < tolerance):
                csl.append([newTop[i][0],newTop[i][1]])

    print len(csl), '# of csl'

    plt.xlim(-a * 2,a * 2)
    plt.ylim(-a * 2,a * 2)
    #plt.gca().set_linewidth(.2)
    plt.gca().set_aspect('equal', adjustable='box')

    for i in range(len(csl)):
        plt.scatter(csl[i][0],csl[i][1],s=100 * w, facecolors='none', edgecolors='red')

    for i in range(ln):
        plt.scatter(newTop[i][0],newTop[i][1],s=20*w,marker='s',color='darkblue')
        plt.scatter(newBot[i][0],newBot[i][1],s=20*w,marker='^',color='green')


    #This loop will join the top lattice in the 1 -1 0 direction
    holder = ll - abs(ll - 1 - nPlane)
    for i in range(ll , 0 , -1):
        plt.plot([newTop[i * holder-1][0],newTop[(i - 1) * holder][0]],[newTop[i * holder-1][1],newTop[(i - 1) * holder][1]], linewidth = w, color = 'darkblue')
        plt.plot([newBot[i * holder-1][0],newBot[(i - 1) * holder][0]],[newBot[i * holder-1][1],newBot[(i - 1) * holder][1]],linewidth = w, color = 'green')

    #This loop will join the top lattice in the 0 0 1 direction
    for i in range(holder):
        plt.plot([newTop[i][0],newTop[i + (holder*(ll-1))][0]],[newTop[i][1],newTop[i + (holder*(ll-1))][1]], linewidth = w, color = 'darkblue')
        plt.plot([newBot[i][0],newBot[i + (holder*(ll-1))][0]],[newBot[i][1],newBot[i + (holder*(ll-1))][1]],linewidth = w, color = 'green')

    plt.show()

def plot111(a,GB, angle,nPlane,w,tolerance):
    GB = [1,1,1]
    xy = [0,0,1]

    #xx,yy,h,xx2,yy2,zz2,zz = data(a,GB, angle)

    plt.xlim(-a * 2,a * 2)
    plt.ylim(-a * 2,a * 2)
    #plt.gca().set_linewidth(.2)
    plt.gca().set_aspect('equal', adjustable='box')

    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #
    #ax.set_xlim(-a , a )
    #ax.set_ylim(-a , a )
    #ax.set_zlim(-a , a )
    #ax.set_xlabel('X axis')
    #ax.set_ylabel('Y axis')
    #ax.set_zlabel('Z axis')

    ll = len(yy[0])

    newTop = []
    newBot = []

    rotatedGB = np.cross(GB, xy)
    phi = getAngle(GB, xy)
    p1 = [0.0,0.0,0.0]

    for i in range(ll):
        for j in range(ll):
            for k in range(ll):
                if (i+j+k == nPlane):
                    p0 = [xx2[i][j][k],yy2[i][j][k],zz2[i][j][k]]
                    newBot.append(PointRotate3D(p1, rotatedGB, p0, phi))
                    p0 = [xx[j][k],yy[j][k],zz[i][j][k]]
                    newTop.append(PointRotate3D(p1, rotatedGB, p0, phi))
                    #ax.scatter(xx2[i][j][k],yy2[i][j][k],zz2[i][j][k],s=5*w,marker='^',color='green')
                    #ax.scatter(xx[j][k],yy[j][k],zz[i][j][k],s=5*w,marker='s',color='darkblue')
                    #plt.plot(xx[i],yy[i],color='darkblue',linewidth = w)
                    #plt.plot(yy[i],xx[i],color='darkblue',linewidth = w)

    ln = len(newTop)


    ##This part will highlight the csls
    csl = []
    for i in range(ln):
        for j in range(ln):
            if (abs(newTop[i][0] - newBot[j][0]) < tolerance and abs(newTop[i][1] - newBot[j][1]) < tolerance):
                csl.append([newTop[i][0],newTop[i][1]])

    for i in range(len(csl)):
        plt.scatter(csl[i][0],csl[i][1],s=100 * w, facecolors='none', edgecolors='red')



    ##This for loop will do a 3D image of the rotated 111 plane
    #for i in range(ln):
    #    ax.scatter(newTop[i][0],newTop[i][1],newTop[i][2],s=20*w,marker='s',color='darkblue')
    #    ax.scatter(newBot[i][0],newBot[i][1],newTop[i][2],s=20*w,marker='^',color='green')

    ##This for loop will do a 2D image of the rotated 111 plane
    for i in range(ln):
        plt.scatter(newTop[i][0],newTop[i][1],s=20*w,marker='s',color='darkblue')
        plt.scatter(newBot[i][0],newBot[i][1],s=20*w,marker='^',color='green')
        #plt.show()
        #plt.pause(.01)


    ##This will connect all the atoms along the 11-2 (y-axis)
    ##After much annoyance, I decided to collect the top dots of the connecting lines in endsTop.
    ##I did the same for endsBot

    endsTop = []
    endsBot = []

    if (nPlane > 2*a and nPlane < 4*a):
        step = ll - nPlane + 2*a
        h = 0
        stop = False
        for i in range(1,step,1):
            h += 1
            endsTop.append(h)

        midLength = (nPlane - 2*a)*2

        for i in range(midLength):
            if (i%2 == 0):
                h += step
                endsTop.append(h)
                if (step < ll and stop == False):
                    step += 1
                else:
                    stop = True
                    step -= 1
            else:
                h += 1
                endsTop.append(h)

        maxx = ll - nPlane + 2*a - 2
        for i in range(maxx):
            if (step < ll and stop == False):
                step += 1
            else:
                stop = True
                step -= 1
            h += step
            endsTop.append(h)


        step = ll - nPlane + 2*a
        h = 0
        stop = False
        maxx = abs(step - ll - 1)
        for i in range(1,maxx,1):
            h += step
            endsBot.append(h)
            if (step < ll and stop == False):
                step += 1
            else:
                stop = True
                step -= 1

        midLength = (4*a - nPlane)*2

        for i in range(midLength):
            if (i%2 == 1):
                if (step < ll and stop == False):
                    step += 1
                else:
                    stop = True
                    step -= 1
                h += step
                endsBot.append(h)

            else:
                h += 1
                endsBot.append(h)

        for i in range(maxx - 2):
            h += 1
            endsBot.append(h)


    elif (nPlane <= 2*a and nPlane >= 2):
        step = nPlane
        h = 0
        maxx = nPlane + 1
        for i in range(2,maxx,1):
            h = i
            endsTop.append(h)
        for i in range(0,maxx - 3,1):
            h += step
            endsTop.append(h)
            step -= 1

        step = nPlane + 1
        h = 0
        maxx = 2 * (nPlane - 2) + 1
        for i in range(maxx):
            if (i%2 ==0):
                h += step
                endsBot.append(h)
                step -= 1

            else:
                endsBot.append(h + 1)

    if (nPlane <= 6*a -2 and nPlane >= 4*a):
        step = 2
        h = 1
        maxx = 2 * (6*a - 2 - nPlane) + 1
        for i in range(maxx):
            if (i%2 ==1):
                h += step
                endsTop.append(h)
                step += 1

            else:
                h += 1
                endsTop.append(h)

        step = 2
        h = 1
        maxx = 6*a - 1  - nPlane
        for i in range(maxx):
            h += step
            endsBot.append(h)
            step += 1

        for i in range(maxx - 1):
            h += 1
            endsBot.append(h)






    if (nPlane >= 2 and nPlane <= 6*a -2):
        for i in range(len(endsTop)):
            plt.plot([newTop[endsTop[i]][0],newTop[endsBot[i]][0]],[newTop[endsTop[i]][1],newTop[endsBot[i]][1]], linewidth = w, color = 'darkblue')
            plt.plot([newBot[endsTop[i]][0],newBot[endsBot[i]][0]],[newBot[endsTop[i]][1],newBot[endsBot[i]][1]],linewidth = w, color = 'green')




    ##This will connect all the atoms along the 1-10 (x-axis)
    if (nPlane <= 2*a):
        step = 1
        holder = ln - 1
        h = 0
        for i in range(nPlane):
            plt.plot([newTop[h][0],newTop[holder][0]],[newTop[h][1],newTop[holder][1]], linewidth = w, color = 'darkblue')
            plt.plot([newBot[h][0],newBot[holder][0]],[newBot[h][1],newBot[holder][1]],linewidth = w, color = 'green')
            h += 1
            holder -= step
            step += 1
    elif (nPlane >= 4*a and nPlane <= 6*a):
        maxx = 6*a - nPlane
        step = 1
        holder = ln - 1
        h = 0
        for i in range(maxx):
            plt.plot([newTop[h][0],newTop[holder][0]],[newTop[h][1],newTop[holder][1]], linewidth = w, color = 'darkblue')
            plt.plot([newBot[h][0],newBot[holder][0]],[newBot[h][1],newBot[holder][1]],linewidth = w, color = 'green')
            h += step
            holder -= 1
            step += 1
    elif (nPlane > 2*a and nPlane < 4*a):
        maxx = 2*a
        step = abs(4*a +1 - nPlane)
        h = 0
        holder = ln - 1
        stop = False
        stop2 = False
        for i in range(maxx):
            plt.plot([newTop[h][0],newTop[holder][0]],[newTop[h][1],newTop[holder][1]], linewidth = w, color = 'darkblue')
            plt.plot([newBot[h][0],newBot[holder][0]],[newBot[h][1],newBot[holder][1]],linewidth = w, color = 'green')
            if (step < ll-1 and stop == False):
                h += step
                holder -= 1
                step += 1
            elif (step == ll-1 and stop == False):
                h += step
                holder -= 1
                stop = True
                plt.plot([newTop[h][0],newTop[holder][0]],[newTop[h][1],newTop[holder][1]], linewidth = w, color = 'darkblue')
                plt.plot([newBot[h][0],newBot[holder][0]],[newBot[h][1],newBot[holder][1]],linewidth = w, color = 'green')
                h = 1
                holder -= 1
                step = nPlane - 2 * a + 2

            else:
                h += 1
                holder -= step
                if (step < ll-1 and stop2 == False):
                    step += 1
                else:
                    step -= 1
                    stop2 = True





    plt.show()

def find_csl001(a,GB, angle,tolerance):
    GB = [0,0,1]
    xx,yy,h,xx2,yy2,zz2,zz = data(a,GB, angle)

    ll = len(yy[0])

    newTop = []
    newBot = []

    for i in range(ll):
        for j in range(ll):
            newTop.append([xx[i][j],yy[i][j]])
            newBot.append([xx2[nPlane][i][j],yy2[nPlane][i][j]])


    ln = len(newTop)

    ##This part will find the csls

    csl = []
    for i in range(ln):
        for j in range(ln):
            if (abs(newTop[i][0] - newBot[j][0]) < tolerance and abs(newTop[i][1] - newBot[j][1]) < tolerance):
                csl.append([newTop[i][0],newTop[i][1]])
    return len(csl)




def find_csl110(a,GB, angle,tolerance):
    GB = [1,1,0]
    xy = [0,0,1]
    xx,yy,h,xx2,yy2,zz2,zz = data(a,GB, angle)

    ll = len(yy[0])

    newTop = []
    newBot = []

    rotatedGB = np.cross(GB,xy)
    phi = getAngle(GB,xy)
    p1 = [0.0,0.0,0.0]

    for i in range(ll):
        for j in range(ll):
            for k in range(ll):
                if (j+k == 2*a):
                        p0 = [xx2[i][j][k],yy2[i][j][k],zz2[i][j][k]]
                        newBot.append(PointRotate3D(p1, rotatedGB, p0, phi))
                        p0 = [xx[j][k],yy[j][k],zz[i][j][k]]
                        newTop.append(PointRotate3D(p1, rotatedGB, p0, phi))

    ln = len(newTop)


    ##This part will find the csls

    csl = []
    for i in range(ln):
        for j in range(ln):
            if (abs(newTop[i][0] - newBot[j][0]) < tolerance and abs(newTop[i][1] - newBot[j][1]) < tolerance):
                csl.append([newTop[i][0],newTop[i][1]])
    return len(csl)


def find_csl111(a,GB, angle,tolerance):
    GB = [1,1,1]
    xy = [0,0,1]
    xx,yy,h,xx2,yy2,zz2,zz = data(a,GB, angle)

    ll = len(yy[0])

    newTop = []
    newBot = []

    rotatedGB = np.cross(GB, xy)
    phi = getAngle(GB, xy)
    p1 = [0.0,0.0,0.0]

    for i in range(ll):
        for j in range(ll):
            for k in range(ll):
                if (i+j+k == 3*a):
                    p0 = [xx2[i][j][k],yy2[i][j][k],zz2[i][j][k]]
                    newBot.append(PointRotate3D(p1, rotatedGB, p0, phi))
                    p0 = [xx[j][k],yy[j][k],zz[i][j][k]]
                    newTop.append(PointRotate3D(p1, rotatedGB, p0, phi))

    ln = len(newTop)


    ##This part will find the csls

    csl = []
    for i in range(ln):
        for j in range(ln):
            if (abs(newTop[i][0] - newBot[j][0]) < tolerance and abs(newTop[i][1] - newBot[j][1]) < tolerance):
                csl.append([newTop[i][0],newTop[i][1]])
    return len(csl)

def search():
    a = 15
    GB = np.array([1,1,1],dtype=np.float)
    angle = 55.0 * np.pi / 180.0
    step = .1
    counts = 4
    tolerance = 5e-3
    n = 3
    c = 0
    maximum = 50
    while(counts or c > maximum):
        angle += step * np.pi / 180.0
        csl = find_csl111(a,GB, angle,tolerance)
        if (csl > n):
            print angle * 180.0 / np.pi, csl
            counts -= 1
            angle -= step * np.pi / 180.0
            step *= .5
            tolerance *= .9
        c += 1
    tolerance = 1e-3
    return angle,tolerance


a = 25
GB = np.array([1,1,1],dtype=np.float)

#nPlane has to be ,between 0 and 2*a for 001, between 0 and 4*a for 110 and between 0 and 6*a for 111
nPlane = 3*a

angle = 60 * np.pi / 180.0
tolerance = 0.01

xx,yy,h,xx2,yy2,zz2,zz = data(a,GB, angle)

#for i in range(3*a,4*a + 1,1):
#
#    plt.figure(i)
#    nPlane = i
#    plot111(a,GB, angle,nPlane,.5,tolerance)
#    plt.scatter(0.0,0.0,s=200, marker='*', color='k')

#for i in range(4*a +1):
#
#    plt.figure(i)
#    nPlane = i
#    plot110(a,GB, angle,nPlane,.5,tolerance)

#for i in range(2*a +1):
#
#    plt.figure(i)
#    nPlane = i
#    plot001(a,GB, angle,nPlane,.5,tolerance)


#
#
#angle,tolerance = search()
#print angle  * 180/np.pi,tolerance

#angle,tolerance = 27.7957725 * np.pi / 180.0, 0.001
#plot111(a,GB, angle,2,1,.5,tolerance)

#y = float(1)
#x = float(2)
#angle = 2*np.arctan(y/x)

#angle = 50 * np.pi / 180.0
#tolerance = 0.01

# plot001(a,GB, angle,nPlane,.5,tolerance)
# plot110(a,GB, angle,nPlane,.5,tolerance)
plot111(a,GB, angle,nPlane,.5,tolerance)
#plot3D(a,GB, angle,nPlane,.5,tolerance)
