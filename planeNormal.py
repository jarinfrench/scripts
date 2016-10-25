from __future__ import division, print_function
import numpy as np
from matplotlib import cm # For changing the colormap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch

# This class was found at: http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


# For creating the x = 0 boundary plane
x100, z100 = np.meshgrid(range(-5,6), range(-5,6))
y100 = x100*0

fig1 = plt.figure(1)
ax1 = fig1.gca(projection='3d') # Makes the plot 3D
ax1.plot_surface(x100,y100,z100, cmap=cm.binary) # Plots a binary colored surface
ax1.set_xlabel("x", fontsize=28)
ax1.set_ylabel("y", fontsize=28)
ax1.set_zlabel("z", fontsize=28)
ax1.set_xlim3d(-5,5)
ax1.set_ylim3d(-5,5)
ax1.set_zlim3d(-5,5)
#Draws the vector.
norm100 = Arrow3D([0,0],[0,-2],[0,0], mutation_scale=20, lw=2, arrowstyle="-|>", color="k")
ax1.add_artist(norm100)
plt.title("<100> Tilt Boundary Plane\n with Normal", fontsize=30)
# Saves the figure as both a jpg and an eps for use in LaTex
plt.savefig("../myPapers/Senior Thesis/Images/Tilt100Plane.jpg")
plt.savefig("../myPapers/Senior Thesis/Images/Tilt100Plane.eps")

# Creates the y=x boundary plane
x110,z110 = np.meshgrid(range(-5,6), range(-5,6))
fig2 = plt.figure(2)
ax2 = fig2.gca(projection='3d')
ax2.plot_surface(x110,x110,z110, cmap=cm.binary)
ax2.set_xlabel("x", fontsize=28)
ax2.set_ylabel("y", fontsize=28)
ax2.set_zlabel("z", fontsize=28)
ax2.set_xlim3d(-5,5)
ax2.set_ylim3d(-5,5)
ax2.set_zlim3d(-5,5)
norm110 = Arrow3D([0,2],[0,-2],[0,0], mutation_scale=20, lw=2, arrowstyle="-|>", color="k")
ax2.add_artist(norm110)
plt.title("<110> Tilt Boundary Plane\n with Normal", fontsize=30)
plt.savefig("../myPapers/Senior Thesis/Images/Tilt110Plane.jpg")
plt.savefig("../myPapers/Senior Thesis/Images/Tilt110Plane.eps")

# Creates the y = 0 boundary plane
fig3 = plt.figure(3)
ax3 = fig3.gca(projection='3d')
ax3.plot_surface(y100,x100,z100, cmap=cm.binary)
ax3.set_xlabel("x", fontsize=28)
ax3.set_ylabel("y", fontsize=28)
ax3.set_zlabel("z", fontsize=28)
ax3.set_xlim3d(-5,5)
ax3.set_ylim3d(-5,5)
ax3.set_zlim3d(-5,5)
norm100Twist = Arrow3D([0,2],[0,0],[0,0], mutation_scale=20, lw=2, arrowstyle="-|>", color="k")
ax3.add_artist(norm100Twist)
plt.title("<100> Twist Boundary Plane\n with Normal", fontsize=30)
plt.savefig("../myPapers/Senior Thesis/Images/Twist100Plane.jpg")
plt.savefig("../myPapers/Senior Thesis/Images/Twist100Plane.eps")

# creates the y=-x boundary plane
fig4 = plt.figure(4)
ax4 = fig4.gca(projection='3d')
ax4.plot_surface(x110,-x110,z100, cmap=cm.binary)
ax4.set_xlabel("x", fontsize=28)
ax4.set_ylabel("y", fontsize=28)
ax4.set_zlabel("z", fontsize=28)
ax4.set_xlim3d(-5,5)
ax4.set_ylim3d(-5,5)
ax4.set_zlim3d(-5,5)
norm110Twist = Arrow3D([0,1],[0,1],[0,0], mutation_scale=20, lw=2, arrowstyle="-|>", color="k")
ax4.add_artist(norm110Twist)
plt.title("<110> Twist Boundary Plane\n with Normal", fontsize=30)
plt.savefig("../myPapers/Senior Thesis/Images/Twist110Plane.jpg")
plt.savefig("../myPapers/Senior Thesis/Images/Twist110Plane.eps")

# Creates the 111 boundary plane
# Code modified from http://stackoverflow.com/questions/3461869/plot-a-plane-based-on-a-normal-vector-and-a-point-in-matlab-or-matplotlib
point  = np.array([0, 0, 0])
normal = np.array([1, 1, 1])

# a Boundary Plane is a*x+b*y+c*z+d=0
# [a,b,c] is the normal. Thus, we have to calculate
# d and we're set
d = -point.dot(normal)

# create x,y
x111, y111 = np.meshgrid(range(-4,5), range(-4,5))

# calculate corresponding z
z111 = (-normal[0] * x111 - normal[1] * y111 - d) * 1. /normal[2]

fig5 = plt.figure(5)
ax5 = fig5.gca(projection='3d')
ax5.plot_surface(x111,y111,z111, cmap=cm.binary)
ax5.set_xlabel("x", fontsize=28)
ax5.set_ylabel("y", fontsize=28)
ax5.set_zlabel("z", fontsize=28)
ax5.set_xlim3d(-5,5)
ax5.set_ylim3d(-5,5)
ax5.set_zlim3d(-5,5)
norm111Twist = Arrow3D([0,1],[0,1],[0,1], mutation_scale=20, lw=2, arrowstyle="-|>", color="k")
ax5.add_artist(norm111Twist)
plt.title("<111> Twist Boundary Plane\n with Normal", fontsize=30)
plt.savefig("../myPapers/Senior Thesis/Images/Twist111Plane.jpg")
plt.savefig("../myPapers/Senior Thesis/Images/Twist111Plane.eps")

plt.show()
