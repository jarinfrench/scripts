import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from math import sin, cos, pi

def pointAtAngle(theta, r=1):
    return (r * cos(theta * pi / 180), r * sin(theta * pi / 180))

data = np.loadtxt("/media/jarinf/Research1/uo2/grain_growth/cylindrical/100/45degree/Basak/UO2_minimized_100_45.00degree_r100A.dat",
    dtype={'names': ('id', 'type', 'q', 'x', 'y', 'z'),
           'formats': ('int', 'int', 'float', 'float', 'float', 'float')},
    skiprows = 10)

selected = (data['x'] > max(data['x'])/2 + 1) & (data['y'] > max(data['y'])/2 + 1) & \
           (data['x'] < max(data['x']) - 1) & (data['y'] < max(data['y']) - 1) & \
           (data['type'] == 1)
xs = data['x'][selected]
ys = data['y'][selected]
# Note that getting this singular position out still returns an array, which FancyArrowPatch has an issue with
atom1 = (data['x'][data['id'] == 43326][0],data['y'][data['id'] == 43326][0])
atom2 = (data['x'][data['id'] == 46727][0],data['y'][data['id'] == 46727][0])
atom3 = (data['x'][data['id'] == 46140][0],data['y'][data['id'] == 46140][0])
atom4 = (data['x'][data['id'] == 48688][0],data['y'][data['id'] == 48688][0])

center = (max(data['x'])/2, max(data['y'])/2)
c2 = (0.75 * max(data['x']), 0.75 * max(data['y']))
pos_ins = [tuple(map(lambda i,j: i + j, center, pointAtAngle(a, 90))) for a in [80, 60, 30, 10]]
pos_outs = [tuple(map(lambda i,j: i + j, center, pointAtAngle(a, 110))) for a in [77, 58, 32, 13]]
arcs = [mpatches.FancyArrowPatch(a, b, connectionstyle = 'arc3,rad=-0.5', color = 'k', arrowstyle = '-|>', mutation_scale = 20) for a,b in zip(pos_ins, pos_outs)]

fig = plt.figure(figsize=(6,5))
gb = plt.Circle(center, 100, color = 'k', fill = False, linestyle = (0,(6,6)), lw = 2)
ax = plt.subplot(111)
ax.scatter(xs,ys, c = 'mediumaquamarine', s = 4, ) # the atoms
ax.add_patch(gb) # the GB
for arc in arcs:
    ax.add_patch(arc)
ax.arrow(*c2, *tuple(map(lambda i, j: (j - i)/3, c2, center)), lw = 3, head_width = 6, color = 'k', zorder = 4.9) # GB motion
ax.set_xlim([center[0], max(data['x'])])
ax.set_ylim([center[1], max(data['y'])])
ax.set_aspect('equal')
ax.set_axis_off()
plt.show()
