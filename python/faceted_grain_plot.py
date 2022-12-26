#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from fractions import Fraction

# The following is from https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c
from matplotlib.text import Annotation
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.patches import FancyArrowPatch

class Annotation3D(Annotation):

    def __init__(self, text, xyz, *args, **kwargs):
        super().__init__(text, xy=(0, 0), *args, **kwargs)
        self._xyz = xyz

    def draw(self, renderer):
        x2, y2, z2 = proj_transform(*self._xyz, self.axes.M)
        self.xy = (x2, y2)
        super().draw(renderer)

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)

# Adds it to the Axes3D class


def _annotate3D(ax, text, xyz, *args, **kwargs):
    '''Add anotation `text` to an `Axes3d` instance.'''

    annotation = Annotation3D(text, xyz, *args, **kwargs)
    ax.add_artist(annotation)

setattr(Axes3D, 'annotate3D', _annotate3D)


def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


setattr(Axes3D, 'arrow3D', _arrow3D)

# files = ["lower_facet_atoms.dat", "lower_left_facet_atoms.dat", "upper_left_facet_atoms.dat", "upper_right_facet_atoms.dat"]
#
# for file in files:
#     xs = np.array([])
#     ys = np.array([])
#     zs = np.array([])
#     with open(file) as f:
#         for line in f:
#             if line.startswith("#"):  # ignore comment lines
#                 continue
#             l = [float(i) for i in line.split()]
#             xs = np.append(xs, l[3])
#             ys = np.append(ys, l[4])
#             zs = np.append(zs, l[5])
#     m, b = np.polyfit(xs, ys, 1)
    # print(f"{file}: y = {m}x + {b} or y = {Fraction(m).limit_denominator(10)}x + {b}")
    # print(f"   {min(xs)} <= x <= {max(xs)}, {min(ys)} <= y <= {max(ys)}")

def calc_final_point(p1, p2_1, p2_1_coord, slope):
    if p2_1_coord not in ["x", "y"]:
        print("Specify whether the x or y coordinate of the second point has been given")
        exit(1)
    if not len(p1) == 2:
        print("Points should be given in 2D format")
        exit(1)
    try:
        if not len(p2_1) == 1:
            print("Length of second argument should be 1")
    except TypeError:
        pass

    b = p1[1] - slope * p1[0]
    if p2_1_coord == "x":
        return (p2_1, slope * p2_1 + b)
    if p2_1_coord == "y":
        return ((p2_1 - b) / slope, p2_1)


atom_data = np.loadtxt("5000000_interface.dat",dtype={"names": ("id", "type", "charge", "x", "y", "z", "grain_num", "op"), "formats": ("int", "int", "float", "float", "float", "float", "int", "float")}, skiprows=1)
neg_vec_slope = np.sqrt(2) # slope for the UL and LR vectors
pos_vec_slope = -np.sqrt(2) # slope fo the UR and LL vectors
# upper left
UL_p1 = (0.15 * max(atom_data["x"]), 0.75 * max(atom_data["y"])) # 2/16 - 5/16 in x, 12/16 - (fit) in y
UL_p2 = calc_final_point(UL_p1, 0.25 * max(atom_data["x"]), "x", 1/neg_vec_slope)
upper_left_111_xs = np.linspace(UL_p1[0], UL_p2[0], 1000)
upper_left_111_ys = np.linspace(UL_p1[1], UL_p2[1], 1000)

# upper right
UR_p1 = (0.75 * max(atom_data["x"]), UL_p2[1]) # 11/16 - 14/16 in x, (fit) - 12/16 in y
UR_p2 = calc_final_point(UR_p1, 0.75 * max(atom_data["y"]), "y", 1/pos_vec_slope)
upper_right_111_xs = np.linspace(UR_p1[0], UR_p2[0], 1000)
upper_right_111_ys = np.linspace(UR_p1[1], UR_p2[1], 1000)

# lower right
LR_p1 = (UR_p1[0], 0.125 * max(atom_data["y"]))
LR_p2 = calc_final_point(LR_p1, UR_p2[0], "x", 1/neg_vec_slope)
lower_right_111_xs = np.linspace(LR_p1[0], LR_p2[0], 1000)
lower_right_111_ys = np.linspace(LR_p1[1], LR_p2[1], 1000)

# lower left
LL_p1 = (UL_p1[0], LR_p2[1])
LL_p2 = (UL_p2[0], LR_p1[1]) #calc_final_point(p1, 0.1875 * max(atom_data["y"]), "y", 1/pos_vec_slope)
lower_left_111_xs = np.linspace(LL_p1[0], LL_p2[0], 1000)
lower_left_111_ys = np.linspace(LL_p1[1], LL_p2[1], 1000)

center_x = (max(atom_data["x"]) - min(atom_data["x"]))/2
center_y = (max(atom_data["y"]) - min(atom_data["y"]))/2
center_z = (max(atom_data["z"]) - min(atom_data["z"]))/2
# upper left, upper right, lower right, lower left
normals_text_position=[(-110 + center_x,97 + center_y,0), (75 + center_x,97 + center_y,0), (75 + center_x,-121 + center_y,0), (-110 + center_x,-121 + center_y,0)]
normals_text = [r"[1$\bf{\bar{1}\bar{1}}$]", r"[1$\bf{\bar{1}}}$1]", r"[$\bf{\bar{1}}}$11]", r"[$\bf{\bar{1}}$1$\bf{\bar{1}}$]"]
normals=[(-48,48*np.sqrt(2),0), (48,48*np.sqrt(2),0), (40,-40*np.sqrt(2),0), (-40,-40*np.sqrt(2),0)]
x_means = [np.mean(upper_left_111_xs), np.mean(upper_right_111_xs), np.mean(lower_right_111_xs), np.mean(lower_left_111_xs)]
y_means = [np.mean(upper_left_111_ys), np.mean(upper_right_111_ys), np.mean(lower_right_111_ys), np.mean(lower_left_111_ys)]
facet_zs = np.ones(np.shape(upper_left_111_xs)) * (max(atom_data["z"]) + 10)
z_mean = np.mean(facet_zs)

# lower_facet_xs = np.linspace(113.568, 206.966, 1000)
# lower_facet_ys = np.linspace(49.7875, 72.5879, 1000)
# facet_zs = np.ones(np.shape(lower_facet_xs)) * (max(atom_data["z"]) + 10)
# lower_left_facet_xs = np.linspace(68.1769, 115.502, 1000)
# lower_left_facet_ys = np.linspace(93.9952, 48.2821, 1000)
# upper_left_facet_xs = np.linspace(65.7378, 150.362, 1000)
# upper_left_facet_ys = np.linspace(213.46, 258.052, 1000)
# upper_right_facet_xs = np.linspace(186.446, 237.537, 1000)
# upper_right_facet_ys = np.linspace(246.264, 202.203, 1000)
# center_x = (max(atom_data["x"]) - min(atom_data["x"]))/2
# center_y = (max(atom_data["y"]) - min(atom_data["y"]))/2
#
# # lower, lower left, upper left, upper right
# normals_text_position=[(-5 + center_x,-110 + center_y,0), (-80 + center_x,-99 + center_y,0), (-60 + center_x,90 + center_y,0), (40 + center_x,85 + center_y,0)]
# normals_text = [r"[$\bf{\bar{4}}41]$", r"[$\bf{\bar{1}}1\bf{\bar{1}}]$", r"[2$\bf{\bar{2}}\bf{\bar{1}}]$", r"[8$\bf{\bar{8}}$7]"]
# normals=[(17,-68,0), (-80,-80,0), (-38,76,0), (70,80,0)]
# x_means = [np.mean(lower_facet_xs), np.mean(lower_left_facet_xs), np.mean(upper_left_facet_xs), np.mean(upper_right_facet_xs)]
# y_means = [np.mean(lower_facet_ys), np.mean(lower_left_facet_ys), np.mean(upper_left_facet_ys), np.mean(upper_right_facet_ys)]
# z_mean = np.mean(facet_zs)

fig = plt.figure(figsize=(10,8))
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = plt.axes(projection = "3d", proj_type = 'ortho', computed_zorder = False)

# ax.scatter3D(atom_data["x"], atom_data["y"], atom_data["z"], c=atom_data["grain_num"], cmap = 'Dark2', s = 2, zorder = 4.1)
ax.scatter3D(atom_data["x"], atom_data["y"], atom_data["z"], c = atom_data["type"], cmap = 'Dark2', s = 2, zorder = 4.1) # used to double check that the vectors are parallel
ax.plot3D(upper_left_111_xs, upper_left_111_ys, facet_zs, "k-", lw = 3, zorder = 4.8)
ax.plot3D(upper_right_111_xs, upper_right_111_ys, facet_zs, "k-", lw = 3, zorder = 4.8)
ax.plot3D(lower_right_111_xs, lower_right_111_ys, facet_zs, "k-", lw = 3, zorder = 4.8)
ax.plot3D(lower_left_111_xs, lower_left_111_ys, facet_zs, "k-", lw = 3, zorder = 4.8)
for i,n in enumerate(normals):
    ax.arrow3D(x_means[i], y_means[i], z_mean, n[0], n[1], n[2], mutation_scale = 20, lw = 3, arrowstyle = "-|>", color = "k", zorder = 4.9)
    # ax.arrow3D(center_x, center_y, center_z, n[0], n[1], n[2], mutation_scale = 20, lw = 3, arrowstyle = "-|>", color = "k", zorder = 4.9)
    ax.annotate3D(normals_text[i], n, xytext = normals_text_position[i][0:2], textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)

# ax.plot3D(lower_facet_xs, lower_facet_ys, facet_zs, "k-", lw = 3, zorder=4.8)
# ax.plot3D(upper_left_facet_xs, upper_left_facet_ys, facet_zs, "k-", lw = 3, zorder=4.7)
# ax.plot3D(upper_right_facet_xs, upper_right_facet_ys, facet_zs, "k-", lw = 3, zorder=4.6)
# ax.plot3D(lower_left_facet_xs, lower_left_facet_ys, facet_zs, "k-", lw = 3, zorder=4.5)
# for i,n in enumerate(normals):
#     ax.arrow3D(x_means[i], y_means[i], z_mean, n[0], n[1], n[2], mutation_scale = 20, lw = 3, arrowstyle = "-|>", color = "k", zorder = 4.9)
#     ax.annotate3D(normals_text[i], n, xytext = normals_text_position[i][0:2], textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)
#
ax.arrow3D(20+max(atom_data["x"]), max(atom_data["y"])/2, max(atom_data["z"]), 40, 0, 0, mutation_scale = 20, lw = 3, arrowstyle = "->", color = "k", zorder = 4.9)
ax.arrow3D(22+max(atom_data["x"]), max(atom_data["y"])/2, max(atom_data["z"]), 0, 40, 0, mutation_scale = 20, lw = 3, arrowstyle = "->", color = "k", zorder = 4.9)
ax.annotate3D(r"[001]", (max(atom_data["x"]), max(atom_data["y"])/2, max(atom_data["z"])), xytext = (60,-5), textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)
ax.annotate3D(r"[1$\bf\bar{1}$0]", (max(atom_data["x"]), max(atom_data["y"])/2, max(atom_data["z"])), xytext = (5,42), textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)
# ax.annotate3D(r"(1)", (x_means[2], y_means[2], z_mean), xytext = (-35,0), textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)
# ax.annotate3D(r"(2)", (x_means[3], y_means[3], z_mean), xytext = (10,-5), textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)
# ax.annotate3D(r"(3)", (x_means[1], y_means[1], z_mean), xytext = (-35,-5), textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)
# ax.annotate3D(r"(4)", (x_means[0], y_means[0], z_mean), xytext = (10,-15), textcoords = 'offset points', zorder = 4.9, weight = 'bold', fontsize = 14)
ax.view_init(azim = -90, elev = 90) # top down view
ax.set_axis_off()
plt.draw()
plt.savefig("Basak_110_20degree_T2400_10ns.png", bbox_inches = fig.bbox_inches.from_bounds(1,1,8,6), facecolor="#AAAAAA")

plt.show()
