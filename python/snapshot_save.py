import numpy as np
import matplotlib.pyplot as plt

# The following is from https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c
from matplotlib.text import Annotation
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.collections import PatchCollection

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

# This comes from https://stackoverflow.com/a/7968690
def forceAspect(ax,aspect=1.0):
    xl, xr = ax.get_xlim()
    yb, yt = ax.get_ylim()
    ax.set_aspect(abs((xr-xl)/(yb-yt))*aspect)

# f_list = ["0_interface.dat", "300000_interface.dat", "600000_interface.dat"] # 100
# f_list = ["0_interface.dat", "1000000_interface.dat", "2000000_interface.dat"] # 110
f_list = ["0_interface.dat", "150000_interface.dat", "300000_interface.dat"] # 111
times = [0.002*int(i.split("_")[0])/1000 for i in f_list]
# fl = [f"({i})" for i in "abc"] # 100
# fl = [f"({i})" for i in "def"] # 110
fl = [f"({i})" for i in "ghi"] # 111

data = [np.loadtxt(f,
    dtype={'names': ('ID', 'type', 'q', 'x', 'y', 'z', 'grain_num', 'op'),
           'formats': ('int', 'int', 'float', 'float', 'float', 'float', 'int', 'float')},
    skiprows = 1) for f in f_list]

fig = plt.figure(figsize=(10,8), constrained_layout = True)
ax1 = plt.subplot(131)
ax2 = plt.subplot(132)
ax3 = plt.subplot(133)

ax1.scatter(data[0]["x"], data[0]["y"], c = data[0]["grain_num"], cmap = 'Dark2', s = 0.1, zorder = 4.1)
ax2.scatter(data[1]["x"], data[1]["y"], c = data[1]["grain_num"], cmap = 'Dark2', s = 0.1, zorder = 4.1)
ax3.scatter(data[2]["x"], data[2]["y"], c = data[2]["grain_num"], cmap = 'Dark2', s = 0.1, zorder = 4.1)
min_x = min(data[0]["x"])
min_y = min(data[0]["y"])
x_tot_range = max(data[0]["x"]) - min_x
y_tot_range = max(data[0]["y"]) - min_y
frac = 0.05
axes_origin = ((2 * frac * x_tot_range) + min_x, (2 * frac * y_tot_range) + min_y)
z_dir_rad = 5
arrow_length = 2 * frac * x_tot_range
z_dir = []
z_dir.append(Circle(axes_origin, z_dir_rad, color = 'black', fill = False)) # outer edge of circle
z_dir.append(Circle(axes_origin, 2, color = 'black', fill = True))
coll = PatchCollection(z_dir, match_original = True, zorder = 4.9)
ax1.arrow(axes_origin[0] + z_dir_rad, axes_origin[1], arrow_length, 0, lw = 2, head_width = 6, color = "k", zorder = 4.9)
ax1.arrow(axes_origin[0], axes_origin[1] + z_dir_rad, 0, arrow_length, lw = 2, head_width = 6, color = "k", zorder = 4.9)
ax1.add_collection(coll)

x_dir = (1.5 * arrow_length + min_x, -0.025 * y_tot_range + min_y)
y_dir = (-2 * frac * x_tot_range + min_x, 1.6 * arrow_length + min_y)
z_dir = (-2 * frac * x_tot_range + min_x, -1.6 * frac * y_tot_range + min_y)
# Create the labels
# 100
# ax1.annotate(r"[100]", tuple(sum(x) for x in zip(axes_origin, x_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)
# ax1.annotate(r"[010]", tuple(sum(x) for x in zip(axes_origin, y_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)
# ax1.annotate(r"[001]", tuple(sum(x) for x in zip(axes_origin, z_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)

# 110
# ax1.annotate(r"[001]", tuple(sum(x) for x in zip(axes_origin, x_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)
# ax1.annotate(r"[1$\bf{\bar{1}}$0]", tuple(sum(x) for x in zip(axes_origin, y_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)
# ax1.annotate(r"[110]", tuple(sum(x) for x in zip(axes_origin, z_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)

# 111
ax1.annotate(r"[1$\bf{\bar{1}}$0]", tuple(sum(x) for x in zip(axes_origin, x_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)
ax1.annotate(r"[11$\bf{\bar{2}}$]", tuple(sum(x) for x in zip(axes_origin, y_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)
ax1.annotate(r"[111]", tuple(sum(x) for x in zip(axes_origin, z_dir)), weight = 'bold', fontsize = 14, zorder = 4.9)

# add the time (in ns) at the top of each subplot
ax1.annotate(f"{times[0]} ns", (max(data[0]["x"])/2 - 30, max(data[0]["y"] - 30.0)), weight = 'bold', fontsize = 14, zorder = 4.9)
ax2.annotate(f"{times[1]} ns", (max(data[1]["x"])/2 - 30, max(data[1]["y"] - 30.0)), weight = 'bold', fontsize = 14, zorder = 4.9)
ax3.annotate(f"{times[2]} ns", (max(data[2]["x"])/2 - 30, max(data[2]["y"] - 30.0)), weight = 'bold', fontsize = 14, zorder = 4.9)

# add the subfigure letters
ax1.annotate(fl[0], (0.02*max(data[0]["x"]), 0.9*max(data[0]["y"])), weight = 'bold', fontsize = 14, zorder = 4.9)
ax2.annotate(fl[1], (0.02*max(data[1]["x"]), 0.9*max(data[1]["y"])), weight = 'bold', fontsize = 14, zorder = 4.9)
ax3.annotate(fl[2], (0.02*max(data[2]["x"]), 0.9*max(data[2]["y"])), weight = 'bold', fontsize = 14, zorder = 4.9)

# view settings
ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

plt.show()
