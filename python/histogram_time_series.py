#! /usr/bin/env python3.9

import argparse, itertools
import matplotlib.pyplot as plt
import numpy as np
import regex as re
import warnings
from math import pi
from matplotlib.widgets import Slider, Button

# Based roughly on https://stackoverflow.com/a/42764815 and https://matplotlib.org/stable/gallery/widgets/slider_snap_demo.html
# See also https://stackoverflow.com/questions/44412296/update-bar-chart-using-a-slider-in-matplotlib
# and https://stackoverflow.com/questions/42764278/how-to-update-a-histogram-when-a-slider-is-used
def update(val):
    ts = stime.val
    ax.clear()

    if -1 in [args.lx, args.ly, args.lz]:
        bars = ax.bar(histogram_all_data[ts]['Bin Number'], histogram_all_data[ts]['Count Normalized Concentration'])
        ax.set_title(f"Total Concentration (by count) = {total_concentration_by_count:.4f}")
        ax.set_ylim([0, 1.1 * max_y_by_count])
        ax.set_ylabel("Count Normalized Concentration")
    else:
        bars = ax.bar(histogram_all_data[ts]['Bin Number'], histogram_all_data[ts]['Volume Normalized Concentration'])
        ax.set_title(f"Total Concentration (by volume) = {total_concentration_by_volume:.4f}")
        ax.set_ylim([0, 1.1 * max_y_by_volume])
        ax.set_ylabel("Volume Normalized Concentration")
    ax.set_xlabel(f"Radius (x {binsize})")

    grain_radius = np.sqrt(area_data["cylinder"][np.where(area_data["time"] == ts * args.dt)][0] / pi)
    grain_radius_bin_num = -1

    for i,j in zip(pairwise(histogram_all_data[ts]['Bin Number']), pairwise(histogram_all_data[ts]['Bin Radius'])): # i will store the bin numbers, j will store the bin radius
        if grain_radius >= j[0] and grain_radius < j[1]:
            grain_radius_bin_num = int(i[0])
    if grain_radius_bin_num != -1:
        ax.get_children()[grain_radius_bin_num].set_color('r')

    legend_colors = {'Estimated Grain Radius': 'red'}
    legend_labels = list(legend_colors.keys())
    legend_handles = [plt.Rectangle((0,0), 1, 1, color = legend_colors[label]) for label in legend_labels]
    ax.legend(legend_handles, legend_labels)

    plt.draw()

def forward(event):
    pos = stime.val
    idx = allowed_timesteps.index(pos)
    if idx != len(allowed_timesteps) - 1:
        stime.set_val(allowed_timesteps[idx + 1])

def backward(event):
    pos = stime.val
    idx = allowed_timesteps.index(pos)
    if idx != 0:
        stime.set_val(allowed_timesteps[idx - 1])

def forward5(event):
    pos = stime.val
    idx = allowed_timesteps.index(pos)
    if not idx >= len(allowed_timesteps) - 5:
        stime.set_val(allowed_timesteps[idx + 5])
    else:
        stime.set_val(allowed_timesteps[-1])

def backward5(event):
    pos = stime.val
    idx = allowed_timesteps.index(pos)
    if not idx < 5:
        stime.set_val(allowed_timesteps[idx - 5])
    else:
        stime.set_val(allowed_timesteps[0])

def save_snapshot(event):
    t = stime.val
    fig.savefig(f"histogram_t{t}.png")

# def swap(event):


# from https://stackoverflow.com/a/5434936
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

# modified slightly from https://stackoverflow.com/a/61306932
def Cartesian2Polar(x, y, x0 = 0, y0 = 0):
    '''
    cartisian to polar coordinate system with origin shift to x0,y0
    '''
    x1 = x - x0
    y1 = y - y0
    #print('(x0,y0)sort',x0,y0)
    r = np.sqrt(x1**2 + y1**2)
    t = np.arctan2(y1,x1) * 180 / pi
    if y1 < 0:
        t = 360 + t
    return r,t

def AntiClockwiseSort(xy_list, x0 = None, y0 = None):
    '''
    Sort points Anti clockwise with x0 y0 as origin
    '''
    if x0 is None and y0 is None:
        (x0,y0) = np.mean(xy_list,axis = 0).tolist()
    elif x0 is None:
        (x0,_) = np.mean(xy_list,axis = 0).tolist()
    elif y0 is None:
        (_,y0) = np.mean(xy_list,axis = 0).tolist()
    #print('origin used:',[x0, y0])

    for i in range(len(xy_list)):
          xy_list[i].append(i)

    xy_list1 = sorted(xy_list, key=lambda a_entry: Cartesian2Polar(a_entry[0], a_entry[1], x0, y0)[1])

    sort_index = []
    for x in xy_list1:
          sort_index.append(x[2])
          del x[2]

    return xy_list1, sort_index

def getIntersectionPoints(r, lx, ly, center):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        points = [[0, ly * center[1] + np.sqrt(r**2 - (lx * center[0]) ** 2)], # left boundary
            [0, ly * center[1] - np.sqrt(r**2 - (lx * center[0]) ** 2)], # left boundary
            [lx, ly * center[1] + np.sqrt(r**2 - (lx * (1 - center[0])) ** 2)], # right boundary
            [lx, ly * center[1] - np.sqrt(r**2 - (lx * (1 - center[0])) ** 2)], # right boundary
            [lx * center[0] + np.sqrt(r**2 - (ly * center[1]) ** 2), 0], # bottom boundary
            [lx * center[0] - np.sqrt(r**2 - (ly * center[1]) ** 2), 0], # bottom boundary
            [lx * center[0] + np.sqrt(r**2 - (ly * (1 - center[1])) ** 2), ly], # top boundary
            [lx * center[0] - np.sqrt(r**2 - (ly * (1 - center[1])) ** 2), ly] # top boundary
        ]
    points = [[i,j] for i,j in points if not np.isnan(i) and not np.isnan(j)]
    points, _ = AntiClockwiseSort(points, lx * center[0], ly * center[1])
    return [[i,j] for i,j in points if i >= 0 and j >= 0]


def calcWedgeVols(r, lx, ly, lz, center):
    # first we need to determine the inersection points
    points = getIntersectionPoints(r, lx, ly, center)
    wedges_vol = []

    # for each pair of points, calculate the wedge volume
    if len(points) == 1: # in the case of only one point being valid, this means that the disk hits the box edge exactly once
        # the only time this will happen is when the disk fits within the box, and the box is tangent to the disk at this point
        # So, we add the full volume of the disk.
        wedges_vol.append(pi * r * r)
    else:
        points.append(points[0]) # we append the first point to the end to allow for pairwise iteration over all wedges
        # whichever side the pair of points are on (e.g. if both are on the left boundary)
        # the circle extends past the domain boundary on that side, so we would only need
        # the area of the triangle (volume of the prism) on that side.
        for i,j in pairwise(points):
            # if the pair of points are on the same side (e.g. i[0] == j[0], or i[1] == j[1])
            # then we need to add the volume of the triangle the wedge makes with the boundary
            if abs(i[0] - j[0]) < 1e-8: # left or right boundary
                chord = abs(i[1] - j[1])
                wedges_vol.append(lz * 0.5 * chord * np.sqrt(r * r - (0.5 * chord) ** 2))
            elif abs(i[1] - j[1]) < 1e-8: # top or bottom boundary
                chord = abs(i[0] - j[0])
                wedges_vol.append(lz * 0.5 * chord * np.sqrt(r * r - (0.5 * chord) ** 2))
            else: # if they are not on the same side, we add the full volume of the wedge
                v1 = [i[0] - lx * center[0], i[1] - ly * center[1]]
                v2 = [j[0] - lx * center[0], j[1] - ly * center[1]]
                v1_mag = np.sqrt(sum([i*i for i in v1]))
                v2_mag = np.sqrt(sum([i*i for i in v2]))
                theta = np.arccos((v1[0]*v2[0] + v1[1] * v2[1])/(v1_mag * v2_mag))
                wedges_vol.append(theta * 0.5 * r * r * lz)

    return wedges_vol

def calcDiskVol(r, lx, ly, lz, center):
    if r <= lx * center[0] and r <= lx * (1 - center[0]) and r <= ly * center[1] and r <= ly * (1 - center[1]):
        return pi * r * r * lz
    else:
        return sum(calcWedgeVols(r, lx, ly, lz, center))

def calcRingVol(r1, r2, lx, ly, lz, center = (0.5, 0.5)):
    assert len(center) == 2, "center should be a tuple of length 2"
    assert center[0] <=1 and center[0] >= 0 and center[1] <= 1 and center[1] >= 0, "center should be given in fractional coordinates (0 <= value <= 1)"
    if -1 in [lx, ly, lz]:
        return 1
    if r1 > r2: # make sure r1 is smaller than r2
        r1, r2 = r2, r1

    return calcDiskVol(r2, lx, ly, lz, center) - calcDiskVol(r1, lx, ly, lz, center)


parser = argparse.ArgumentParser(usage = '%(prog)s [-h] file', description = "Visualize the histogram time series of a set of data")
parser.add_argument('file', help = "The file containing the time series histogram data")
parser.add_argument('-x', '--lx', type = float, default = -1, help = "The length of the x dimension of the simulation")
parser.add_argument('-y', '--ly', type = float, default = -1, help = "The length of the y dimension of the simulation")
parser.add_argument('-z', '--lz', type = float, default = -1, help = "The length of the z dimension of the simulation")
parser.add_argument('-a', dest = 'a0', type = float, default = 1, help = "The lattice parameter of the material")
parser.add_argument('-t', '--dt', type = float, default = 0.002, help = "The timestep used for these dump files (default = 0.002)")
parser.add_argument('-c', '--structure', choices = ["bcc", "fcc", "sc"], default = "bcc", help = "The crystal structure of the system (default: bcc)")
parser.add_argument('--area-file', default = "area_data.txt", help = "The datafile containing the number of atoms in the grain over time. (default = area_data.txt)")
parser.add_argument('-s', '--save', nargs = '?', default = "histogram_time_series.mp4", help = "The file name of the video file")

args = parser.parse_args()

structure_factor = {"bcc": 2, "fcc": 4, "sc": 1}
histogram_all_data = dict()
histogram_single_data = dict()
first = True
first_timestep = 0
max_y_by_count = -1
max_y_by_volume = -1
binsize = 0
ring_vol = 0

area_data = np.loadtxt(args.area_file, dtype = {'names': ('time', 'matrix', 'cylinder'), 'formats': ('int', 'float', 'float')})

with open(args.file) as f:
    for line in f:
        if line.startswith('#'): # ignore comment lines
            l = [e for e in re.split('\(|\)|\=| ', line) if e]
            if "file" in line:
                filename = l[l.index("file") + 1]
                timestep = int(filename.split(".")[0])
                if first:
                    first_timestep = timestep
                    first = False
            if "binsize" in line:
                binsize = float(l[l.index("binsize") + 1])
            continue
        elif line == '\n': # the line is simply '\n',
            histogram_all_data[timestep] = histogram_single_data
            histogram_single_data = dict()
            continue
        else: # the actual data
            bin_num, type_count, all_count = [int(i) for i in line.split()]
            if bin_num == 0:
                ring_vol = calcRingVol(0, binsize, args.lx, args.ly, args.lz) * structure_factor[args.structure] / args.a0 ** 3
                normalized_concentration = type_count / all_count
                normalized_concentration_by_volume = type_count / ring_vol
                histogram_single_data['Bin Number'] = [bin_num]
                histogram_single_data['Bin Radius'] = [binsize]
                histogram_single_data['Type Count'] = [type_count]
                histogram_single_data['Total Count'] = [all_count]
                histogram_single_data['Count Normalized Concentration'] = [normalized_concentration]
                histogram_single_data['Volume Normalized Concentration'] = [normalized_concentration_by_volume]
            else:
                ring_vol = calcRingVol(bin_num * binsize, (bin_num + 1) * binsize, args.lx, args.ly, args.lz) * structure_factor[args.structure] / args.a0 ** 3
                normalized_concentration = type_count / all_count #* ring_vol / all_count
                normalized_concentration_by_volume = type_count / ring_vol
                histogram_single_data['Bin Number'].append(bin_num)
                histogram_single_data['Bin Radius'].append((bin_num + 1) * binsize)
                histogram_single_data['Type Count'].append(type_count)
                histogram_single_data['Total Count'].append(all_count)
                histogram_single_data['Count Normalized Concentration'].append(normalized_concentration)
                histogram_single_data['Volume Normalized Concentration'].append(normalized_concentration_by_volume)


# set up the plot stuff - based on https://matplotlib.org/stable/gallery/widgets/slider_snap_demo.html
fig, ax = plt.subplots()
plt.subplots_adjust(bottom = 0.25) # my guess (haven't checked this) is that this adds space to the bottom of the plot figure
total_type = sum(histogram_all_data[first_timestep]['Type Count'])
total_count = sum(histogram_all_data[first_timestep]['Total Count'])
total_concentration_by_count = total_type / total_count
total_concentration_by_volume = total_type / (args.lx * args.ly * args.lz / structure_factor[args.structure])
for key in histogram_all_data.keys():
    histogram_all_data[key]['Count Normalized Concentration'] = [i / total_concentration_by_count for i in histogram_all_data[key]['Count Normalized Concentration']]
    histogram_all_data[key]['Volume Normalized Concentration'] = [i / total_concentration_by_volume for i in histogram_all_data[key]['Volume Normalized Concentration']]
    if max(histogram_all_data[key]['Count Normalized Concentration']) > max_y_by_count:
        max_y_by_count = max(histogram_all_data[key]['Count Normalized Concentration'])
    if max(histogram_all_data[key]['Volume Normalized Concentration']) > max_y_by_volume:
        max_y_by_volume = max(histogram_all_data[key]['Volume Normalized Concentration'])
if -1 in [args.lx, args.ly, args.lz]:
    h = plt.bar(histogram_all_data[first_timestep]['Bin Number'], histogram_all_data[first_timestep]['Count Normalized Concentration'])
    ax.set_title(f"Total Concentration (by count) = {total_concentration_by_count:.4f}")
    ax.set_ylim([0, 1.1 * max_y_by_count])
    ax.set_ylabel("Count Normalized Concentration")
else:
    h = plt.bar(histogram_all_data[first_timestep]['Bin Number'], histogram_all_data[first_timestep]['Volume Normalized Concentration'])
    ax.set_title(f"Total Concentration (by volume) = {total_concentration_by_volume:.4f}")
    ax.set_ylim([0, 1.1 * max_y_by_volume])
    ax.set_ylabel("Volume Normalized Concentration")

ax.set_xlabel(f"Radius (x {binsize})")
grain_radius = np.sqrt(area_data["cylinder"][np.where(area_data["time"] == first_timestep * args.dt)][0] / pi)
grain_radius_bin_num = -1

for i,j in zip(pairwise(histogram_all_data[first_timestep]['Bin Number']), pairwise(histogram_all_data[first_timestep]['Bin Radius'])): # i will store the bin numbers, j will store the bin radius
    if grain_radius >= j[0] and grain_radius < j[1]:
        grain_radius_bin_num = int(i[0])
if grain_radius_bin_num != -1:
    ax.get_children()[grain_radius_bin_num].set_color('r')

legend_colors = {'Estimated Grain Radius': 'red'}
legend_labels = list(legend_colors.keys())
legend_handles = [plt.Rectangle((0,0), 1, 1, color = legend_colors[label]) for label in legend_labels]
ax.legend(legend_handles, legend_labels)

slider_bkd_color = "white"
ax_timestep = plt.axes([0.15,0.1,0.65,0.03], facecolor = slider_bkd_color)
ax_btnforward = plt.axes([0.35, 0.025, 0.05, 0.04])
ax_btnbackward = plt.axes([0.25, 0.025, 0.05, 0.04])
ax_btnforward5 = plt.axes([0.45, 0.025, 0.05, 0.04])
ax_btnbackward5 = plt.axes([0.15, 0.025, 0.05, 0.04])
ax_btnSave = plt.axes([0.60, 0.025, 0.1, 0.04])
# ax_swap = plt.axes([0.65, 0.025, 0.1, 0.04])
allowed_timesteps = [*histogram_all_data] # allowed values for snapping

# create the slider
stime = Slider(
    ax_timestep,
    "Timestep",
    min(allowed_timesteps),
    max(allowed_timesteps),
    valinit = first_timestep,
    valstep = allowed_timesteps,
    color = "black",
    initcolor = 'none'
)

bforward = Button(
    ax_btnforward,
    '+1',
    color = slider_bkd_color,
    hovercolor = 'lightgoldenrodyellow'
)

bbackward = Button(
    ax_btnbackward,
    '-1',
    color = slider_bkd_color,
    hovercolor = 'lightgoldenrodyellow'
)

bforward5 = Button(
    ax_btnforward5,
    '+5',
    color = slider_bkd_color,
    hovercolor = 'lightgoldenrodyellow'
)

bbackward5 = Button(
    ax_btnbackward5,
    '-5',
    color = slider_bkd_color,
    hovercolor = 'lightgoldenrodyellow'
)

bSave = Button(
    ax_btnSave,
    'save',
    color = slider_bkd_color,
    hovercolor = 'lightgoldenrodyellow'
)

# bswap = Button(
#     ax_swap,
#     'Swap',
#     color = slider_bkd_color,
#     hovercolor = 'lightgoldenrodyellow'
# )

stime.on_changed(update)
bforward.on_clicked(forward)
bbackward.on_clicked(backward)
bforward5.on_clicked(forward5)
bbackward5.on_clicked(backward5)
bSave.on_clicked(save_snapshot)
# bswap.on_clicked(swap)
plt.show()
