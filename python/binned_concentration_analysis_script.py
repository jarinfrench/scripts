#! /usr/bin/env python3

import plotly.graph_objects as go
import natsort
import pandas as pd
import numpy as np
from glob import glob
import os
from convertData import convertData

if not os.path.exists("c_vs_t_images"):
    os.mkdir('c_vs_t_images')

cdir = os.getcwd().split('/')
impurity_dir = cdir[6]
if impurity_dir == 'moly_effect':
    cMo = float(cdir[11].split('a')[0])
    cXe = 0
elif impurity_dir == 'xenon_effect':
    cMo = 0
    cXe = float(cdir[11].split('a')[0])
else:
    cMo = float(cdir[11].split('a')[0])
    cXe = float(cdir[12].split('a')[0])
ctotal = cMo + cXe
yrange = {0.25: 1, 0.5: 1.5, 0.75: 2, 1: 3, 1.25: 4, 1.5: 5, 1.75: 3,
          3: 5, 3.25: 5, 3.5: 5, 3.75: 5, 5: 7, 5.25: 7, 5.5: 7, 5.75: 7,
          10: 15, 14: 20, 17: 22, 21: 25}
def deg2rad(degree):
    return degree * np.pi / 180.0

deg2radv = np.vectorize(deg2rad)
angles = deg2radv(np.linspace(0,360,1000))
cosangles = np.cos(angles)
sinangles = np.sin(angles)

cmap = ['rgb(99,110,250)', 'rgb(239, 85, 59)']
cmap_fill = ['rgba(99,110,250, 0.2)', 'rgba(239, 85, 59, 0.2)']

files = natsort.natsorted(glob('*.dump'))
ring_width = 15.0
# files = files[0::5] # truncates the total amount of data to speed up computations

# we assume that these values are consistent between each data file, so only take
# them from the first one
x_bounds = []
y_bounds = []
z_bounds = []
with open(files[0]) as f:
    for line in f:
        if line.startswith("ITEM: BOX BOUNDS"):
            x_bounds = [float(i) for i in next(f).split()]
            y_bounds = [float(i) for i in next(f).split()]
            z_bounds = [float(i) for i in next(f).split()]

center = [(x_bounds[1] - x_bounds[0])/2, (y_bounds[1] - y_bounds[0])/2, (z_bounds[1] - z_bounds[0])/2]
n_rings = int(np.ceil(max([x_bounds[1] - center[0], y_bounds[1] - center[1]]) / ring_width))
ddtype = {"id": "int", "type": "int", "x": "float", "y": "float", "z": "float"}
n_imps = len(pd.read_csv(files[0], sep = ' ', names = tuple(ddtype.keys()), dtype = ddtype, skiprows = 9).type.unique()) - 1
timesteps = [int(file.split('.')[0]) for file in files]
raw_data = pd.concat(
    (pd.read_csv(file, sep = ' ',
                 names = tuple(ddtype.keys()),
                 dtype = ddtype,
                 skiprows = 9)\
    .assign(timestep = int(file.split('.')[0])) for file in files),
    ignore_index = True) \
    .assign(radius = lambda x: np.sqrt(np.sum((x.iloc[:,2:4] - np.array(center[0:2]))**2, axis = 1)),
        ring_num = lambda x: np.digitize(x.radius / ring_width, np.array(range(n_rings))))

area_data = pd.read_csv('area_data.txt', sep = ' ',
    names = ['time', 'area1', 'area2'],
    dtype = {'time': 'int', 'area1': 'float', 'area2': 'float'},
    skiprows = 1).assign(timestep = lambda x: x.time / 0.002,
        radius = lambda x: np.sqrt(x.area2 / np.pi),
        ring_num = lambda x: np.digitize(x.radius / ring_width, np.array(range(n_rings))))

with open('slope_calc.txt') as f:
    for line in f:
        pass
    gb_end = int(line.split()[-1].split(':')[1]) - 1

try:
    gb_end_time = area_data.iloc[gb_end]['time']
except IndexError:
    gb_end_time = area_data.iloc[gb_end - 1]['time']

multi_index = pd.MultiIndex.from_product([timesteps, range(1, n_rings + 1)], names = ["timestep", "ring_num"])
concentration_data = pd.DataFrame(convertData(raw_data, multi_index).values(), multi_index, columns = ["N_ring", "N_imp1", "N_imp2"])
all_fig = go.Figure()
all_fig.update_layout(autosize = False, width = 1000, height = 600, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4), title = {'text': f'Concentration over time with ringwidth = {ring_width}'})
zmax = 0

for nring in range(1, n_rings + 1):
    ymax = 0
    fig = go.Figure()
    fig.update_layout(autosize = False, width = 1000, height = 600, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4),
        title = {'text': f'Concentration over time for ring {nring}'})
    fig.update_xaxes(title = 'Time (ps)', range = [0, 8000])
    fig.update_yaxes(title = 'Impurity concentration (percent)', range = [0, yrange[ctotal]])
    selected_keys = [key for key in concentration_data.index.values if key[1] == nring]
    selected_times = area_data[area_data.ring_num == nring].time.values
    x = [i[0] * 0.002 for i in selected_keys]
    xerrs = x + x[::-1]
    for imp in range(n_imps):
        label = "Xe" if impurity_dir == 'xenon_effect' else "Mo" if impurity_dir == 'moly_effect' else "Mo" if imp == 0 else "Xe"
        y = [concentration_data.loc[i][f'N_imp{imp + 1}'] / concentration_data.loc[i]['N_ring'] * 100 for i in selected_keys]
        yerr = [concentration_data.loc[i][f'N_imp{imp+1}'] / concentration_data.loc[i]['N_ring'] * (1 / concentration_data.loc[i][f'N_imp{imp+1}'] + 1 / concentration_data.loc[i]['N_ring']) * 100 for i in selected_keys]
        yerrs = [i+j for i,j in zip(y,yerr)] + [i-j if i-j >= 0 else 0 for i,j in zip(y,yerr)][::-1]
        ymax = max(ymax,max(y))
        fig.add_trace(go.Scatter(x = x, y = y, mode = 'lines+markers', name = label, line = dict(color = cmap[imp]))) #, error_y = dict(type = 'data', array = yerrs, visible = True)))
        fig.add_trace(go.Scatter(x = xerrs, y = yerrs, fill = 'toself', fillcolor = cmap_fill[imp], line = dict(color = 'rgba(255,255,255,0)'), hoverinfo = 'skip', showlegend = False))
        fig.add_hline(y = cMo if label == 'Mo' else cXe, line = dict(color = 'black', dash = 'dash'))
        fig.add_vline(x = gb_end_time, line = dict(color = 'black', dash = 'dot'))

        all_fig.add_trace(go.Scatter3d(x = x, y = [nring for i in selected_keys], z = y, mode = 'lines+markers', name = label))
    if np.any(selected_times):
        median_time = np.median(selected_times)
        low_time = min(selected_times)
        high_time = max(selected_times)
        vline_min = min(y) - 10 if min(y) - 10 > 0 else 0
        vline_max = max(y) + 10
        fig.add_vrect(x0 = low_time, x1 = high_time, line_width = 0,
            fillcolor = "red", opacity = 0.2, annotation_text = 'GB within specified ring',
            annotation_position = 'top left')
    fig.write_image(f'c_vs_t_images/ring{nring}_ringwidth-{ring_width}.png')
