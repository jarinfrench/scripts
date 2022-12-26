#! /usr/bin/env python3

from dash import Dash, html, dcc, Input, Output
import plotly.graph_objects as go
import natsort
import pandas as pd
import numpy as np
from glob import glob
from copy import deepcopy
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

cmap = ['#66c2a5', '#fc8d62', '#8da0cb']

files = natsort.natsorted(glob('*.dump'))
files = files[0::5] # truncates the total amount of data to speed up computations

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
    .assign(radius = lambda x: np.sqrt(np.sum((x.iloc[:,2:4] - np.array(center[0:2]))**2, axis = 1)))

area_data = pd.read_csv('area_data.txt', sep = ' ',
    names = ['time', 'area1', 'area2'],
    dtype = {'time': 'int', 'area1': 'float', 'area2': 'float'},
    skiprows = 1).assign(timestep = lambda x: x.time / 0.002, radius = lambda x: np.sqrt(x.area2 / np.pi))

app = Dash(__name__)

# Store the snapshot figures so we aren't recalculating them all the time
snapshot_figs_all = {}
snapshot_figs_no_U = {}
for i in timesteps:
    snapshot_figs_all[i] = go.Figure()
    snapshot_figs_no_U[i] = go.Figure()
    snapshot_figs_all[i].update_layout(autosize = False, width = 600, height = 600, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
    snapshot_figs_all[i].update_xaxes(showgrid = False, showline = False, showticklabels = False, zeroline = False)
    snapshot_figs_all[i].update_yaxes(showgrid = False, showline = False, showticklabels = False, zeroline = False)
    snapshot_figs_no_U[i].update_layout(autosize = False, width = 600, height = 600, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
    snapshot_figs_no_U[i].update_xaxes(showgrid = False, showline = False, showticklabels = False, zeroline = False)
    snapshot_figs_no_U[i].update_yaxes(showgrid = False, showline = False, showticklabels = False, zeroline = False)
    selected_data_no_U = raw_data[(raw_data.timestep == i) & (raw_data.type != 1)]
    selected_data_all = raw_data[(raw_data.timestep == i) & (((raw_data.type == 1) & ((raw_data.id % 3) == 0)) | (raw_data.type != 1))]
    snapshot_figs_all[i].add_trace(go.Scatter(x = selected_data_all['x'],
                             y = selected_data_all['y'],
                             mode = 'markers',
                             marker_color = [cmap[i-1] for i in selected_data_all.type], marker_size = [3 if i == 1 else 6 for i in selected_data_all.type]))
    snapshot_figs_no_U[i].add_trace(go.Scatter(x = selected_data_no_U['x'],
                             y = selected_data_no_U['y'],
                             mode = 'markers',
                             marker_color = [cmap[i-1] for i in selected_data_no_U.type], marker_size = [3 if i == 1 else 6 for i in selected_data_no_U.type]))

default_fig = go.Figure()
default_fig.update_layout(autosize = False, width = 1000, height = 600, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
default_fig.update_xaxes(title = 'Time (ps)', range = [0, 8000])
default_fig.update_yaxes(title = 'Impurity concentration (percent)')
default_fig2 = go.Figure()
default_fig2.update_layout(autosize = False, width = 600, height = 600, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
default_fig2.update_xaxes(showgrid = False, showline = False, showticklabels = False, zeroline = False)
default_fig2.update_yaxes(showgrid = False, showline = False, showticklabels = False, zeroline = False)

app.layout = html.Div([
    html.H1('Impurity Concentration vs Time in the UMoXe System'),
    html.H2("Select the data to plot"),
    html.Div([
        html.Div([html.Label('Ring width:'), dcc.Input(
            id = 'ring-width-input',
            type = 'number',
            placeholder = 'Enter a number (minimum = 0.1)',
            debounce = True,
            min = 0.1,
            style = {'width': '260px'}
        )], style = {'width': '300px', 'float': 'left'}),
        html.Div([html.Label('Ring number:'), dcc.Input(
            id = 'ring-number-input',
            type = 'number',
            min = 1,
            step = 1
            ), html.Label('Max value: ', id = 'max-ring-num-label')], style = {'width': '240px', 'float': 'left', 'padding-left': '10px'}),
        html.Div([html.Label('Impurity type:'), dcc.RadioItems(
            [{"label": i, "value": f"N_imp{i}"} for i in range(1, n_imps + 1)],
            id = 'impurity-radio',
            value = 'N_imp1', # 1st impurity by default
            inline = True
        )], style = {'float': 'left'}, hidden = True if n_imps == 1 else False)
    ]),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div([
        html.Div([dcc.Loading(id = 'concentration-graph-loading', children = [html.Div(dcc.Graph(figure = default_fig, id = 'concentration-graph'))], type = 'circle'),
                  ],
                            style = {'float': 'left'}),
        html.Div([dcc.Graph(figure = default_fig2, id = 'snapshot-graph'),
                  dcc.Checklist([{'label': 'Hide U atoms', 'value': True}],
                                 id = 'hide-U-checklist')],
                style = {'float': 'left', 'padding-left': '10px', 'padding-top': '9px'}),

    ]),
    dcc.Store(id = 'concentration-data'),
    dcc.Store(id = 'concentration-figures'),
    dcc.Store(id = 'area-data')
])

@app.callback(
    Output('area-data', 'data'),
    Input('ring-width-input', 'value')
)
def update_area_data_rings(ring_width):
    if not ring_width:
        return None
    d = deepcopy(area_data)
    n_rings = int(np.ceil(max([x_bounds[1] - center[0], y_bounds[1] - center[1]]) / ring_width))
    d = d.assign(ring_num = lambda x: np.digitize(x.radius / ring_width, np.array(range(n_rings))))

    return d.reset_index().to_json(orient = 'split')

@app.callback(
    Output('concentration-data', 'data'),
    Output('ring-number-input', 'max'),
    Output('max-ring-num-label', 'children'),
    Output('concentration-figures', 'data'),
    Input('ring-width-input', 'value')
)
def update_data(ring_width):
    if not ring_width:
        return None, None, '', default_fig
    n_rings = int(np.ceil(max([x_bounds[1] - center[0], y_bounds[1] - center[1]]) / ring_width))
    # np.digitize is a fast way of binning the values x into the defined array
    tmp_data = raw_data.assign(ring_num = lambda x: np.digitize(x.radius / ring_width, np.array(range(n_rings))))
    multi_index = pd.MultiIndex.from_product([timesteps, range(1, n_rings + 1)], names = ["timestep", "ring_num"])
    concentration_data = pd.DataFrame(convertData(tmp_data, multi_index).values(), multi_index, columns = ["N_ring", "N_imp1", "N_imp2"])
    concentration_figs = {}
    for i in range(1,ring_width + 1):
        concentration_figs[i] = {}
        selected_keys = [key for key in concentration_data.index.values if key[1] == i]
        x = [j[0] * 0.002 for j in selected_keys]
        for imp_type in ['N_imp1', 'N_imp2']:
            y = [concentration_data.loc[i][imp_type] / concentration_data.loc[i]['N_ring'] * 100 for i in selected_keys]
            concentration_figs[i][imp_type] = go.Figure()
            concentration_figs[i][imp_type].update_layout(autosize = False, width = 1000, height = 600, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
            concentration_figs[i][imp_type].update_xaxes(title = 'Time (ps)', range = [0, 8000])
            concentration_figs[i][imp_type].update_yaxes(title = 'Impurity concentration (percent)', range = [0, yrange[ctotal]])
            concentration_figs[i][imp_type].add_trace(go.Scatter(x = x, y = y, mode = 'lines+markers'))

    return concentration_data.reset_index().to_json(orient = 'split'), n_rings, f'Max value: {n_rings}', concentration_figs

@app.callback(
    Output('concentration-graph', 'figure'),
    Output('concentration-graph', 'hoverData'),
    Input('concentration-data', 'data'),
    Input('area-data', 'data'),
    Input('ring-number-input', 'value'),
    Input('impurity-radio', 'value'),
    Input('concentration-figures', 'data')
)
def update_concentration_figure(c_data, a_data, ring_num, imp_type, c_figs):
    if not c_data or not a_data or not ring_num:
        return default_fig, None

    concentration_data = pd.read_json(c_data, orient = 'split')
    concentration_data = concentration_data.set_index(['timestep', 'ring_num'])
    area_info = pd.read_json(a_data, orient = 'split')
    selected_keys = [key for key in concentration_data.index.values if key[1] == ring_num]
    selected_times = area_info[area_info.ring_num == ring_num]['time'].values

    fig = go.Figure(c_figs[str(ring_num)][imp_type])
    y = fig.data[0].y

    if np.any(selected_times):
        median_time = np.median(selected_times)
        low_time = min(selected_times)
        high_time = max(selected_times)
        vline_min = min(y) - 10 if min(y) - 10 > 0 else 0
        vline_max = max(y) + 10
        fig.add_vrect(x0 = low_time, x1 = high_time, line_width = 0,
            fillcolor = "red", opacity = 0.2, annotation_text = 'GB within specified ring',
            annotation_position = 'top left')

    return fig, {'points': [{'curveNumber': 0, 'pointNumber': 0, 'pointIndex': 0, 'x': 0, 'y': 0}]}

@app.callback(
    Output('snapshot-graph', 'figure'),
    Input('concentration-graph', 'hoverData'),
    Input('ring-number-input', 'value'),
    Input('ring-width-input', 'value'),
    Input('hide-U-checklist', 'value')
)
def update_snapshot(hoverData, ring_num, ring_width, hide_U):
    if not hoverData or not ring_num or not ring_width:
        return default_fig2

    n_rings = int(np.ceil(max([x_bounds[1] - center[0], y_bounds[1] - center[1]]) / ring_width))
    timestep = int(hoverData['points'][0]['x'] / 0.002)
    if hide_U:
        fig = deepcopy(snapshot_figs_no_U[timestep])
    else:
        fig = deepcopy(snapshot_figs_all[timestep])

    low_rad = ring_width * (ring_num - 1)
    high_rad = ring_width * ring_num
    small_circle_x = center[0] + low_rad * cosangles
    small_circle_y = center[1] + low_rad * sinangles
    large_circle_x = center[0] + high_rad * cosangles
    large_circle_y = center[1] + high_rad * sinangles

    fig.add_trace(go.Scatter(x = small_circle_x, y = small_circle_y,
                             line = dict(color = 'Gray', dash = 'dash')))
    if not ring_num == n_rings:
        fig.add_trace(go.Scatter(x = large_circle_x, y = large_circle_y,
                                 line = dict(color = 'Gray', dash = 'dash')))
    fig.update_layout(showlegend = False, scene = dict(aspectmode = 'data'),
        title = {'text': f'Snapshot at t = {timestep * 0.002}ps',
                 'x': 0.5, 'y': 0.9,
                 'xanchor': 'center',
                 'yanchor': 'top'})
    return fig

if __name__ == '__main__':
    app.run_server(debug = True)
