#! /usr/bin/env python
'''This script counts the percentages of impurities within 2nn
distance to generate an 'averaged' histogram'''

from collections import defaultdict
import pandas as pd
# import plotly.express as px
import plotly.graph_objects as go


files = [f'dir_{i}/nearest_neighbors.txt' for i in range(1,4)]

data = pd.concat((pd.read_csv(file, sep=' ',
                              names=('timestep', 'id', 'n_impurities', 'n_neighbors'),
                              dtype={'timestep': 'int',
                                     'id': 'int',
                                     'n_impurities': 'int',
                                     'n_neighbors': 'int'}).assign(run = i + 1) for i, file in enumerate(files)),
                 ignore_index = True).sort_values(['run', 'timestep', 'id'])

nn_counts = {}
for run in data.run.unique():
    nn_counts[run] = {}
    for t in data.timestep.unique():
        nn_counts[run][t] = {}
        for n in data.n_impurities.unique():
            nn_counts[run][t][n] = sum(
                data.query("run == @run and timestep == @t")['n_impurities'] == n)

avg_counts = {}
for run, d1 in nn_counts.items():
    d = defaultdict(lambda: 0)
    for timestep, d2 in d1.items():
        for n_neigh, count in d2.items():
            d[n_neigh] += count
    avg_counts[run] = {k: v / len(d1) for k, v in d.items()}

avg_data = pd.DataFrame.from_dict(avg_counts).sort_index()
x_labels = [f'{i} Impurity' if i == 1 else f'{i} Impurities' for i in avg_data.index.values]
fig = go.Figure(data = go.Bar(x = x_labels,
                              y = avg_data[avg_data.keys()[0]],
                              name = 'Run 1'))
for i in avg_data.keys()[1:]:
    fig.add_trace(go.Bar(x = x_labels, y = avg_data[i], name = f'Run {i}'))
fig.update_layout(yaxis = dict(title = 'Average number of 2nn Spheres across all timesteps',
                               type = 'log'))
# fig = go.Figure(data = go.Histogram(x = data[data['run'] == 1]['n_impurities'], name = 'Run 1'))
# fig.add_trace(go.Histogram(x = data[data['run'] == 2]['n_impurities'], name = 'Run 2'))
# fig.add_trace(go.Histogram(x = data[data['run'] == 3]['n_impurities'], name = 'Run 3'))
# fig.update_xaxes(tickvals = [0, 1/15, 2/15, 3/15, 4/15],
#                  ticktext = ['0 impurities', '1 impurity', '2 impurities',
#                          '3 impurities', '4 impurities'], range = (-1/50,1/3))
# fig.show()
fig.write_image('nearest_neighbor_histogram.png')
