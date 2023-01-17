#! /usr/bin/env python3

import numpy as np
import plotly.express as px

data = np.loadtxt('nearest_neighbors.txt', dtype = {'names': ('id', 'fraction'), 'formats': ('int', 'float')})
fig = px.scatter(data, x = 'id', y = 'fraction', marginal_y = 'histogram')
fig.write_image('nearest_neighbors.png')
