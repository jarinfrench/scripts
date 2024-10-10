from math import atan2

import numpy as np
import plotly.graph_objects as go

npoints = 200
centerpoint = (50, 50)
grain_radius_initial = 35
grain_radius_final = 25
radius_fudge_factor = 2
xdata_initial = [np.random.randint(0, 100) for i in range(npoints)]
ydata_initial = [np.random.randint(0, 100) for i in range(npoints)]
data_initial = np.stack((xdata_initial, ydata_initial), axis=1)
radii = [
    np.sqrt((i[0] - centerpoint[0]) ** 2 + (i[1] - centerpoint[1]) ** 2)
    for i in data_initial
]
angles = [
    atan2(i[1] - centerpoint[1], i[0] - centerpoint[0]) - atan2(0, centerpoint[0])
    for i in data_initial
]
adj_angles = [i if i >= 0 else i + 2 * np.pi for i in angles]
point_in_growth = [
    True
    if i >= grain_radius_final - radius_fudge_factor
    and i <= grain_radius_initial + radius_fudge_factor
    else False
    for i in radii
]
shifted_points_idx = [idx for idx, val in enumerate(point_in_growth) if val]

xdata_frame2 = [
    i + np.random.random() * np.random.choice([-1, 1])
    if not j
    else centerpoint[0]
    + grain_radius_final * np.cos(k)
    + np.random.random() * np.random.choice([-1, 1])
    for i, j, k in zip(xdata_initial, point_in_growth, adj_angles)
]
ydata_frame2 = [
    i + np.random.random() * np.random.choice([-1, 1])
    if not j
    else centerpoint[1]
    + grain_radius_final * np.sin(k)
    + np.random.random() * np.random.choice([-1, 1])
    for i, j, k in zip(ydata_initial, point_in_growth, adj_angles)
]
data_frame2 = np.stack((xdata_frame2, ydata_frame2), axis=1)

shifted_points_x = [xdata_frame2[i] for i in shifted_points_idx]
shifted_points_y = [ydata_frame2[i] for i in shifted_points_idx]
shifted_points = np.stack((shifted_points_x, shifted_points_y), axis=1)

# estimate the clustering behavior
GenerateClusters(shifted_points)

fig = go.Figure(
    go.Scatter(x=xdata_initial, y=ydata_initial, mode="markers", marker_color="red")
)
fig.update_layout(width=800, height=800, plot_bgcolor="White")
fig.update_xaxes(
    tickvals=[],
    range=[0, 100],
    showline=True,
    linewidth=2,
    linecolor="black",
    mirror=True,
)
fig.update_yaxes(
    tickvals=[],
    range=[0, 100],
    showline=True,
    linewidth=2,
    linecolor="black",
    mirror=True,
)

fig.add_shape(
    type="circle",
    xref="x",
    yref="y",
    x0=centerpoint[0] - grain_radius_initial,
    y0=centerpoint[1] - grain_radius_initial,
    x1=centerpoint[0] + grain_radius_initial,
    y1=centerpoint[1] + grain_radius_initial,
)

fig.show()

fig2 = go.Figure(
    go.Scatter(x=xdata_frame2, y=ydata_frame2, mode="markers", marker_color="red")
)
fig2.update_layout(width=800, height=800, plot_bgcolor="White")
fig2.update_xaxes(
    tickvals=[],
    range=[0, 100],
    showline=True,
    linewidth=2,
    linecolor="black",
    mirror=True,
)
fig2.update_yaxes(
    tickvals=[],
    range=[0, 100],
    showline=True,
    linewidth=2,
    linecolor="black",
    mirror=True,
)

fig2.add_shape(
    type="circle",
    xref="x",
    yref="y",
    x0=centerpoint[0] - grain_radius_final,
    y0=centerpoint[1] - grain_radius_final,
    x1=centerpoint[0] + grain_radius_final,
    y1=centerpoint[1] + grain_radius_final,
)

fig2.show()
