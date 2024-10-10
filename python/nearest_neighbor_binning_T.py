#! /usr/bin/env python
"""This script counts the percentages of impurities within 2nn
distance to generate an 'averaged' histogram over all temperatures"""

from itertools import product
from collections import defaultdict
import pandas as pd
import plotly.graph_objects as go

Temperatures = [1050, 1200, 1400]
files = [f"dir_{i}/nearest_neighbors.txt" for i in range(1, 4)]
all_files = ["/".join(["T" + str(i), j]) for i, j in product(Temperatures, files)]
df = pd.concat(
    (
        pd.read_csv(
            file,
            sep=" ",
            names=("timestep", "id", "n_impurities", "n_neighbors"),
            dtype={
                "timestep": "int",
                "id": "int",
                "n_impurities": "int",
                "n_neighbors": "int",
            },
        ).assign(temperature=Temperatures[int(i / 3) - 1], run=(i % 3) + 1)
        for i, file in enumerate(all_files)
    ),
    ignore_index=True,
    axis=0,
)

nn_counts = {}
for T in df.temperature.unique():
    nn_counts[T] = {}
    for run in df.run.unique():
        nn_counts[T][run] = {}
        for t in df.timestep.unique():
            nn_counts[T][run][t] = {}
            for n in df.n_impurities.unique():
                nn_counts[T][run][t][n] = sum(
                    df.query("temperature == @T and run == @run and timestep == @t")[
                        "n_impurities"
                    ]
                    == n
                )

avg_counts = {}
for T, d1 in nn_counts.items():
    d = defaultdict(lambda: 0)
    for run, d2 in d1.items():
        for timestep, d3 in d2.items():
            for n_neigh, count in d3.items():
                d[n_neigh] += count
    avg_counts[T] = {k: v / sum(len(d1[i]) for i in d1) for k, v in d.items()}
avg_data = pd.DataFrame.from_dict(avg_counts).sort_index()
avg_data = avg_data.reindex(sorted(avg_data.columns), axis=1)
x_labels = [
    f"{i} Impurity" if i == 1 else f"{i} Impurities" for i in avg_data.index.values
]
fig = go.Figure()
for i in avg_data.keys():
    fig.add_trace(go.Bar(x=x_labels, y=avg_data[i], name=f"{i} K"))
fig.update_layout(
    yaxis=dict(
        title="Average number of 2nn Spheres across all timesteps and runs", type="log"
    )
)
# fig.update_xaxes(tickvals=[0, 1/15, 2/15, 3/15, 4/15],
#                  ticktext=['0 impurities', '1 impurity', '2 impurities',
#                            '3 impurities', '4 impurities'], range=(-1/50, 9/30))
# fig.update_layout(barmode = 'overlay')
# fig.update_traces(opacity = 0.25)
# fig.show()
fig.write_image("nearest_neighbor_histogram_T.png")
