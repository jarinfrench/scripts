#! /usr/bin/env python
"""This script counts the percentages of impurities within 2nn
distance to generate an 'averaged' histogram over all temperatures,
and comparing all concentrations"""

from collections import defaultdict
from os import getcwd
import pandas as pd
import plotly.graph_objects as go

Mo_imp = [f"{i}at%" for i in [1, 3, 5, 7, 10, 12, 14, 17, 19, 21]]
Xe_imp = [f"{i:.2f}at%" for i in [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.0, 2.5, 5, 7.5]]
Temperatures = [1050, 1200, 1400]
files = [f"dir_{i}/nearest_neighbors.txt" for i in range(1, 4)]

base_dir = getcwd().split("/")[-1]
if base_dir == "Mo":
    impurity = Mo_imp
elif base_dir == "Xe":
    impurity = Xe_imp
else:
    raise Exception(
        f"Cannot identify impurity type from ./{base_dir}" + "(expected 'Mo' or 'Xe')"
    )

avg_counts = {}
for imp in impurity:
    for T in Temperatures:
        df = pd.concat(
            (
                pd.read_csv(
                    "/".join([imp, "T" + str(T), file]),
                    sep=" ",
                    names=("step", "id", "n_impurities", "n_neighbors"),
                    dtype={
                        "step": "int",
                        "id": "int",
                        "n_impurities": "int",
                        "n_neighbors": "int",
                    },
                ).assign(
                    temperature=T, concentration=float(imp.split("a")[0]), run=i + 1
                )
                for i, file in enumerate(files)
            ),
            ignore_index=True,
            axis=0,
        )
        df = df.set_index(["concentration", "temperature", "run", "step"])
        df = df.sort_index()
        keys = df.index.unique()
        current_c = df.index.get_level_values(0).unique()[0]
        current_T = df.index.get_level_values(1).unique()[0]
        number_of_impurities = sorted(df.n_impurities.unique())
        nn_counts = {}
        nn_counts = {
            key: {
                n_imps: df.loc[key][df.loc[key]["n_impurities"] == n_imps].size
                for n_imps in number_of_impurities
            }
            for key in keys
        }
        d = defaultdict(lambda: 0)
        for key, d1 in nn_counts.items():
            for n_neigh, count in d1.items():
                d[n_neigh] += count
        avg_counts[current_c, current_T] = {
            k: v / len(nn_counts.keys()) for k, v in d.items()
        }

avg_data = pd.DataFrame.from_dict(avg_counts).sort_index()
ws = "       "
x_labels = [
    f"{ws}{i} Impurity" if i == 1 else f"{ws}{i} Impurities"
    for i in avg_data.index.values
]
# Modified from https://community.plotly.com/t/grouped-stacked-bar-chart/60805/5
fig = go.Figure(
    layout=go.Layout(
        barmode="relative",
        yaxis2=go.layout.YAxis(visible=False, matches="y", overlaying="y", anchor="x"),
        yaxis3=go.layout.YAxis(visible=False, matches="y", overlaying="y", anchor="x"),
        legend_x=0.86,
        legend_y=1,
        legend_orientation="h",
        # legend_traceorder="reversed",
        hovermode="x unified",
        margin=dict(b=0, t=10, l=0, r=10),
    )
)

# chosen with the help of https://coolors.co/1b9e77-d95f02-7570b3-f9f9f9-4e3d42
colors = {
    impurity[0]: {
        Temperatures[0]: "#06231A",
        Temperatures[1]: "#3D1B01",
        Temperatures[2]: "#1E1C35",
    },
    impurity[1]: {
        Temperatures[0]: "#0C4634",
        Temperatures[1]: "#652C01",
        Temperatures[2]: "#2D2A50",
    },
    impurity[2]: {
        Temperatures[0]: "#12694F",
        Temperatures[1]: "#8D3E01",
        Temperatures[2]: "#3C386B",
    },
    impurity[3]: {
        Temperatures[0]: "#188C69",
        Temperatures[1]: "#B65002",
        Temperatures[2]: "#4B4686",
    },
    impurity[4]: {
        Temperatures[0]: "#1EAE83",
        Temperatures[1]: "#D95F02",
        Temperatures[2]: "#5A54A0",
    },
    impurity[5]: {
        Temperatures[0]: "#23D19D",
        Temperatures[1]: "#FD750D",
        Temperatures[2]: "#716CB2",
    },
    impurity[6]: {
        Temperatures[0]: "#3FDEAF",
        Temperatures[1]: "#FD8C35",
        Temperatures[2]: "#8B87C0",
    },
    impurity[7]: {
        Temperatures[0]: "#62E4BD",
        Temperatures[1]: "#FDA35D",
        Temperatures[2]: "#A4A1CE",
    },
    impurity[8]: {
        Temperatures[0]: "#85EACC",
        Temperatures[1]: "#FEBA86",
        Temperatures[2]: "#BEBCDC",
    },
    impurity[9]: {
        Temperatures[0]: "#A8F0DA",
        Temperatures[1]: "#FED1AE",
        Temperatures[2]: "#D8D7EA",
    },
}
for i, imp in enumerate(colors):
    imp_val = float(imp.split("a")[0])
    for j, temp in enumerate(avg_data[imp_val].columns):
        if (avg_data[imp_val][temp] == 0).all():
            continue
        fig.add_bar(
            x=x_labels,
            y=avg_data[imp_val][temp],
            yaxis=f"y{j+1}",
            offsetgroup=str(j),
            offset=(j - 1) * 1 / 4,
            width=1 / 4,
            legendgroup=temp,
            legendgrouptitle_text=temp,
            legendrank=1 / imp_val,
            name=imp,
            marker_color=colors[imp][temp],
            marker_line=dict(width=2, color="#333"),
            hovertemplate="%{y}<extra></extra>",
        )

fig.update_layout(
    yaxis=dict(
        title="Average # of 2nn Spheres across all steps, runs, and concentrations",
        type="log",
        tickvals=[],
    )
)

fig.write_html("nearest_neighbor_histogram_T_and_c.html")
fig.write_image("nearest_neighbor_histogran_T_and_c.png")
