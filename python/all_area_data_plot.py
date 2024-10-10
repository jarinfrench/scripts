# /usr/bin/env python3

from dash import Dash, html, dcc, Input, Output
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from copy import deepcopy
from itertools import product


def get_selection(dataframe, match_criteria, key):
    tmp = [None] * len(match_criteria)
    for idx, val in enumerate(match_criteria):
        tmp[idx] = dataframe[key] == val

    return [any(i) for i in zip(*tmp)]


ddtype = {
    "type": "str",
    "axis": "int",
    "c": "float",
    "theta": "float",
    "T": "int",
    "run_index": "int",
    "time": "int",
    "grain1_area": "float",
    "grain2_area": "float",
    "c_label": "str",
}

data = pd.read_csv(
    "/media/jarinf/Research2/tmp/U/area_data_raw_v3_extended_simplified.txt",
    sep=" ",
    names=tuple(ddtype.keys()),
    dtype=ddtype,
    skiprows=1,
).set_index(["type", "axis", "c_label", "theta", "T", "run_index", "c"])

conc_dropdown_dict = [
    {"label": "None", "value": "0Mo0.00Xe", "disabled": True},
    {"label": "1Mo", "value": "1Mo0.00Xe", "disabled": True},
    {"label": "3Mo", "value": "3Mo0.00Xe", "disabled": True},
    {"label": "5Mo", "value": "5Mo0.00Xe", "disabled": True},
    {"label": "10Mo", "value": "10Mo0.00Xe", "disabled": True},
    {"label": "14Mo", "value": "14Mo0.00Xe", "disabled": True},
    {"label": "17Mo", "value": "17Mo0.00Xe", "disabled": True},
    {"label": "21Mo", "value": "21Mo0.00Xe", "disabled": True},
    {"label": "0.25Xe", "value": "0Mo0.25Xe", "disabled": True},
    {"label": "0.50Xe", "value": "0Mo0.50Xe", "disabled": True},
    {"label": "0.75Xe", "value": "0Mo0.75Xe", "disabled": True},
    {"label": "1.00Xe", "value": "0Mo1.00Xe", "disabled": True},
    {"label": "1.25Xe", "value": "0Mo1.25Xe", "disabled": True},
    {"label": "1.50Xe", "value": "0Mo1.50Xe", "disabled": True},
    {"label": "1Mo0.25Xe", "value": "1Mo0.25Xe", "disabled": True},
    {"label": "1Mo0.50Xe", "value": "1Mo0.50Xe", "disabled": True},
    {"label": "1Mo0.75Xe", "value": "1Mo0.75Xe", "disabled": True},
    {"label": "3Mo0.25Xe", "value": "3Mo0.25Xe", "disabled": True},
    {"label": "3Mo0.50Xe", "value": "3Mo0.50Xe", "disabled": True},
    {"label": "3Mo0.75Xe", "value": "3Mo0.75Xe", "disabled": True},
    {"label": "5Mo0.25Xe", "value": "5Mo0.25Xe", "disabled": True},
    {"label": "5Mo0.50Xe", "value": "5Mo0.50Xe", "disabled": True},
    {"label": "5Mo0.75Xe", "value": "5Mo0.75Xe", "disabled": True},
]

none_labels = ["None"]
Mo_labels = [f"{i}Mo" for i in [1, 3, 5, 10, 14, 17, 21]]
Xe_labels = [f"{i:.2f}Xe" for i in [0.25, 0.50, 0.75, 1.00, 1.25, 1.50]]
MoXe_labels = [f"{i}Mo{j:.2f}Xe" for i in [1, 3, 5] for j in [0.25, 0.50, 0.75]]
conc_labels = {
    "None": none_labels,
    "Mo": Mo_labels,
    "Xe": Xe_labels,
    "MoXe": MoXe_labels,
}
system_select_dict = [
    {"label": i, "value": i} if not i == "None" else {"label": "Pure U", "value": i}
    for i in data.reset_index().type.unique()
]
run_index_select_dict = [
    {"label": i, "value": i, "disabled": False}
    for i in np.sort(data.reset_index().run_index.unique())
]

app = Dash(__name__)
default_fig = go.Figure()
default_fig.update_layout(
    autosize=False, width=1500, height=550, margin=dict(l=5, r=10, b=50, t=50, pad=4)
)
ytickvals = [i * 5000 for i in range(15)]
yticktext = [str(i) if i in [i * 10000 for i in range(8)] else "" for i in ytickvals]
default_fig.update_yaxes(
    title="Area (\U0000212b<sup>2</sup>)",
    range=[-10, 70000],
    tickvals=ytickvals,
    ticktext=yticktext,
)
default_fig.update_xaxes(title="Time (ps)", showgrid=False, range=[0, 8000])

app.layout = html.Div(
    [
        html.H1(
            "Grain Size vs Time in the UMoXe System",
            style={"margin": "auto", "width": "50%"},
        ),
        html.H2("Select the data to plot"),
        html.Div(
            [
                html.Div(
                    [
                        html.Label("System:"),
                        dcc.Checklist(
                            system_select_dict,
                            id="impurity-checklist",
                            labelStyle={
                                "display": "inline-block",
                                "padding-left": "10px",
                                "width": "6em",
                            },
                        ),
                    ],
                    style={"width": "280px", "float": "left"},
                ),
                html.Div(
                    [
                        html.Label("Rotation axis:"),
                        dcc.Checklist(
                            np.sort(data.reset_index().axis.unique()),
                            id="axis-checklist",
                            labelStyle={
                                "display": "inline-block",
                                "padding-left": "10px",
                            },
                        ),
                    ],
                    style={"width": "280px", "float": "left", "padding-left": "5px"},
                ),
                html.Div(
                    [
                        html.Label("Misorientation:"),
                        dcc.Checklist(
                            np.sort(data.reset_index().theta.unique()),
                            id="misorientation-checklist",
                            labelStyle={"display": "inline-block"},
                        ),
                    ],
                    style={"width": "280px", "float": "left", "padding-left": "5px"},
                ),
                html.Div(
                    [
                        html.Label("Concentration:"),
                        dcc.Checklist(
                            conc_dropdown_dict,
                            id="concentration-checklist",
                            labelStyle={
                                "display": "inline-block",
                                "width": "8em",
                                "border": "1px green",
                            },
                        ),
                    ],
                    style={"width": "450px", "float": "left", "padding-left": "5px"},
                ),
                html.Div(
                    [
                        html.Label("Temperature:"),
                        dcc.Checklist(
                            np.sort(data.reset_index()["T"].unique()),
                            id="temperature-checklist",
                            labelStyle={
                                "display": "inline-block",
                                "padding-left": "5px",
                            },
                        ),
                    ],
                    style={"width": "240px", "float": "left", "padding-left": "5px"},
                ),
                html.Div(
                    [
                        html.Label("Run #:"),
                        dcc.Checklist(
                            run_index_select_dict,
                            id="run-index-checklist",
                            labelStyle={"display": "inline-block"},
                        ),
                    ],
                    style={"width": "200px", "float": "left", "padding-left": "5px"},
                ),
            ]
        ),
        html.Div(
            dcc.Graph(figure=default_fig, id="area-data-graph"),
            style={"width": "80%", "float": "left"},
        ),
    ]
)


@app.callback(
    Output("concentration-checklist", "options"), Input("impurity-checklist", "value")
)
def update_available_concentrations(impurity):
    concentration_options = deepcopy(conc_dropdown_dict)
    if not impurity:
        for i, ii in [(j, item) for j, item in enumerate(concentration_options)]:
            concentration_options[i]["disabled"] = True
    else:
        allowed_labels = []
        if "None" in impurity:
            allowed_labels += none_labels
        if "Mo" in impurity:
            allowed_labels += Mo_labels
        if "Xe" in impurity:
            allowed_labels += Xe_labels
        if "MoXe" in impurity:
            allowed_labels += MoXe_labels
        for i, ii in [(j, item) for j, item in enumerate(concentration_options)]:
            if ii["label"] in allowed_labels:
                concentration_options[i]["disabled"] = False
    return concentration_options


@app.callback(
    Output("area-data-graph", "figure"),
    Input("impurity-checklist", "value"),
    Input("axis-checklist", "value"),
    Input("misorientation-checklist", "value"),
    Input("concentration-checklist", "value"),
    Input("temperature-checklist", "value"),
    Input("run-index-checklist", "value"),
)
def update_figure(
    impurities, axes, misorientations, concentrations, temperatures, run_indices
):
    if (
        not impurities
        or not axes
        or not misorientations
        or not concentrations
        or not temperatures
        or not run_indices
    ):
        return default_fig
    keys = []
    if "None" in impurities:
        keys += list(
            product(
                *[
                    ["None"],
                    axes,
                    [i for i in concentrations if i == "0Mo0.00Xe"],
                    misorientations,
                    temperatures,
                    run_indices,
                ]
            )
        )
    if "Mo" in impurities:
        keys += list(
            product(
                *[
                    ["Mo"],
                    axes,
                    [
                        i
                        for i in concentrations
                        if not i[0] == "0" and i[-6:] == "0.00Xe"
                    ],
                    misorientations,
                    temperatures,
                    run_indices,
                ]
            )
        )
    if "Xe" in impurities:
        keys += list(
            product(
                *[
                    ["Xe"],
                    axes,
                    [
                        i
                        for i in concentrations
                        if i[0:3] == "0Mo" and not i[-6:] == "0.00Xe"
                    ],
                    misorientations,
                    temperatures,
                    run_indices,
                ]
            )
        )
    if "MoXe" in impurities:
        keys += list(
            product(
                *[
                    ["MoXe"],
                    axes,
                    [
                        i
                        for i in concentrations
                        if not i[0:3] == "0Mo" and not i[-6:] == "0.00Xe"
                    ],
                    misorientations,
                    temperatures,
                    run_indices,
                ]
            )
        )

    times = []
    areas = []
    labels = []
    for key in keys:
        res = data.query(
            "type == @key[0] and axis == @key[1] and c_label == @key[2] and "
            + "theta == @key[3] and T == @key[4] and run_index == @key[5]"
        )
        if res.empty:
            continue
        times.append(res.time.to_list())
        areas.append(res.grain2_area.to_list())
        simplified_label = next(
            (sub for sub in conc_dropdown_dict if sub["value"] == key[2]), None
        )["label"]
        if simplified_label == "None":
            simplified_label = ""
        labels.append(
            f"U{simplified_label} {key[1]} {key[3]}degree {key[4]}K run{key[5]}"
        )

    if not times or not areas:
        return default_fig
    fig = deepcopy(default_fig)
    for x, y, l in zip(times, areas, labels):
        fig.add_trace(go.Scatter(x=x, y=y, name=l, mode="markers"))
    fig.update_yaxes(range=[-10, 1.2 * max([max(i) for i in areas])])
    return fig


# This would work if it did things in a reasonable amount of time, but currently
# it takes a very long time to update just once.
# @app.callback(
#     Output("run-index-checklist", "options"),
#     Input("impurity-checklist", "value"),
#     Input("axis-checklist", "value"),
#     Input("misorientation-checklist", "value"),
#     Input("concentration-checklist", "value"),
#     Input("temperature-checklist", "value"),
# )
# def update_available_run_indices(
#     impurities, axes, misorientations, concentrations, temperatures
# ):
#     if (
#         not impurities
#         or not axes
#         or not misorientations
#         or not concentrations
#         or not temperatures
#     ):
#         return run_index_select_dict
#     selection = [True] * len(data)
#     selection = [
#         all(i) for i in zip(selection, get_selection(data, impurities, "type"))
#     ]
#     print(selection)
#     selection = [all(i) for i in zip(selection, get_selection(data, axes, "axis"))]
#     selection = [
#         all(i) for i in zip(selection, get_selection(data, misorientations, "theta"))
#     ]
#     selection = [
#         all(i) for i in zip(selection, get_selection(data, concentrations, "c"))
#     ]
#     selection =
# [all(i) for i in zip(selection, get_selection(data, temperatures, "T"))]

#     selected_data = data[selection]
#     print(selected_data.run_index.unique())
#     return [
#         {"label": i, "value": i, "disabled": False}
#         if i in selected_data.run_index.unique()
#         else {"label": i, "value": i, "disabled": True}
#         for i in np.sort(data.run_index.unique())
#     ]


if __name__ == "__main__":
    app.run_server(debug=False, host="0.0.0.0")
