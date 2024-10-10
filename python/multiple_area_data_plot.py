#! /usr/bin/env python3

import argparse

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, Input, Output, State, dcc, html


def addLabel(trace, points, selector):
    pass


def simpleFit(x: pd.Series, y: pd.Series, threshold: int = 100):
    """Fit a linear function to either the whole dataset,
    or until threshold is reached. Returns the model function.
    x: Pandas Series
    y: Pandas Series
    threshold: (optional) positive integer"""
    if threshold < 0:
        raise ValueError("Integer value must be greater than 0")
    if y.iloc[-1] > threshold:
        model = np.polyfit(x, y, 1)
        last_point_index = len(y) - 1
    else:
        last_point_index = next(
            (idx for idx, obj in enumerate(y) if obj < threshold), len(y) - 1
        )
        model = np.polyfit(x[:last_point_index], y[:last_point_index], 1)
    return (np.poly1d(model), last_point_index)


parser = argparse.ArgumentParser(
    usage="%(prog)s [-h] file [file file ...]",
    description="Plot the specified area vs time plots on a single graph",
)
parser.add_argument("files", metavar="file", nargs="+", help="The files to plot")
parser.add_argument("--names", nargs="+", help="Legend entries for the plotted files")

args = parser.parse_args()

if not args.names or not (len(args.names) == len(args.files)):
    args.names = args.files

colors = [
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#666666",
]

# Selected based on the above 'colors' values at 50% shade using
# https://imagecolorpicker.com/color-code
shades = [
    "#0E4F3C",
    "#6D3001",
    "#3B385A",
    "#741545",
    "#33530F",
    "#735601",
    "#533B0F",
    "#333333",
]

# names = [
#     "<100> 20\N{DEGREE SIGN} U3Mo at 1200 K",
#     "<111> 30\N{DEGREE SIGN} U10Mo at 1300 K",
#     "<110> 45\N{DEGREE SIGN} U21Mo at 1150 K",
# ]

ddtype = {
    "time": "float",
    "grain1_area": "float",
    "grain2_area": "float",
}
data = [
    pd.read_csv(file, sep=" ", names=tuple(ddtype.keys()), dtype=ddtype, skiprows=[0])
    for file in args.files
]

font_size = 26
fig = go.Figure()
fig.update_xaxes(
    title="Time (ps)",
    range=[-100, 4100],
    title_font_size=font_size,
    tickfont_size=font_size,
    ticks="inside",
    showline=True,
    linecolor="#000000",
    linewidth=3,
    ticklen=8,
    tickwidth=2,
    mirror="all",
    showgrid=False,
).update_yaxes(
    title="Area (\N{Latin Capital Letter A With Ring Above}<sup>2</sup>)",
    range=[-1000, 33000],
    title_font_size=font_size,
    tickfont_size=font_size,
    ticks="inside",
    showline=True,
    linecolor="#000000",
    linewidth=3,
    ticklen=8,
    tickwidth=2,
    mirror="all",
    showgrid=False,
).update_layout(
    width=800,
    height=800,
    margin=dict(l=5, r=10, b=50, t=50, pad=4),
    legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.58, font_size=font_size),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

for i, df in enumerate(data):
    fit, last_point_index = simpleFit(df["time"], df["grain2_area"])
    fig.add_trace(
        go.Scatter(
            x=df["time"],
            y=df["grain2_area"],
            mode="markers",
            marker_color=colors[i],
            name=args.names[i],
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            x=df["time"][:last_point_index],
            y=fit(df["time"[:last_point_index]]),
            mode="lines",
            line=dict(color=shades[i], dash="20,30", width=2.5),
            showlegend=False,
        )
    )

app = Dash()

app.layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Label("Label:   "),
                        dcc.Input(
                            id="input-annotations", type="text", style={"width": "60%"}
                        ),
                        html.Button(
                            id="btn-submit", type="submit", children="Add annotation"
                        ),
                    ]
                )
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dcc.Graph(
                            id="area-vs-time-graph",
                            figure=fig,
                            config={
                                "editable": True,
                                "edits": {
                                    "shapePosition": True,
                                    "annotationPosition": True,
                                },
                            },
                        )
                    ]
                )
            ]
        ),
    ]
)


@app.callback(
    Output(component_id="area-vs-time-graph", component_property="figure"),
    Input(component_id="area-vs-time-graph", component_property="figure"),
    State(component_id="input-annotations", component_property="value"),
    Input("btn-submit", "n_clicks"),
)
def updateFigure(figure, annotations, n_clicks):
    fig1 = go.Figure(figure)

    try:
        if len(annotations) != 0:
            fig1.add_annotation(
                x=2000, y=1000, text=annotations, showarrow=False, font_size=font_size
            )
        else:
            return fig
    except TypeError:
        return fig
    return fig1


# fig.add_annotation(x=2800, y=5000, text="(a)", showarrow=False, font=dict(size=26))
# fig.add_annotation(x=3000, y=15_000, text="(b)", showarrow=False, font=dict(size=26))
# fig.add_annotation(x=3500, y=24_500, text="(c)", showarrow=False, font=dict(size=26))
# fig.add_annotation(x=5000, y=21_000, text="(d)", showarrow=False, font=dict(size=26))

# fig.show()
# fig.write_image(
#     "/media/jarinf/Research1/Images/"
#     "Manuscript Images/UMoXe paper/UMo_area_vs_time_examples_no_legend.png"
# )

if __name__ == "__main__":
    app.run_server(debug=True)
