#! /usr/bin/env python3

from dash import Dash, html, dcc, Input, Output
import plotly.graph_objects as go
from plotly.graph_objs.layout import YAxis, XAxis, Margin
import plotly.express as px
import pandas as pd
import numpy as np
from copy import deepcopy
from itertools import product
from numbers import Number

color_cycle = px.colors.qualitative.Dark24
symbol_dict = {'None': 'circle', 'Mo': 'square', 'Xe': 'triangle-up', 'MoXe': 'star'}

kB = 8.617e-5
def Arrhenius(val):
    if isinstance(val, list):
        return [1 / (kB * i) for i in val if isinstance(i, Number)]
    elif isinstance(val, Number):
        return 1/ (kB * val)
    else:
        raise TypeError("Invalid type: must be a numeric type")

def update_selection(dataframe, match_criteria, key):
    tmp = [None] * len(match_criteria)

    for idx, val in enumerate(match_criteria):
        tmp[idx] = (dataframe[key] == val)

    return [any(i) for i in zip(*tmp)]

ddtype = {"type": "str", "c": "float", "axis": "int", "theta": "float", "T": "int", "growth": "float", "M": "float", "c_label": "str"}
data = pd.read_csv("/media/jarinf/Research2/tmp/U/raw_individual_mobility_data_with_percent_growth.txt", sep = ',',
       names = tuple(ddtype.keys()), dtype = ddtype, skiprows = 1
       ).set_index(['type', 'axis', 'c_label', 'theta', 'c', 'growth'])

# Generate the app
# Defaults and constants
default_fig = go.Figure()
default_fig.update_layout(autosize = False, width = 1500, height = 500, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
ytickvals = [i*1e-10 for i in range(1,10)] + [i * 1e-9 for i in range(10)] + [i *1e-8 for i in range(10)] + [1e-7]
yticktext = [str(i) if i in [1e-10,1e-9,1e-8,1e-7] else "" for i in ytickvals]
xtickvals = [Arrhenius(i) for i in np.sort(data['T'].unique())[::-1]]
xticktext = [str(i) for i in np.sort(data['T'].unique())[::-1]]
default_fig.update_yaxes(title = r"Reduced Mobility (m<sup>2</sup>/s)", exponentformat = 'e', type = 'log', range = [-10.5, -7], tickvals = ytickvals, ticktext = yticktext)
default_fig.update_xaxes(title = "Temperature (K)", range = [Arrhenius(1450), Arrhenius(1000)], tickvals = xtickvals, ticktext = xticktext)

system_select_dict = [{'label': i, 'value': i} if not i == "None" else {'label': 'Pure U', 'value': i} for i in data.index.get_level_values('type').unique()]
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
    {"label": "5Mo0.75Xe", "value": "5Mo0.75Xe", "disabled": True}
]
none_labels = ["None"]
Mo_labels = [f"{i}Mo" for i in [1,3,5,10,14,17,21]]
Xe_labels = [f"{i:.2f}Xe" for i in [0.25,0.50,0.75,1.00,1.25,1.50]]
MoXe_labels = [f"{i}Mo{j:.2f}Xe" for i in [1,3,5] for j in [0.25, 0.50, 0.75]]
conc_labels = {"None": none_labels, "Mo": Mo_labels, "Xe": Xe_labels, "MoXe": MoXe_labels}

# The actual app
app = Dash(__name__)
app.layout = html.Div([
    html.H1('Raw Mobility Calculations in the UMoXe System', style = {'margin': 'auto', 'width': '50%'}),
    html.H2('Select the data to plot'),
    html.Div([
        html.Div([
            html.Div([html.Label('Select the system(s):'), dcc.Checklist(system_select_dict, id = 'impurity-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px', 'width': '6em'})], style = {'width': "280px", "float": "left"}),
            html.Div([html.Label('Select the rotation axis/axes:'), dcc.Checklist(data.index.get_level_values('axis').unique(), id = 'axis-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px'})], style = {'width': '280px', 'float': 'left', 'padding-left': '15px'}),
            html.Div([html.Label('Select the misorientation(s):'), dcc.Checklist(data.index.get_level_values('theta').unique(), id = 'misorientation-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px'})], style = {'width': "280px", "float": "left", 'padding-left': '15px'}),
            html.Div([html.Label('Select the concentration(s):'), dcc.Checklist(conc_dropdown_dict, id = 'concentration-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px', 'width': '8em', 'border': '1px green'})], style = {'width': "450px", "float": "left", 'padding-left': '15px'}),
            html.Div([html.Label('Select the temperature(s):'), dcc.Checklist(data["T"].unique(), id = 'temperature-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px'})], style = {'width': "280px", "float": "left", 'padding-left': '15px'})
            ]),
    html.Div([
        html.Label('Growth Threshold:'),
        dcc.Input(
        id = 'growth-threshold-input',
        type = 'number',
        placeholder = 'Enter the minimum % growth as a fraction',
        debounce = False,
        min = 0,
        max = 1, step = 0.01,
        style = {'width': '260px'}
        )], style = {'width': '1000px', 'float': 'left'}),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Br(),
    html.Div([
        html.Label('Closed symbols: simulation meets growth threshold')
    ]),
    html.Div([
        html.Label('Open symbols: simulation does not meet growth threshold')
    ]),
    html.Div([
        html.Label('If showing the average, the average value uses all available data (full opacity). The weaker opacity data point(s) are the averages of the respective subset(s).')
    ]),
    html.Div([
    # See https://stackoverflow.com/questions/66301881/plotly-plotting-one-trace-with-two-x-axes-that-are-not-linearly-related
    # for some discussion on non-linearly linked x axes in Plotly
    # See also:
    #   https://stackoverflow.com/questions/61803438/python-plotly-figure-with-secondary-x-axis-linked-to-primary
    #   https://stackoverflow.com/questions/27334585/in-plotly-how-do-i-create-a-linked-x-axis
        dcc.Graph(figure = default_fig, id = 'mobility-plot'),
        dcc.RadioItems(
            options = [{'label': 'Show all', 'value': False},
            {'label': 'Show average only', 'value': True}],
            value = False,
            inline = True,
            id = 'average-mobility-radio')
    ])
    ]),
])

@app.callback(
    Output('concentration-checklist', 'options'),
    Input('impurity-checklist','value')
)
def update_available_concentrations(impurity):
    concentration_options = deepcopy(conc_dropdown_dict)
    if not impurity:
        for i, ii in [(j, item) for j, item in enumerate(concentration_options)]:
            concentration_options[i]["disabled"] = True
    else:
        allowed_labels = []
        if 'None' in impurity:
            allowed_labels += none_labels
        if 'Mo' in impurity:
            allowed_labels += Mo_labels
        if 'Xe' in impurity:
            allowed_labels += Xe_labels
        if 'MoXe' in impurity:
            allowed_labels += MoXe_labels
        for i, ii in [(j, item) for j, item in enumerate(concentration_options)]:
            if ii['label'] in allowed_labels:
                concentration_options[i]["disabled"] = False
    return concentration_options

@app.callback(
    Output('mobility-plot', 'figure'),
    Input('impurity-checklist', 'value'),
    Input('axis-checklist', 'value'),
    Input('misorientation-checklist', 'value'),
    Input('concentration-checklist', 'value'),
    Input('temperature-checklist', 'value'),
    Input('growth-threshold-input', 'value'),
    Input('average-mobility-radio', 'value')
)
def update_figure(impurities, axes, misorientations, concentrations, temperatures,
                  threshold, show_average):
    if not impurities or not axes or not misorientations or not concentrations or not temperatures or not threshold:
        return default_fig
    keys = []
    if 'None' in impurities:
        keys += list(product(*[['None'], axes, [i for i in concentrations if i == '0Mo0.00Xe'], misorientations, temperatures, [threshold]]))
    if 'Mo' in impurities:
        keys += list(product(*[['Mo'], axes, [i for i in concentrations if not i[0] == '0' and i[-6:] == '0.00Xe'], misorientations, temperatures, [threshold]]))
    if 'Xe' in impurities:
        keys += list(product(*[['Xe'], axes, [i for i in concentrations if i[0:3] == '0Mo' and not i[-6:] == '0.00Xe'], misorientations, temperatures, [threshold]]))
    if 'MoXe' in impurities:
        keys += list(product(*[['MoXe'], axes, [i for i in concentrations if not i[0:3] == '0Mo' and not i[-6:] == '0.00Xe'], misorientations, temperatures, [threshold]]))

    Ts = []
    Ms = []
    imps = []
    labels = []
    for key in keys:
        Ttmp = [[],[]]
        Mtmp = [[],[]]
        res = data.query("type == @key[0] and axis == @key[1] and c_label == @key[2] and theta == @key[3] and T == @key[4]")
        if res.empty:
            continue
        enough = res.query("growth >= @key[5]")
        not_enough = res.query("growth < @key[5]")
        Ttmp[0] = enough['T'].to_list()
        Ttmp[1] = not_enough['T'].to_list()
        Mtmp[0] = enough['M'].to_list()
        Mtmp[1] = not_enough['M'].to_list()
        Ts.append(Ttmp)
        Ms.append(Mtmp)
        imps.append(key[0])
        simplified_label = next((sub for sub in conc_dropdown_dict if sub['value'] == key[2]), None)['label']
        if simplified_label == "None":
            simplified_label = ''
        labels.append(f"U{simplified_label} {key[1]} {key[3]}degree {key[4]}K")

    if not Ts or not Ms:
        return default_fig

    fig = deepcopy(default_fig)
    color_idx_dict = {'None': 0, 'Mo': 0, 'Xe': 0, 'MoXe': 0}
    for idx, x_y_imp_label in enumerate(zip(Ts, Ms, imps, labels)):
        solidx = Arrhenius(x_y_imp_label[0][0])
        solidy = [abs(v) for v in x_y_imp_label[1][0]]
        openx = Arrhenius(x_y_imp_label[0][1])
        openy = [abs(v) for v in x_y_imp_label[1][1]]
        imp = x_y_imp_label[2]
        label = x_y_imp_label[3]
        color = color_cycle[color_idx_dict[imp] % len(color_cycle)]
        color_idx_dict[imp] += 1
        if show_average:
            solidx = np.mean(solidx) if solidx else None
            solidy = np.mean(solidy) if solidy else None
            openx = np.mean(openx) if openx else None
            openy = np.mean(openy) if openy else None
            if solidx and openx:
                all_x = np.mean([solidx, openx])
                all_y = np.mean([solidy, openy])
                avg_symbol = symbol_dict[imp]
            elif solidx and not openx:
                all_x = solidx
                all_y = solidy
                avg_symbol = symbol_dict[imp]
            else:
                all_x = openx
                all_y = openy
                avg_symbol = symbol_dict[imp] + '-open'
            fig.add_trace(go.Scatter(x = [all_x], y = [all_y], name = label, mode = 'markers', showlegend = True, marker_color = color, marker_symbol = avg_symbol))
            if solidx and openx:
                fig.add_trace(go.Scatter(x = [solidx], y = [solidy], name = label, showlegend = False, mode = 'markers', marker_color = color, marker_symbol = symbol_dict[imp], opacity = 0.5))
                fig.add_trace(go.Scatter(x = [openx], y = [openy], name = label, showlegend = False, mode = 'markers', marker_color = color, marker_symbol = symbol_dict[imp] + '-open', opacity = 0.5))
        else:
            fig.add_trace(go.Scatter(x = solidx, y = solidy, name = label, showlegend = True, mode = 'markers', marker_color = color, marker_symbol = symbol_dict[imp]))
            fig.add_trace(go.Scatter(x = openx, y = openy, name = label, showlegend = False if solidx else True, mode = 'markers', marker_color = color, marker_symbol = symbol_dict[imp] + '-open'))

    return fig

if __name__ == "__main__":
    app.run_server(debug = True)
