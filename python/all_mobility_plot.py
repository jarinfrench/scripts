# /usr/bin/env python3

from dash import Dash, html, dcc, Input, Output
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from copy import deepcopy

kB = 8.617e-5

def data_as_dict(df):
    d = df.to_dict()
    names = list(df.index.names)
    organize_by = names[-1]
    res = dict()

    for k1 in ["M", "StdErr"]:
        for keys, val in d[k1].items():
            t = keys[names.index('type')] # impurity type
            axis = keys[names.index('axis')] # rotation axis
            c = keys[names.index('c')] # concentration
            mis = keys[names.index('theta')] # misorientation angle (degrees)
            T = keys[names.index('T')] # temperature (K)
            organize_val = keys[names.index(organize_by)]

            if organize_by == "T": # Arrhenius
                tuplekey = (t,axis,c,mis)
            elif organize_by == "theta": # Misorientation
                tuplekey = (t,axis,T,c)
            elif organize_by == "c": # concentration
                tuplekey = (t,axis,mis,T)
            else:
                raise KeyError(f"Specified key {organize_by} not valid")

            if tuplekey not in res.keys():
                res[tuplekey] = {}
            if organize_by not in res[tuplekey].keys():
                res[tuplekey][organize_by] = []
            if k1 not in res[tuplekey].keys():
                res[tuplekey][k1] = []
            res[tuplekey][k1].append(val)
            if k1 == "M":
                res[tuplekey][organize_by].append(organize_val)

    return res, names

def convertToArrhenius(T):
    return [1 / (kB * t) for t in T]

def update_selection(dataframe, match_criteria, key):
    tmp = [None] * len(match_criteria)

    for idx, val in enumerate(match_criteria):
        tmp[idx] = (dataframe[key] == val)

    return [any(i) for i in zip(*tmp)]

def generateLabel(key, plot_type):
    if plot_type == "Arrhenius":
        if key[0] == "None":
            return f"Pure U <{key[1]}> {key[3]}\N{DEGREE SIGN}"
        else:
            if key[0] == "MoXe":
                c = key[2] * 100
                cMo = int(c)
                cXe = c - cMo
                return f"U{cMo}Mo{cXe}Xe <{key[1]}> {key[3]}\N{DEGREE SIGN}"
            else:
                c = int(key[2] * 100) if key[0] == "Mo" else key[2] * 100
                return f"U{c}{key[0]} <{key[1]}> {key[3]}\N{DEGREE SIGN}"
    elif plot_type == "Concentration":
        if key[0] == "None":
            return f"Pure U <{key[1]}> {key[2]}\N{DEGREE SIGN} T{key[3]}"
        else:
            return f"U{key[0]} <{key[1]}> {key[2]}\N{DEGREE SIGN} T{key[3]}"
    elif plot_type == "Misorientation":
        if key[0] == "None":
            return f"Pure U <{key[1]}> T{key[2]}"
        else:
            if key[0] == "MoXe":
                c = key[3] * 100
                cMo = int(c)
                cXe = c - cMo
                return f"U{cMo}Mo{cXe}Xe <{key[1]}> T{key[2]}"
            else:
                c = int(key[2] * 100) if key[0] == "Mo" else key[2] * 100
                return f"U${c}{key[0]} <{key[1]}> T{key[2]}"


ddtype = {"type": "str", "c": "float", "axis": "int", "theta": "float", "T": "int", "M": "float", "StdErr": "float"}

data = pd.read_csv("/media/jarinf/Research2/tmp/U/mobility_raw_data_v2.txt", sep = ' ',
       names = ("type", "c", "axis", "theta", "T", "M", "StdErr"), dtype = ddtype, skiprows = 1
       )

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
c_label = [None] * len(data)
for i in range(len(data)):
    if data.iloc[i]["type"] == "None":
        label = "None"
    elif data.iloc[i]["type"] == "Mo":
        label = f"{data.iloc[i]['c']*100:.0f}Mo"
    elif data.iloc[i]["type"] == "Xe":
        label = f"{data.iloc[i]['c']*100:.2f}Xe"
    else:
        cMo = int(data.iloc[i]['c']*100)
        cXe = data.iloc[i]['c']*100 - cMo
        label = f"{cMo:.0f}Mo{cXe:.2f}Xe"
    c_label[i] = next((sub["value"] for sub in conc_dropdown_dict if sub["label"] == label))
data["c_label"] = c_label
none_labels = ["None"]
Mo_labels = [f"{i}Mo" for i in [1,3,5,10,14,17,21]]
Xe_labels = [f"{i:.2f}Xe" for i in [0.25,0.50,0.75,1.00,1.25,1.50]]
MoXe_labels = [f"{i}Mo{j:.2f}Xe" for i in [1,3,5] for j in [0.25, 0.50, 0.75]]
conc_labels = {"None": none_labels, "Mo": Mo_labels, "Xe": Xe_labels, "MoXe": MoXe_labels}

# see plotly.com and dash.plotly.com for guidance on how to set this up.
app = Dash(__name__)
default_fig = go.Figure()
default_fig.update_layout(autosize = False, width = 1500, height = 550, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
ytickvals = [i*1e-10 for i in range(5,10)] + [i * 1e-9 for i in range(10)] + [i * 1e-8 for i in range(10)] + [1e-7]
yticktext = [str(i) if i in [1e-9, 1e-8, 1e-7] else "" for i in ytickvals]
default_fig.update_yaxes(title = r"Reduced Mobility (m<sup>2</sup>/s)", exponentformat = 'e', type = 'log', range = [-9.5,-7], tickvals = ytickvals, ticktext = yticktext)
default_fig.update_xaxes(showgrid = False, showticklabels = False, range = [0,1])

system_select_dict = [{'label': i, 'value': i} if not i == "None" else {'label': 'Pure U', 'value': i} for i in data.type.unique()]

app.layout = html.Div([
    html.H1('Mobility in the UMoXe System', style = {'margin': 'auto', 'width': '50%'}),
    html.H2("Select the data to plot"),
    html.Div([
        html.Div([html.Label('Select the system(s):'), dcc.Checklist(system_select_dict, id = 'impurity-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px', 'width': '6em'})], style = {'width': "280px", "float": "left"}),
        html.Div([html.Label('Select the rotation axis/axes:'), dcc.Checklist(data.axis.unique(), id = 'axis-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px'})], style = {'width': '280px', 'float': 'left', 'padding-left': '15px'}),
        html.Div([html.Label('Select the misorientation(s):'), dcc.Checklist(data.theta.unique(), id = 'misorientation-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px'})], style = {'width': "280px", "float": "left", 'padding-left': '15px'}),
        html.Div([html.Label('Select the concentration(s):'), dcc.Checklist(conc_dropdown_dict, id = 'concentration-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px', 'width': '8em', 'border': '1px green'})], style = {'width': "450px", "float": "left", 'padding-left': '15px'}),
        html.Div([html.Label('Select the temperature(s):'), dcc.Checklist(data["T"].unique(), id = 'temperature-checklist', labelStyle = {'display': 'inline-block', 'padding-left': '10px'})], style = {'width': "280px", "float": "left", 'padding-left': '15px'}),
    ]),
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
    html.Div([html.Div(dcc.Graph(figure = default_fig, id = "mobility-graph"), style = {'width': '80%', 'float': 'left'}),
    html.Div([html.Br(), html.Br(), html.Br(), html.Label('X-axis:'),
              dcc.Dropdown(["Arrhenius", "Concentration", "Misorientation"],
                            value = "Arrhenius", clearable = False, id = "plot-type-dropdown")],
              style = {'width': "200px", "display": "block", "float": "left", "padding-left": "50px"})]),
    html.Div([html.Label('Color data by:'),dcc.RadioItems([
        {"label": "Axis", "value": "axis"},
        {"label": "Misorientation", "value": "theta"},
        {"label": "Concentration", "value": "c"},
        {"label": "Temperature", "value": "T"}], "axis", inline = True, id = 'radio-group-by', inputStyle = {'margin-left': '20px'})], style = {'float': 'left', 'padding-left': '100px'}, hidden = True)

])

@app.callback(
    Output('concentration-checklist', 'options'),
    Output('concentration-checklist', 'value'),
    Output('misorientation-checklist', 'options'),
    Output('misorientation-checklist', 'value'),
    Output('temperature-checklist', 'options'),
    Output('temperature-checklist', 'value'),
    Output('mobility-graph', 'figure'),
    Input('impurity-checklist', 'value'),
    Input('axis-checklist', 'value'),
    Input('misorientation-checklist', 'value'),
    Input('concentration-checklist', 'value'),
    Input('temperature-checklist', 'value'),
    Input('plot-type-dropdown', 'value')
)
def update_app(impurity, axis, misorientation, concentration, temperature, plot_type):
    concentration_options = deepcopy(conc_dropdown_dict) # copies the list into a separate entity instead of making another reference
    misorientation_options = [{"label": i, "value": i} for i in data.theta.unique()]
    temperature_options = [{"label": i, "value": i} for i in data["T"].unique()]
    temperature_results = [] if not temperature else list(temperature)
    concentration_results = []
    misorientation_results = [] if not misorientation else list(misorientation)
    fig = go.Figure()
    fig.update_layout(autosize = False, width = 1500, height = 550, margin = dict(l = 5, r = 10, b = 50, t = 50, pad = 4))
    fig.update_yaxes(title = r"Reduced Mobility (m<sup>2</sup>/s)", exponentformat = 'e', type = 'log', range = [-9.5,-7], tickvals = ytickvals, ticktext = yticktext)
    fig.update_xaxes(showgrid = False, showticklabels = False, range = [0,1])
    if not impurity:
        for i,ii in [(j, item) for j, item in enumerate(concentration_options)]:
            concentration_options[i]["disabled"] = True
    else:
        for imp in ["None", "Mo", "Xe", "MoXe"]:
            for i,ii in [(j,item) for j, item in enumerate(concentration_options) if item["label"] in conc_labels[imp]]:
                if imp in impurity:
                    concentration_options[i]["disabled"] = False
                    if concentration and ii["value"] in concentration:
                        concentration_results.append(ii["value"])
                else:
                    concentration_options[i]["disabled"] = True

    for idx, _ in enumerate(concentration_options):
        if concentration_options[idx]["disabled"] == True:
            concentration_options[idx]["label"] = html.Div([conc_dropdown_dict[idx]["label"]], style = {'color': 'lightgray','display': 'inline-block'})
        else:
            concentration_options[idx]["label"] = conc_dropdown_dict[idx]["label"]


    selection = [True] * len(data)
    if plot_type == "Arrhenius":
        temperature_results = []
        for idx, _ in enumerate(misorientation_options):
            misorientation_options[idx]["disabled"] = False
        for idx, _ in enumerate(temperature_options):
            temperature_options[idx]["disabled"] = True
            temperature_options[idx]["label"] = html.Div([data["T"].unique()[idx]], style = {'color': 'lightgray', 'display': 'inline-block'})
        if not impurity or not axis or not misorientation or not concentration:
            return concentration_options, concentration_results, misorientation_options, misorientation_results, temperature_options, temperature_results, fig
        selection = [all(i) for i in zip(selection,update_selection(data, impurity, "type"))]
        selection = [all(i) for i in zip(selection,update_selection(data, axis, "axis"))]
        selection = [all(i) for i in zip(selection,update_selection(data, misorientation, "theta"))]
        selection = [all(i) for i in zip(selection,update_selection(data, concentration, "c_label"))]
        selection = [all(i) for i in zip(selection,update_selection(data, data["T"].unique(), "T"))]

        selected_data, key_names = data_as_dict(data[selection].set_index(['type', 'axis', 'c', 'theta', 'T']))
        for key in selected_data.keys():
            selected_data[key]["Arrhenius_T"] = convertToArrhenius(selected_data[key]["T"])
            selected_data[key]["Label"] = generateLabel(key, plot_type)
            fig.add_trace(go.Scatter(x = selected_data[key]["Arrhenius_T"], y = selected_data[key]["M"],
                            error_y = dict(type = 'data',  array = selected_data[key]["StdErr"], visible = True),
                            mode = 'markers', name = selected_data[key]["Label"]
                            # legendgroup = key[key_names.index(color_by)]
                            ))
        fig.update_xaxes(title = r"1/k<sub>B</sub>T (eV<sup>-1</sup>)", range = [8,12], showticklabels = True)
    elif plot_type == "Concentration":
        concentration_results = []
        for idx, _ in enumerate(misorientation_options):
            misorientation_options[idx]["disabled"] = False
        for idx, _ in enumerate(temperature_options):
            temperature_options[idx]["disabled"] = False
        for idx, _ in enumerate(concentration_options):
            concentration_options[idx]["disabled"] = True
            concentration_options[idx]["label"] = html.Div([conc_dropdown_dict[idx]["label"]], style = {'color': 'lightgray', 'display': 'inline-block'})
        if not impurity or not axis or not misorientation or not temperature:
            return concentration_options, concentration_results, misorientation_options, misorientation_results, temperature_options, temperature_results, fig
        selection = [all(i) for i in zip(selection,update_selection(data, impurity, "type"))]
        selection = [all(i) for i in zip(selection,update_selection(data, axis, "axis"))]
        selection = [all(i) for i in zip(selection,update_selection(data, misorientation, "theta"))]
        selection = [all(i) for i in zip(selection,update_selection(data, data["c_label"].unique(), "c_label"))]
        selection = [all(i) for i in zip(selection,update_selection(data, temperature, "T"))]

        selected_data, key_names = data_as_dict(data[selection].set_index(['type', 'axis', 'theta', 'T', 'c']))
        for key in selected_data.keys():
            selected_data[key]["Label"] = generateLabel(key, plot_type)
            fig.add_trace(go.Scatter(x = selected_data[key]["c"], y = selected_data[key]["M"],
                            error_y = dict(type = 'data',  array = selected_data[key]["StdErr"], visible = True),
                            mode = 'markers', name = selected_data[key]["Label"],
                            # legendgroup = key[key_names.index(color_by)]
                            ))
        xtickvals = []
        xticktext = []
        if len(impurity) == 1 or (len(impurity) == 2 and "None" in impurity):
            if "Mo" in impurity:
                r = [-0.005, 0.215]
                xtickvals = [0.01, 0.03, 0.05, 0.10, 0.14, 0.17, 0.21]
                xticktext = [f"{i*100:.0f}Mo" for i in xtickvals]
            elif "MoXe" in impurity:
                r = [-0.005, 0.065]
                xtickvals = [0.0125, 0.0150, 0.0175, 0.0325, 0.0350, 0.0375, 0.0525, 0.0550, 0.0575]
                xticktext = [f"{int(i*100):.0f}Mo{i-int(i*100):.2f}Xe" for i in xtickvals]
            else: # "Xe" in impurity
                r = [-0.005, 0.025]
                xtickvals = [i * 0.0025 for i in range(1,7)]
                xticktext = [f"{i*100:.2f}Xe" for i in xtickvals]
            if "None" in impurity:
                if len(impurity) == 1:
                    r = [-0.005, 0.025]
                xtickvals = [0] + xtickvals
                xticktext = ["Pure U"] + xticktext
        elif len(impurity) == 2 or (len(impurity) == 3 and "None" in impurity):
            if "Mo" in impurity and "Xe" in impurity:
                r = [-0.005, 0.215]
                xtickvals = [0.0025, 0.0050, 0.0075, 0.01, 0.0125, 0.0150, 0.03, 0.05, 0.10, 0.14, 0.17, 0.21]
                xticktext = ["0.25Xe", "0.50Xe", "0.75Xe", "1Mo\n1.00Xe", "1.25Xe", "1.50Xe", "3Mo", "5Mo", "10Mo", "14Mo", "17Mo", "21Mo"]
            elif "Mo" in impurity and "MoXe" in impurity:
                r = [-0.005, 0.215]
                xtickvals = [0.01, 0.0125, 0.0150, 0.0175, 0.03, 0.0325, 0.0350, 0.0375, 0.0525, 0.0550, 0.0575, 0.10, 0.14, 0.17, 0.21]
                xticktext = ["1Mo", "1Mo0.25Xe", "1Mo0.50Xe", "1Mo0.75Xe", "3Mo", "3Mo0.25Xe", "3Mo0.50Xe", "3Mo0.75Xe", "5Mo", "5Mo0.25Xe", "5Mo0.50Xe", "5Mo0.75Xe", "10Mo", "14Mo", "17Mo", "21Mo"]
            else: # "Xe" in impurity and "MoXe" in impurity:
                r = [-0.005, 0.065]
                xtickvals = [0.0025, 0.0050, 0.0075, 0.01, 0.0125, 0.0150, 0.0175, 0.03, 0.0325, 0.0350, 0.0375, 0.0525, 0.0550, 0.0575, 0.10, 0.14, 0.17, 0.21]
                xticktext = ["0.25Xe", "0.50Xe", "0.75Xe", "1.00Xe", "1Mo0.25Xe\n1.25Xe", "1Mo0.50Xe\n1.50Xe", "1Mo0.75Xe", "3Mo0.25Xe", "3Mo0.50Xe", "3Mo0.75Xe", "5Mo0.25Xe", "5Mo0.50Xe", "5Mo0.75Xe"]

            if "None" in impurity:
                xtickvals = [0] + xtickvals
                xticktext = ["Pure U"] + xticktext
        else: # len(impurity) == 4
            r = [-0.005, 0.215]
            xtickvals = [0, 0.0025, 0.0050, 0.0075, 0.01, 0.0125, 0.0150, 0.0175, 0.03, 0.0325, 0.0350, 0.0375, 0.05, 0.0525, 0.0550, 0.0575, 0.10, 0.14, 0.17, 0.21]
            xticktext = ["Pure U", "0.25Xe", "0.50Xe", "0.75Xe", "1Mo\n1.00Xe", "1Mo0.25Xe\n1.25Xe", "1Mo0.50Xe\n1.50Xe", "1Mo0.75Xe", "3Mo", "3Mo0.25Xe", "3Mo0.50Xe", "3Mo0.75Xe", "5Mo", "5Mo0.25Xe", "5Mo0.50Xe", "5Mo0.75Xe", "10Mo", "14Mo", "17Mo", "21Mo"]

        fig.update_xaxes(title = r"Concentration", range = r, showticklabels = True, tickvals = xtickvals, ticktext = xticktext, showgrid = True)
    elif plot_type == "Misorientation":
        misorientation_results = []
        for idx, _ in enumerate(misorientation_options):
            misorientation_options[idx]["disabled"] = True
            misorientation_options[idx]["label"] = html.Div([data.theta.unique()[idx]], style = {'color': 'lightgray', 'display': 'inline-block'})
        for idx, _ in enumerate(temperature_options):
            temperature_options[idx]["disabled"] = False
        if not impurity or not axis or not temperature or not concentration:
            return concentration_options, concentration_results, misorientation_options, misorientation_results, temperature_options, temperature_results, fig
        selection = [all(i) for i in zip(selection,update_selection(data, impurity, "type"))]
        selection = [all(i) for i in zip(selection,update_selection(data, axis, "axis"))]
        selection = [all(i) for i in zip(selection,update_selection(data, data["theta"].unique(), "theta"))]
        selection = [all(i) for i in zip(selection,update_selection(data, concentration, "c_label"))]
        selection = [all(i) for i in zip(selection,update_selection(data, temperature, "T"))]

        selected_data, key_names = data_as_dict(data[selection].set_index(['type', 'axis', 'T', 'c', 'theta']))
        for key in selected_data.keys():
            selected_data[key]["Label"] = generateLabel(key, plot_type)
            fig.add_trace(go.Scatter(x = selected_data[key]["theta"], y = selected_data[key]["M"],
                            error_y = dict(type = 'data',  array = selected_data[key]["StdErr"], visible = True),
                            mode = 'markers', name = selected_data[key]["Label"],
                            # legendgroup = key[key_names.index(color_by)]
                            ))
        xtickvals = [20,30,45]
        xticktext = [str(i) for i in xtickvals]
        fig.update_xaxes(title = r"Misorientation (degrees)", showticklabels = True, range = [10,55], tickvals = xtickvals, ticktext = xticktext)
    fig.update_yaxes(title = r"Reduced Mobility (m<sup>2</sup>/s)", exponentformat = 'e', type = 'log', range = [-9.5,-7])
    return concentration_options, concentration_results, misorientation_options, misorientation_results, temperature_options, temperature_results, fig

# @app.callback(
#     Output('radio-group-by', 'options'),
#     Output('radio-group-by', 'value'),
#     Input('plot-type-dropdown', 'value'),
#     Input('radio-group-by', 'value')
# )
def update_color_by_radio(plot_type, color_by):
    opts = [
        {"label": "Axis", "value": "axis"},
        {"label": "Misorientation", "value": "theta"},
        {"label": "Concentration", "value": "c"},
        {"label": "Temperature", "value": "T"},
    ]
    val = color_by
    if plot_type == "Arrhenius":
        opts[3]['disabled'] = True
        if color_by == opts[3]["value"]:
            val = opts[0]["value"]
    elif plot_type == "Concentration":
        opts[2]['disabled'] = True
        if color_by == opts[2]["value"]:
            val = opts[0]["value"]
    elif plot_type == "Misorientation":
        opts[1]['disabled'] = True
        if color_by == opts[1]["value"]:
            val = opts[0]["value"]
    return opts, val

if __name__ == '__main__':
    app.run_server(debug = True)
