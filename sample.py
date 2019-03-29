# -*- coding: utf-8 -*-
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph.objs as go

import pandas as pd
import json

# Import json
input_data = pd.read_json('metoncofit.json', orient='columns')
data = []

# Create map that will reference values in dataset and what will be shown in the actual website
cancer_dict = {
    "Breast Cancer":"Breast",
    "Glioma":"CNS",
    "Colorectal Cancer":"Colorectal",
    "Lung Cancer":"Lung",
    "Melanoma":"Melanoma",
    "Renal Cancer":"Renal",
    "Prostate Cancer":"Prostate",
    "Ovarian Cancer":"Ovarian",
    "B-cell Lymphoma":"Leukemia",
    "Pan Cancer":"Pan"
};

# This will be the options show in the dropdown menus for both cancers and targets
cancers = ["B-Cell Lymphoma", "Breast Cancer", "Colorectal Cancer", "Glioma", "Lung Cancer", "Melanoma", "Prostate Cancer", "Ovarian Cancer", "Renal Cancer" "Pan Cancer"]
targets = ["Differential Expression", "Copy Number Variation", "Patient Survival"]

# Set default values to show Pan Cancer and Differential Expression Values
cancer_default = cancers[9]
df_default = input_data[(input_data['Cancer'] == cancer_default) && (input_data['Target'] == target_default)]

app = dash.Dash()
app.layout = html.Div([
    # Create first dropdown
    html.Div([
        html.H4('Select Cancer Type:'),
        dcc.Dropdown(
            id='cancer_dropdown',
            options=[{'label':i,"value":i} for i in lst], value=df_default
        )
    ],
    style = {'width':'50%', 'display':'inline-block'}),

    # Create second dropdown
    html.Div([
        html.H4('Select MetOncoFit Target:')
        dcc.Dropdown(
            id='target_dropdown',
            value='default'
        ),
        ],
        style={'width':'50%', 'float':'right', 'display':'inline-block'}),

    # Heatmap variable
    dcc.Graph(
        id='heatmap',
        figure={
            'data': [go.Heatmap(
                x=df_default['features'],
                y=df_default['Gene'],
                z=df_default['value'],
                name='',
                colorscale='RdBu'
                )],
            'layout': go.Layout(
                xaxis=dict('title':'Genes'),
                yaxis=dict('title':'Features'),
                )})
    ]),
])

# Callback decorator
@app.callback(
    dash.dependencies.Output(component_id='cancer_dropdown', component_property='options')
    [dash.dependencies.Input(component_id='target_dropdown', component_property='value')]

    Output('output', 'children'),
    [Input('heatmap', 'hoverData'),
     Input('heatmap', 'clickData')])
     )

# Update target type:
def update_targets():
    return [{'label':i, 'value':i} for i in input_data[input_data['Cancer']==selected]['Item'].unique()]

# Update heatmap
def update_heatmap(cancer_dropdown, target_dropdown):
    hm_data = input_data[(input_data["Cancer"]==)]
    hm)data = pd.merge(data)

# Display the hover data
def display_hoverdata(hoverData, clickData):
    return[
        json.dump(hoverData, indent=2),
        html.Br(),
        json.dumps(clickData, indent=2)
        ]


if __name__ == '__main__':
    app.run_server(debug=True)
