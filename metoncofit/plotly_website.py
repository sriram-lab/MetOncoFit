#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Dash website
@author: Scott Campit
"""
import pandas as pd
import numpy as np

from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.tools import FigureFactory as FF

from ipywidgets import widgets
from IPython.display import display, clear_output, Image
from plotly.widgets import GraphWidget
import cufflinks as cf

import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html

# Get data from MetOncoFit database
df = pd.read_json("metoncofit.json", orient='columns')
up = df.loc[(df["type"] == "UPREG") | (df["type"] == "GAIN")]
neut = df.loc[(df["type"] == "NEUTRAL") | (df["type"] == "NEUT")]
down = df.loc[(df["type"] == "DOWNREG") | (df["type"] == "LOSS")]

# get values that will be used in the widget
num_uniq_genes = df["Gene"].nunique()
cancer_type = df["Cancer"].unique()
prediction_type = df["Target"].unique()

# Start the application using a css file
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    # Create the drop down menus
    html.Div([
        dcc.Dropdown(
            id='cancer-type',
            options=[{'label':i, 'value':i} for i in cancer_type],
            value='Pan'
        ),
        dcc.Dropdown(
            id='prediction-type',
            options=[{'label':i,'value':i} for i in prediction_type],
            value='TCGA'
        )
    ],
    style={'width': '10%', 'display': 'inline-block'}
    ),

    # Create the slider
    html.Div(dcc.Slider(
        id='gene-slider',
        min=1,
        max=num_uniq_genes,
        step=1,
        value=5,
        marks={
            5: {'label':'5'},
            25: {'label':'25'},
            100: {'label':'100'},
            250: {'label':'250'},
            500: {'label':'500'},
            750: {'label':'750'},
            num_uniq_genes: {'label':str(num_uniq_genes)}
        },
        updatemode='drag'
    ),
    style={'width':'49%', 'padding':'0px 20px 20px 20px'},
    ),
    html.Div(id='updatemode-output-container', style={'margin-top':20}),

    # Finally, create the heatmap vessels
    html.Div(dcc.Graph(
        id='up-heatmap',
        figure={
            'data': [(
                go.Heatmap(
                    x=up['Gene'],
                    y=up['feature'],
                    z=up['value'],
                    name='up',
                    colorscale='RdBu')
                    )],
            'layout':go.Layout(
                title='Increased',
                xaxis=dict(title='Genes'),
                yaxis=dict(title='Features', automargin=True)
            )
        }
    ),
    style={'width':'75%', 'padding':'0px 20px 20px 20px'},
    ),
    html.Div(dcc.Graph(
        id='neut-heatmap',
        figure={
            'data': [(
                go.Heatmap(
                    x=neut['Gene'],
                    y=neut['feature'],
                    z=neut['value'],
                    name='neut',
                    colorscale='RdBu')
                    )],
            'layout':go.Layout(
                xaxis=dict(title='Genes'),
                yaxis=dict(title='Features', automargin=True)
            )
        }
    ),
    style={'width':'75%', 'padding':'0px 20px 20px 20px'},
    ),
    html.Div(dcc.Graph(
        id='down-heatmap',
        figure={
            'data': [(
                go.Heatmap(
                    x=down['Gene'],
                    y=down['feature'],
                    z=down['value'],
                    name='down',
                    colorscale='RdBu')
                    )],
            'layout':go.Layout(
                xaxis=dict(title='Genes'),
                yaxis=dict(title='Features', automargin=True)
            )
        }
    ),
    style={'width':'75%', 'padding':'0px 20px 20px 20px'},
    )
])

# Create callback decorator
@app.callback(
    dash.dependencies.Output('up-heatmap', 'figure'),
    [dash.dependencies.Input('cancer-type', 'value'),
    dash.dependencies.Input('prediction-type', 'value'),
    dash.dependencies.Input('gene-slider', 'value')])

def update_up(cancer_choice, prediction_choice, slider_choice):
    up_df = up[up['Target'] == prediction_choice]
    up_df = up_df[up_df['Cancer'] == cancer_choice]
    up_tmp = up_df['Gene'].unique().tolist()
    up_tmp = up_tmp[0:slider_choice]
    up_df = up_df.loc[up_df['Gene'].isin(up_tmp)]

    return {
        'data': [(
            go.Heatmap(
                x=up_df['Gene'],
                y=up_df['feature'],
                z=up_df['value'],
                name='up-heatmap',
                colorscale='RdBu')
                )],
        'layout':go.Layout(
            title="Increased",
            xaxis=dict(title='Genes'),
            yaxis=dict(title='Features')
        )
    }

@app.callback(
    dash.dependencies.Output('neut-heatmap', 'figure'),
    [dash.dependencies.Input('cancer-type', 'value'),
    dash.dependencies.Input('prediction-type', 'value'),
    dash.dependencies.Input('gene-slider', 'value')])

def update_neut(cancer_choice, prediction_choice, slider_choice):
    neut_df = neut[neut['Target'] == prediction_choice]
    neut_df = neut_df[neut_df['Cancer'] == cancer_choice]
    neut_tmp = neut_df['Gene'].unique().tolist()
    neut_tmp = neut_tmp[0:slider_choice]
    neut_df = neut_df.loc[neut_df['Gene'].isin(neut_tmp)]

    return {
        'data': [(
            go.Heatmap(
                x=neut_df['Gene'],
                y=neut_df['feature'],
                z=neut_df['value'],
                name='neut-heatmap',
                colorscale='RdBu')
                )],
        'layout':go.Layout(
            title="Neutral",
            xaxis=dict(title='Genes'),
            yaxis=dict(title='Features')
        )
    }

@app.callback(
    dash.dependencies.Output('down-heatmap', 'figure'),
    [dash.dependencies.Input('cancer-type', 'value'),
    dash.dependencies.Input('prediction-type', 'value'),
    dash.dependencies.Input('gene-slider', 'value')])

def update_down(cancer_choice, prediction_choice, slider_choice):
    down_df = down[down['Target'] == prediction_choice]
    down_df = down_df[down_df['Cancer'] == cancer_choice]
    down_tmp = down_df['Gene'].unique().tolist()
    down_tmp = down_tmp[0:slider_choice]
    down_df = down_df.loc[down_df['Gene'].isin(down_tmp)]

    return {
        'data': [(
            go.Heatmap(
                x=down_df['Gene'],
                y=down_df['feature'],
                z=down_df['value'],
                name='down-heatmap',
                colorscale='RdBu')
                )],
        'layout':go.Layout(
            title="Decreased",
            xaxis=dict(title='Genes'),
            yaxis=dict(title='Features')
        )
    }

@app.callback(Output('updatemode-output-container', 'children'),
    [dash.dependencies.Input('gene-slider', 'value')])

def display_value(value):
    return 'Maximum number of genes displayed: {}'.format(value, value)

if __name__ == '__main__':
    app.run_server(debug=True)
