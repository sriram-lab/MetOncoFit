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

# get the unique values that will be used in the widget
num_uniq_genes = df["Gene"].nunique()
cancer_type = df["Cancer"].unique()
prediction_type = df["Target"].unique()

# Start the application using a css file from someone
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
    style={'width': '49%', 'display': 'inline-block'}
    ),

    # Create the slider
    html.Div(dcc.Slider(
        id='gene-slider',
        min=1,
        max=df.index.max(),
        value=5,
        #marks={gene:gene for gene in df['Gene'].unique()}
    ),
    style={'width':'49%', 'padding':'0px 20px 20px 20px'}
    ),

    # Finally, create the heatmap vessel
    dcc.Graph(
        id='heatmap',
        figure={
            'data': [(
                go.Heatmap(
                    x=df['Gene'],
                    y=df['feature'],
                    z=df['value'],
                    name='heatmap-test',
                    colorscale='RdBu')
                    )],
            'layout':go.Layout(
                xaxis=dict(title='Genes'),
                yaxis=dict(title='Features')
            )
        }
    )
])

# Create callback decorator
@app.callback(
    dash.dependencies.Output('heatmap', 'figure'),
    [dash.dependencies.Input('cancer-type', 'value'),
    dash.dependencies.Input('prediction-type', 'value'),
    dash.dependencies.Input('gene-slider', 'value')])

def update_graph(cancer_choice, prediction_choice, slider_amount):
    dff = df[df['Target'] == prediction_choice]
    dff = dff[dff['Cancer'] == cancer_choice]
    dff = dff[dff.index < slider_amount]
    return {
        'data': [(
            go.Heatmap(
                x=dff['Gene'],
                y=dff['feature'],
                z=dff['value'],
                name='heatmap-test',
                colorscale='RdBu')
                )],
        'layout':go.Layout(
            xaxis=dict(title='Genes'),
            yaxis=dict(title='Features')
        )
    }

if __name__ == '__main__':
    app.run_server(debug=True)

"""

# Slider widget that controls number of genes shown for all dataframes
gene_number=widgets.FloatSlider(
    value=5,
    min=1,
    max=num_uniq_genes,
    step=1,
    description='Number of genes to display',
    continuous_update=True
)
# Dropdown widget that will show several cancer type options
tissue = widgets.Dropdown(
    options=list(cancer_type),
    description="Cancer Tissue:",
    value = 'Pan'
)
# Dropdown widget that will show several target prediction options
target = widgets.Dropdown(
    options=list(prediction_type),
    description="MetOncoFit Target:",
    value = "TCGA"
)

# Create the empty figures that will have the traces and the responding widget
up = [go.Heatmap(z=df['value'].to_list(),
    colorscale='RdBu',
    x=df['Gene'],
    y=df['feature'])]
neut = [go.Heatmap(z=df['value'].to_list(),
    colorscale='RdBu',
    x=df['Gene'],
    y=df['feature'])]
down = [go.Heatmap(z=df['value'].to_list(),
    colorscale='RdBu',
    x=df['Gene'],
    y=df['feature'])]

#heatmap_up = up.iplot(kind='heatmap', colorscale='spectral', filename='test')
#heatmap_neut = neut.iplot(kind='heatmap', colorscale='spectral', filename='test')
#heatmap_down = down.iplot(kind='heatmap', colorscale='spectral', filename='test')

g = GraphWidget()

# write functions that will handle the inputs from the widget and alter the graphs:

def response(change):
    filter_list = [i and j for i, j in zip(df['Cancer'] == tissue.value, df['Target'] == target.value)]
    new_df = df[filtered_list]
    g.restyle({'x':[]})
    g.layout.xaxis.title='Genes'
    g.layout.yaxis.title='Features'

tissue.observe(response, names="value")
target.observe(response, names="value")

#fig = tools.make_subplots(rows=1, cols=3)
#fig.append_trace(up, 1, 1)
#fig.append_trace(neut,1, 2)
#fig.append_trace(down,1,3)

container1 = widgets.HBox(children=[gene_number])
container2 = widgets.HBox(children=[tissue, target])
all_widgets = widgets.VBox([container1, container2])

# Edit the hovertext
#hovertext = list()
#for fidx, feat in enumerate(features):
#    hovertext.append(list())
#    for gidx, gene in enumerate(genes):
#        hovertext[-1].append('Gene: {}</br />Feature: {}<br />Value: {}'.format(gene, feat, value[fidx][gidx]))

# Graph params
#layout=go.Layout(
#    autosize=True,
#    margin=dict(t=0, b=0, l=0, r=0)
#    )
#)
display(all_widgets)
display(g)

updatemenus=list([
    dict(
        buttons=list([
            dict(
                args=['Cancer Tissue', 'Breast']
                label='Breast Cancer'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'CNS']
                label='Glioma'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Colorectal']
                label='Colorectal Cancer'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Lung']
                label='Lung Cancer'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Melanoma']
                label='Melanoma'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Renal']
                label='Renal Cancer'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Prostate']
                label='Prostate Cancer'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Ovarian']
                label='Ovarian Cancer'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Leukemia']
                label='B-cell Lymphoma'
                method='restyle'
            ),
            dict(
                args=['Cancer Tissue', 'Pan']
                label='Pan Cancer'
                method='restyle'
            ),
        ]),
        direction='down',
        pad={'r':10,'t'=10},
        showactive=True,
        x=0.3,
        xanchor='left',
        y=1.12,
        yanchor='top'
    ),
    dict(
        buttons=list([
            dict(
                args=['MetOncoFit Target', 'TCGA']
                label='Differential Expression'
                method='restyle'
            ),
            dict(
                args=['MetOncoFit Target', 'CNV']
                label='Copy Number Variation'
                method='restyle'
            ),
            dict(
                args=['MetOncoFit Target', 'SURV']
                label='Cancer Patient Survival'
                method='restyle'
            ),
        ]),
        direction='down',
        pad={'r':10,'t'=10},
        showactive=True,
        x=0.3,
        xanchor='left',
        y=1.0,
        yanchor='top'
    )
])

annotations = list([
    dict(text='Cancer<br>Tissue', x=0, y=1.11, yref='paper', align='left', showarrow=False),
    dict(text='MetOncoFit<br>Target', x=0.25, y=1.11, yref='paper', align='left', showarrow=False)
])

layout['updatemenus'] = updatemenus
layout['annotations'] = annotations

fig['layout'] = layout
py.iplot(fig, filename='')
"""
