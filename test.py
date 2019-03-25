#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import pandas as pd
import numpy as np

from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.tools import FigureFactory as FF
from ipwidgets import widgets
from IPython.display import display, clear_output, Image
from plotly.widgets import GraphWidget
plotly.offline.init_notebook_mode()

#import dash
#from dash.dependencies import Input, Output, Event
#import dash_core_components as dcc
#import dash_html_components as html
#import itertools

# Get data from db
df = pd.read_json("metoncofit.json", orient='columns')

# create the data structures for the three heatmaps
up = [go.Heatmap(z=df['values'].to_list(),
    colorscale='RdBu',
    x=df.index,
    y=df.columns)]
neut = [go.Heatmap(z=df['values'].to_list(),
    colorscale='RdBu',
    x=df.index,
    y=df.columns)]
down = [go.Heatmap(z=df['values'].to_list(),
    colorscale='RdBu',
    x=df.index,
    y=df.columns)]

g = go.FigureWidget(data=[up, neut, down])

# Assign widget values that will be used in app
#gene_number=widgets.FloatSlider(
#    value=5,
#    min=1,
#    max=100,
#    step=1,
#    description='Number of genes to display'
#    continuous_update=True
#)
tissue = widgets.Dropdown(
    description="Cancer Tissue:",
    value = 'Pan Cancer'
)
target = widgets.Dropdown(
    description="MetOncoFit Target:",
    value = "Differential Expression"
)
#container1 = widgets.HBox(children=[gene_number])
container2 = widgets.HBox(children=[tissue, target])

# Create the sliders
#steps = []
#for i in range(len(data)):
#    step = dict(
#        method='restyle'
#        args=['visible', [False]*len(data)],
#    )
#    step['args'][1][i] = True
#    steps.append(step)

#sliders=[dict(
#    active=5,
#    currentvalue={"prefix":"Number of genes: "},
#    pad={"t":50}
#    steps=steps
#)]

# write functions that will handle the inputs from the widgets:
def response(tissue, target):
    filter_list = [i and j for i, j in zip(df['Cancer']== tissue.value, df['Target']==target.value)]
    tmp = df[filtered_list]

    # Turn to 3 dfs
    up = tmp.loc[(tmp["type"] == "GAIN") | (tmp["type"] == "UPREG"))]
    neut = tmp.loc[(tmp["type"] == "NEUT") | (tmp["type"] == "NEUTRAL"))]
    down = tmp.loc[(tmp["type"] == "LOSS") | (tmp["type"] == "DOWNREG"))]

tissue.observe(response, names="value")
target.observe(response, names="value")

fig = tools.make_subplots(rows=1, cols=3)
fig.append_trace(up, 1, 1)
fig.append_trace(neut,1, 2)
fig.append_trace(down,1,3)

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
