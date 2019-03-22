#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""


@author: Scott Campit
"""

import numpy as np
import pandas as pd
from math import pi

from bokeh.io import output_file, show, curdoc
from bokeh.models import BasicTicker, ColorBar, ColumnDataSource, LinearColorMapper, PrintfTickFormatter, CustomJS, BoxSelectTool
from bokeh.plotting import figure
from bokeh.transform import transform
from bokeh.layouts import column, row, widgetbox, gridplot
from bokeh.models.widgets import Button, RadioButtonGroup, Select, Slider, TextInput
from bokeh.palettes import brewer
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.models.callbacks import CustomJS

output_file('test.html')

df = pd.read_json("metoncofit.json", orient='columns')
up = df.loc[(df["type"] == "UPREG" | df["type"] == "GAIN")]
neut = df.loc[(df["type"] == "NEUTRAL" | df["type"] == "NEUT")]
down = df.loc[(df["type"] == "DOWNREG" | df["type"] == "LOSS")]

# Figure toolbar functions for interactions
tools_in_figure = ["hover, save, pan, box_zoom, reset, wheel_zoom"]
TOOLTIPS = [('Feature', '@feature'),('Gene', '@Gene'),('Value', '@value')]

# Set up heat map figure spaces to be filled in
hm1 = figure(x_axis_location='above', plot_height=400, plot_width=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
hm2 = figure(x_axis_location='above', plot_height=400, plot_width=10000, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
hm3 = figure(x_axis_location='above', plot_height=400, plot_width=3000, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)

source_up = columnDataSource(data=dict(Gene=up["Gene"], feature=up["feature"], value=up["value"]))
source_neut = columnDataSource(data=dict(Gene=neut["Gene"], feature=neut["feature"], value=neut["value"]))
source_down = columnDataSource(data=dict(Gene=down["Gene"], feature=down["feature"], value=down["value"]))

# The actual figure vessels
hm1.rect(x="Gene", y="feature", width=1, height=1, source=source_up, line_color=None, fill_color=transform('value', mapper))
hm2.rect(x="Gene", y="feature", width=1, height=1, source=source_neut, line_color=None, fill_color=transform('value', mapper))
hm3.rect(x="Gene", y="feature", width=1, height=1, source=source_down, line_color=None, fill_color=transform('value', mapper))

# Color bar and params
colors = brewer["RdBu"][8]
mapper = LinearColorMapper(palette=colors, low=0, high=1)
color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7pt",  border_line_color=None, location=(0,0))

# Figure params
hm1.grid.grid_line_color=None
hm2.grid.grid_line_color=None
hm3.grid.grid_line_color=None
hm1.axis.axis_line_color=None
hm2.axis.axis_line_color=None
hm3.axis.axis_line_color=None
hm1.axis.major_tick_line_color=None
hm2.axis.major_tick_line_color=None
hm3.axis.major_tick_line_color=None
hm1.axis.major_label_text_font_size='7pt'
hm2.axis.major_label_text_font_size='7pt'
hm3.axis.major_label_text_font_size='7pt'
hm1.yaxis.visible=True
hm2.yaxis.visible=False
hm3.yaxis.visible=False
hm1.axis.major_label_standoff=0
hm2.axis.major_label_standoff=0
hm3.axis.major_label_standoff=0
hm1.xaxis.major_label_orientation = pi/3
hm2.xaxis.major_label_orientation = pi/3
hm3.xaxis.major_label_orientation = pi/3

hm1.add_layout(color_bar, 'left')

# Custom subfunctions that will only be used in the scope of heatmap_html
def selection_changes():
    """
    This module will update the number of columns shown in the plot for a single heatmap based on user input via the button widget. If there are genes in the gene list field, it will take a list of genes from user input to show where on the heatmap they are. If there is no match, it needs to print somewhere, and just not show the gene that isn't in the dataset.
    """

    cancer_map = {
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
    }

    target_map = {
    "Differential Expression":"TCGA",
    "Copy Number Variation":"CNV",
    "Cancer Patient Survival":"SURV"
    }

    # Number of genes selected (Default = 5)
    up_gene_val = up_gene_select.value
    neut_gene_val = neut_gene_select.value
    down_gene_val = down_gene_select.value

    # The cancer type and target prediction
    canc_select = cancer_type.value.strip()
    target_select = target_type.value.strip()
    # First get the values corresponding to the cancer and target the user wants to query
    selected = df.loc[(df['Target'] == target_map.get(target_select)) & (df['Cancer'] == cancer_map.get(canc_select))]

    # Second, separate into 3 dataframes that will be used to make the 3 heat maps
    up_select = selected.loc[(selected["type"] == targ_labels[0])]
    neut_select = selected.loc[(selected["type"] == targ_labels[1])]
    down_select = selected.loc[(selected["type"] == targ_labels[2])]

    # Third, get the number of genes they want to display.
    up_select = up_select.sample(n=up_gene_val)
    neut_select = neut_select.sample(n=neut_gene_val)
    down_select = down_select.sample(n=down_gene_val)

    # Fourth, if they specified some genes:
    #if gene_list != ("Gene symbols (ie: IDH2, TDO2, PDGDH, ...)" or ""):
    #    up_select = up_select.loc[up_select["Gene"].str.contains(gene_list)==True]
    #    neut_select = neut_select.loc[neut_select["Gene"].str.contains(gene_list)==True]
    #    down_select = down_select.loc[down_select["Gene"].str.contains(gene_list)==True]

    return up_select, neut_select, down_select

callback = CustomeJS(args=dict(up=up, neut=neut, down=down),
code="""
    var up_data =

""")

# Input controls:
# Slider to choose number of genes to show
up_gene_select = Slider(start=0, end=len(mask["type"] == targ_labels[0]), value=5, step=5, title="Number of upregulated genes to display")
up_gene_select.js_on_change('value', callback)

neut_gene_select = Slider(start=0, end=len(mask["type"] == targ_labels[1]), value=5, step=1, title="Number of unregulated genes to display")
neut_gene_select.js_on_change('value', callback)

down_gene_select = Slider(start=0, end=len(mask["type"] == targ_labels[2]), value=5, step=1, title="Number of downregulated genes to display")
down_gene_select.js_on_change('value', callback)

# Drop down menu to choose cancer and target
cancer_type = Select(title="Cancer type:", value="Pan Cancer", options=["Breast Cancer", "Glioma", "Colorectal Cancer", "Lung Cancer", "Melanoma", "Renal Cancer", "Prostate Cancer", "Ovarian Cancer", "B-cell Lymphoma", "Pan Cancer"])
cancer_type.json_on_change('value', callback)

target_type = Select(title="MetOncoFit Predictions:", value="Differential Expression", options=["Differential Expression", "Copy Number Variation", "Cancer Patient Survival"])
target_type.json_on_change('value', callback)

# Text input box for selecting a custom list of genes from the user
gene_list = TextInput(value="Gene symbols (ie: IDH2, TDO2, PDGDH, ...)", title="Enter a list of gene symbols to query")
gene_list.json_on_change('value', callback)

# Additional plotting parameters
plots = gridplot([[hm1, hm2, hm3]], sizing_mode='fixed')
layout = column()
show(plots)
