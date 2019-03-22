#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_html - Make the html file containing all MetOncoFit data.

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
from bokeh.models.callbacks import CustomJS

# HTML file that will be outputted
output_file('test.html')

# Get 3 dataframes that will be turned into heatmaps
df = pd.read_json("metoncofit.json", orient='columns')
up = df.loc[(df["type"] == "UPREG") | (df["type"] == "GAIN")]
neut = df.loc[(df["type"] == "NEUTRAL") | (df["type"] == "NEUT")]
down = df.loc[(df["type"] == "DOWNREG") | (df["type"] == "LOSS")]

# Figure toolbar functions for interactions
tools_in_figure = ["hover, save, pan, box_zoom, reset, wheel_zoom"]
TOOLTIPS = [('Feature', '@feature'),('Gene', '@Gene'),('Value', '@value')]

# Set up heat map figure spaces to be filled in
hm1 = figure(x_axis_location='above', plot_height=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
hm2 = figure(x_axis_location='above', plot_height=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
hm3 = figure(x_axis_location='above', plot_height=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)

# Make source dicts that will be referenced in the callback functions. I added the cancer and target strings that will be referenced for callback. Finally, the empty source dicts contain empty values that will be filled in at call back
source_up = ColumnDataSource(data=dict(Gene=up["Gene"], feature=up["feature"], value=up["value"], Cancer=up["Cancer"], Target=up["Target"]))
final_up = ColumnDataSource(data=dict(Gene=up["Gene"], feature=up["feature"], value=up["value"], Cancer=up["Cancer"], Target=up["Target"]))

source_neut = ColumnDataSource(data=dict(Gene=neut["Gene"], feature=neut["feature"], value=neut["value"], Cancer=neut["Cancer"], Target=neut["Target"]))
final_neut = ColumnDataSource(data=dict(Gene=neut["Gene"], feature=neut["feature"], value=neut["value"], Cancer=neut["Cancer"], Target=neut["Target"]))

source_down = ColumnDataSource(data=dict(Gene=down["Gene"], feature=down["feature"], value=down["value"], Cancer=down["Cancer"], Target=down["Target"]))
final_down = ColumnDataSource(data=dict(Gene=down["Gene"], feature=down["feature"], value=down["value"], Cancer=down["Cancer"], Target=down["Target"]))

# Color bar and params
colors = brewer["RdBu"][8]
mapper = LinearColorMapper(palette=colors, low=0, high=1)
color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7pt",  border_line_color=None, location=(0,0))

# The actual figure vessels
hm1.rect(x="Gene", y="feature", width=1, height=1, source=final_up, line_color=None, fill_color=transform('value', mapper))
hm2.rect(x="Gene", y="feature", width=1, height=1, source=final_neut, line_color=None, fill_color=transform('value', mapper))
hm3.rect(x="Gene", y="feature", width=1, height=1, source=final_down, line_color=None, fill_color=transform('value', mapper))

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
up_callback = CustomJS(args=dict(source=source_up, out=final_up), code="""
// initialize some params
var data = source.data;
var gene = data['Gene'];
var feature = data['feature']
var value = data['value']
var target = data['Target'];
var cancer = data['Cancer'];

// create maps
var canc_dict = {
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

var target_dict = {
    "Differential Expression":"TCGA",
    "Copy Number Variation":"CNV",
    "Cancer Patient Survival":"SURV"
};

//var gene_no = up_gene_select.value;
function returnCancer(key){
    return canc_dict[key]
}

var cancer_type = returnCancer(cancer_type.value)

function returnType(key){
    return target_dict[key]
}

var target_type = returnType(target_type.value)

// these values will be filled after filtering for the target
var new_gene = [];
var new_feat = [];
var new_value = [];
var new_targ = [];
var new_canc = [];

// first get the Target type selected
for(var i=0; i<gene.length; i++){
    if(target[i]==targ_type){
        new_gene.push(gene[i]);
        new_feat.push(feature[i])
        new_value.push(value[i])
        new_targ.push(target[i])
        new_canc.push(cancer[i])
    }
}

// create array new array that looks like the original and capture the target
var newData = data.slice()
newData['Gene']=new_gene
newData['feature']=new_feat
newData['value']=new_value
newData['Target']=new_targ
newData['Cancer']=new_canc

// these values will be filled after filtering for the cancer
var final_gene = [];
var final_feat = [];
var final_value = [];
var final_targ = [];
var final_canc = [];

// now get the Cancer type selected
for (var j = 0; j <= new_canc.get_length(); j++){
    if (new_canc][j] == canc_type){
        final_gene.push(new_gene[j]);
        final_feat.push(new_feat[j])
        final_value.push(new_value[j])
        final_targ.push(new_targ[j])
        final_canc.push(new_canc[j])
    }
}
// make array that is filtered from both target and gene query
var finalData = newData.slice()
finalData['Gene']=final_gene
finalData['feature']=final_feat
finalData['value']=final_value
finalData['Target']=final_targ
finalData['Cancer']=final_canc

// these values will be filled after filtering for the number of genes
var forreal_gene = [];
var forreal_feat = [];
var forreal_value = [];
var forreal_targ = [];
var forreal_canc = [];

// Finally, let's return a specific number of genes specified by user input
//for (var k = 0; k <= final_gene.get_length(); k++){
//    if (final_gene.get_length()[k] < gene_no){
//        forreal_gene.push(final_gene[k]);
//        forreal_feat.push(final_feat[k])
//        forreal_value.push(final_value[k])
//        forreal_targ.push(final_targ[k])
//        forreal_canc.push(final_canc[k])
//    }
//}
// get it in df
//var forrealData = finaData.slice()
//forrealData['Gene']=forreal_gene
//forrealData['feature']=forreal_feat
//forrealData['value']=forreal_value

// return df
var out = finalData.slice()
out.change.emit()
""")

neut_callback = CustomJS(args=dict(source=source_neut, out=final_neut), code="""
// initialize some params
var data = source.data;

var gene = data['Gene'];
var feature = data['feature']
var value = data['value']
var target = data['Target'];
var cancer = data['Cancer'];

// create maps
var canc_dict = {
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
var target_dict = {
    "Differential Expression":"TCGA",
    "Copy Number Variation":"CNV",
    "Cancer Patient Survival":"SURV"
};

//var gene_no = up_gene_select.value;
function returnCancer(key){
    return canc_dict[key]
}
var cancer_type = returnCancer(cancer_type.value)

function returnType(key){
    return target_dict[key]
}
var target_type = returnType(target_type.value)

// these values will be filled after filtering for the target
var new_gene = [];
var new_feat = [];
var new_value = [];
var new_targ = [];
var new_canc = [];

// these values will be filled after filtering for the cancer
var final_gene = [];
var final_feat = [];
var final_value = [];
var final_targ = [];
var final_canc = [];

// these values will be filled after filtering for the number of genes
var forreal_gene = [];
var forreal_feat = [];
var forreal_value = [];
var forreal_targ = [];
var forreal_canc = [];

// first get the Target type selected
for(var i=0; i<gene.length; i++){
    if(target[i]==targ_type){
        new_gene.push(gene[i]);
        new_feat.push(feature[i])
        new_value.push(value[i])
        new_targ.push(target[i])
        new_canc.push(cancer[i])
    }
}

// create array new array that looks like the original and capture the target
var newData = data.slice()
newData['Gene']=new_gene
newData['feature']=new_feat
newData['value']=new_value
newData['Target']=new_targ
newData['Cancer']=new_canc

// now get the Cancer type selected
for (var j = 0; j <= new_canc.get_length(); j++){
    if (new_canc][j] == canc_type){
        final_gene.push(new_gene[j]);
        final_feat.push(new_feat[j])
        final_value.push(new_value[j])
        final_targ.push(new_targ[j])
        final_canc.push(new_canc[j])
    }
}
// make array that is filtered from both target and gene query
var finalData = newData.slice()
finalData['Gene']=final_gene
finalData['feature']=final_feat
finalData['value']=final_value
finalData['Target']=final_targ
finalData['Cancer']=final_canc

// Finally, let's return a specific number of genes specified by user input
//for (var j = 0; j <= final_gene.get_length(); j++){
//    if (final_gene.get_length()[j] < gene_no){
//        forreal_gene.push(final_gene[j]);
//        forreal_feat.push(final_feat[j])
//        forreal_value.push(final_value[j])
//        forreal_targ.push(final_targ[j])
//        forreal_canc.push(final_canc[j])
//    }
//}
// get it in df
//var forrealData = finaData.slice()
//forrealData['Gene']=forreal_gene
//forrealData['feature']=forreal_feat
//forrealData['value']=forreal_value

// return df
var out = finalData.slice()
out.change.emit()
""")

down_callback = CustomJS(args=dict(source=source_down, out=final_down), code="""
// initialize some params
var data = source.data;
var gene = data['Gene'];
var feature = data['feature']
var value = data['value']
var target = data['Target'];
var cancer = data['Cancer'];

// create maps
var canc_dict = {
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
var target_dict = {
    "Differential Expression":"TCGA",
    "Copy Number Variation":"CNV",
    "Cancer Patient Survival":"SURV"
};

//var gene_no = up_gene_select.value;
function returnCancer(key){
    return canc_dict[key]
}
var cancer_type = returnCancer(cancer_type.value)

function returnType(key){
    return target_dict[key]
}
var target_type = returnType(target_type.value)

// these values will be filled after filtering for the target
var new_gene = [];
var new_feat = [];
var new_value = [];
var new_targ = [];
var new_canc = [];

// these values will be filled after filtering for the cancer
var final_gene = [];
var final_feat = [];
var final_value = [];
var final_targ = [];
var final_canc = [];

// these values will be filled after filtering for the number of genes
var forreal_gene = [];
var forreal_feat = [];
var forreal_value = [];
var forreal_targ = [];
var forreal_canc = [];

// first get the Target type selected
for(var i=0; i<gene.length; i++){
    if(target[i]==targ_type){
        new_gene.push(gene[i]);
        new_feat.push(feature[i])
        new_value.push(value[i])
        new_targ.push(target[i])
        new_canc.push(cancer[i])
    }
}

// create array new array that looks like the original and capture the target
var newData = data.slice()
newData['Gene']=new_gene
newData['feature']=new_feat
newData['value']=new_value
newData['Target']=new_targ
newData['Cancer']=new_canc

// now get the Cancer type selected
for (var j = 0; j <= new_canc.get_length(); j++){
    if (new_canc][j] == canc_type){
        final_gene.push(new_gene[j]);
        final_feat.push(new_feat[j])
        final_value.push(new_value[j])
        final_targ.push(new_targ[j])
        final_canc.push(new_canc[j])
    }
}
// make array that is filtered from both target and gene query
var finalData = newData.slice()
finalData['Gene']=final_gene
finalData['feature']=final_feat
finalData['value']=final_value
finalData['Target']=final_targ
finalData['Cancer']=final_canc

// Finally, let's return a specific number of genes specified by user input
//for (var j = 0; j <= final_gene.get_length(); j++){
//    if (final_gene.get_length()[j] < gene_no){
//        forreal_gene.push(final_gene[j]);
//        forreal_feat.push(final_feat[j])
//        forreal_value.push(final_value[j])
//        forreal_targ.push(final_targ[j])
//        forreal_canc.push(final_canc[j])
//    }
//}
// get it in df
//var forrealData = finaData.slice()
//forrealData['Gene']=forreal_gene
//forrealData['feature']=forreal_feat
//forrealData['value']=forreal_value
//forrealData['Target']=forreal_targ
//forrealData['Cancer']=forreal_canc

// return df
var out = finalData.slice()
out.change.emit()
""")

# Input controls

# Drop down menu to choose cancer and target
cancer_type = Select(title="Cancer type:", value="Pan Cancer", options=["Breast Cancer", "Glioma", "Colorectal Cancer", "Lung Cancer", "Melanoma", "Renal Cancer", "Prostate Cancer", "Ovarian Cancer", "B-cell Lymphoma", "Pan Cancer"], callback=neut_callback)
up_callback.args['cancer_type'] = cancer_type
down_callback.args['cancer_type'] = cancer_type
neut_callback.args['cancer_type'] = cancer_type

target_type = Select(title="MetOncoFit Predictions:", value="Differential Expression", options=["Differential Expression", "Copy Number Variation", "Cancer Patient Survival"], callback=neut_callback)
up_callback.args['target_type'] = target_type
down_callback.args['target_type'] = target_type
neut_callback.args['target_type'] = target_type

# Slider to choose number of genes to show
#if target_type.value == 'Copy Number Variation':
#    targ_labels = ["GAIN", "NEUT", "LOSS"]
#else:
#    targ_labels = ["UPREG", "NEUTRAL", "DOWNREG"]

#up_gene_select = Slider(start=0, end=100, value=5, step=5, title="Number of up-regulated genes to display", callback=up_callback)
#up_callback.args['up_gene_select'] = up_gene_select

#neut_gene_select = Slider(start=0, end=100, value=5, step=5, title="Number of neutrally-regulated genes to display", callback=neut_callback)
#neut_callback.args['neut_gene_select'] = neut_gene_select

#down_gene_select = Slider(start=0, end=100, value=5, step=5, title="Number of down-regulated genes to display", callback=down_callback)
#down_callback.args['down_gene_select'] = down_gene_select

# Text input box for selecting a custom list of genes from the user
#gene_list = TextInput(value="Gene symbols (ie: IDH2, TDO2, PDGDH, ...)", title="Enter a list of gene symbols to query")
#gene_list.json_on_change('value', callback)

# Additional plotting parameters
#plots = gridplot([[hm1, hm2, hm3]], sizing_mode='fixed')
widgets = widgetbox(cancer_type, target_type)#up_gene_select, neut_gene_select, down_gene_select,)
layout = row(widgets, hm1, hm2, hm3)
show(layout)
