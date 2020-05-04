#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_html - Make the html file containing all MetOncoFit data.

@author: Scott Campit
"""
import copy
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
#up = df.loc[(df["type"] == "UPREG") | (df["type"] == "GAIN")]
#neut = df.loc[(df["type"] == "NEUTRAL") | (df["type"] == "NEUT")]
#down = df.loc[(df["type"] == "DOWNREG") | (df["type"] == "LOSS")]

# Color bar and params
colors = brewer["RdBu"][8]
mapper = LinearColorMapper(palette=colors, low=0, high=1)
color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7pt",  border_line_color=None, location=(0,0))

# Figure toolbar functions for interactions
tools_in_figure = ["hover, save, pan, box_zoom, reset, wheel_zoom"]
TOOLTIPS = [('Feature', '@feature'),('Gene', '@Gene'),('Value', '@value')]

# Set up heat map figure spaces to be filled in
hm1 = figure(x_axis_location='above', plot_height=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
hm2 = figure(x_axis_location='above', plot_height=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
hm3 = figure(x_axis_location='above', plot_height=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)

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

# Make source dicts that will be referenced in the callback functions. I added the cancer and target strings that will be referenced for callback. Finally, the empty source dicts contain empty values that will be filled in at call back
data=dict(Gene=df["Gene"], feature=df["feature"], value=df["value"], Cancer=df["Cancer"], Target=df["Target"])
source = ColumnDataSource(copy.deepcopy(data))
up = ColumnDataSource(copy.deepcopy(data))
neut = ColumnDataSource(copy.deepcopy(data))
down = ColumnDataSource(copy.deepcopy(data))

# The actual figure vessels corresponding to 3 heatmaps
hm1.rect(x="Gene", y="feature", width=1, height=1, source=up, line_color=None, fill_color=transform('value', mapper))
hm2.rect(x="Gene", y="feature", width=1, height=1, source=neut, line_color=None, fill_color=transform('value', mapper))
hm3.rect(x="Gene", y="feature", width=1, height=1, source=down, line_color=None, fill_color=transform('value', mapper))

# the callback function to make this thing responsive
callback = CustomJS(args=dict(source=source, up=up, neut=neut, down=down), code="""
    // create cancer and target dictionaries to map names back to values
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

    // some functions I need for mapping
    //var gene_no = up_gene_select.value;
    function returnCancer(key){
        return canc_dict[key]
    }
    var cancer_type = returnCancer(cancer_type.value)

    function returnType(key){
        return target_dict[key]
    }
    var target_type = returnType(target_type.value)
    console.log(target_type)

    // Get Target and Cancer specifications
    function filter_by_target(in, target, cancer){
        for(var i=0; i<in.get_length(); i++){
            if(in.data['Target'][i] == target_type && in.data['Cancer'][j] == cancer_type){
                out.data['Gene'].push(in.data['Gene'][[i]);
                out.data['feature'].push(in.data['feature'][i])
                out.data['value'].push(in.data['value'][i])
                out.data['Target'].push(in.data['Target'][i])
                out.data['Cancer'].push(in.data['Cancer'][i])
                out.data['type'].push(in.data['type'][i])
            }
        }
        return out
    }

    source = filter_by_target(source, target_type, cancer_type)

    function get_up_arr(out, in)
        for(var i=0, i<in.get_length(), i++){
            if((in.data['type'] == 'UPREG') || (in.data['type'] == 'GAIN')){
                out.data['Gene'].push(in.data['Gene'][[i]);
                out.data['feature'].push(in.data['feature'][i])
                out.data['value'].push(in.data['value'][i])
                out.data['Target'].push(in.data['Target'][i])
                out.data['Cancer'].push(in.data['Cancer'][i])
                out.data['type'].push(in.data['type'][i])
            }
        }
        return out
    }

    function get_neut_arr(in)
        for(var i=0, i<in.get_length(), i++){
            if((in.data['type'] == 'NEUT') || (in.data['type'] == 'NEUTRAL')){
                out.data['Gene'].push(in.data['Gene'][[i]);
                out.data['feature'].push(in.data['feature'][i])
                out.data['value'].push(in.data['value'][i])
                out.data['Target'].push(in.data['Target'][i])
                out.data['Cancer'].push(in.data['Cancer'][i])
                out.data['type'].push(in.data['type'][i])
            }
        }
        return out
    }

    function get_down_arr(in)
        for(var i=0, i<in.get_length(), i++){
            if((in.data['type'] == 'DOWNREG') || (in.data['type'] == 'LOSS')){
                out.data['Gene'].push(in.data['Gene'][[i]);
                out.data['feature'].push(in.data['feature'][i])
                out.data['value'].push(in.data['value'][i])
                out.data['Target'].push(in.data['Target'][i])
                out.data['Cancer'].push(in.data['Cancer'][i])
                out.data['type'].push(in.data['type'][i])
            }
        }
        return out
    }

    up = get_up_arr(source)
    neut = get_neut_arr(source)
    down = get_down_arr(source)

    up.change.emit()
    neut.change.emit()
    down.change.emit()

    // these values will be filled after filtering for the number of genes
    //var forreal_gene = [];
    //var forreal_feat = [];
    //var forreal_value = [];
    //var forreal_targ = [];
    //var forreal_canc = [];

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
    //var out = finalData.slice()
    //out.change.emit()
""")

# Drop down menu to choose cancer and target
cancer_type = Select(title="Cancer type:", value="Pan Cancer", options=["Breast Cancer", "Glioma", "Colorectal Cancer", "Lung Cancer", "Melanoma", "Renal Cancer", "Prostate Cancer", "Ovarian Cancer", "B-cell Lymphoma", "Pan Cancer"], callback=callback)
callback.args['cancer_type'] = cancer_type

target_type = Select(title="MetOncoFit Predictions:", value="Differential Expression", options=["Differential Expression", "Copy Number Variation", "Cancer Patient Survival"], callback=callback)
callback.args['target_type'] = target_type

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
