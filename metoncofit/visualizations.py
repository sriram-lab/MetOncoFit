#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The Visualizations module contains the code that can generate the confusion matrix, heatmap, clustermap, empirical cumulative distribution (ECD) plots, and dotplot that are in the manuscript

@author: Scott Campit
"""

import sys, operator, copy, itertools
from random import shuffle
from itertools import chain
from math import pi

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from statsmodels.distributions.empirical_distribution import ECDF

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

def conf_matr(orig_classes, pred_class, cv_acc, pval, zscore, targ, canc, normalize=True, cmap=False, savepath=False, filename=False):
    """
    The code to generate the confusion matrix was adapted from https://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html#sphx-glr-auto-examples-model-selection-plot-confusion-matrix-py
    """

    if cmap == False:
        cmap = ListedColormap(sns.color_palette("Blues", 1000).as_hex())

    if(targ == "CNV"):
        targ_labels = ["GAIN","NEUT","LOSS"]
    else:
        targ_labels = ["UPREG","NEUTRAL","DOWNREG"]

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    orig_classes = np.array(orig_classes)
    cm = confusion_matrix(orig_classes, pred_class)

    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    fig, ax = plt.subplots(figsize=(5,5))
    im = plt.imshow(cm, interpolation='nearest', cmap=cmap, vmin=0, vmax=1)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    tick_marks = np.arange(len(targ_labels))
    plt.xticks(tick_marks, targ_labels)
    plt.yticks(tick_marks, targ_labels)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.0

    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    plt.ylabel('True Class')
    plt.xlabel('Predicted Class')
    ax.xaxis.set_label_position('top')

    plt.text(-0.5, 2.9, "Cross validation accuracy (%): "+str(cv_acc), size=12, ha="left")
    plt.text(-0.5, 3.1, "p-value of accuracy: "+str(pval), size=12, ha="left")
    plt.text(-0.5, 3.3, "Z-score of accuracy: "+str(zscore), size=12, ha="left")
    plt.tight_layout()
    #plt.savefig(savepath+'/'+filename+'.svg', dpi=600)

    return cm

def plot_heatmap(df, one_gene_class, targ, canc, savepath=False, filename=False):
    """
    Create a heatmap that shows all of the genes and their corresponding feature values.
    """

    if targ == False:
        print("ERROR: Need to input targ set name.")

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    # Map each class to a specific color
    if targ == 'CNV':
        class_col = {"NEUT":"#808080", "LOSS": "#E69F00", "GAIN":"#009E73"}
    else:
        class_col = {"NEUTRAL":"#808080", "DOWNREG": "#E69F00", "UPREG":"#009E73"}
    col_heat = sns.color_palette("Reds")

    one_gene_class["Gene Class"] = one_gene_class[targ].map(class_col)
    one_gene_class = one_gene_class.drop(targ, axis=1).set_index("Gene")

    sns.heatmap(df, square=True, cmap = col_map, cbar=True)
    #plt.savefig(savepath+'/'+filename+'.svg', dpi=600)

def plot_clustermap(df, one_gene_class, targ, canc, method=False, metric=False, savepath=False, filename=False):
    """
    Create a Seaborn clustermap that shows all of the genes and their corresponding feature values. The method and metric parameters are optional arguments.
    """

    if method == False:
        method = 'single'

    if metric == False:
        metric = 'euclidean'

    if targ == False:
        print("ERROR: Need to input targ set name.")

    if canc == "breast":
        canc = "Breast"
    elif canc == "cns":
        canc = "CNS"
    elif canc == "colon":
        canc = "Colorectal"
    elif canc == "complex":
        canc = "Pan"
    elif canc == "leukemia":
        canc = "Leukemia"
    elif canc == "melanoma":
        canc = "Melanoma"
    elif canc == "nsclc":
        canc = "Lung"
    elif canc == "ovarian":
        canc = "Ovarian"
    elif canc == "prostate":
        canc = "Prostate"
    elif canc == "renal":
        canc = "Renal"
    else:
        print("ERROR: Need to input canc tissue name.")

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    # Map each class to a specific color
    if targ == 'CNV':
        class_col = {"NEUT":"#808080", "LOSS": "#E69F00", "GAIN":"#009E73"}
    else:
        class_col = {"NEUTRAL":"#808080", "DOWNREG": "#E69F00", "UPREG":"#009E73"}
    col_heat = sns.color_palette("Reds")

    one_gene_class["Gene Class"] = one_gene_class[targ].map(class_col)
    one_gene_class = one_gene_class.drop(targ, axis=1).set_index("Gene")

    sns.clustermap(df, method=method, metric=metric, figsize=(100,20), row_cluster=False, col_colors=one_gene_class, cmap = col_heat, robust=True)
    #plt.savefig(savepath+'/'+filename+'.svg', dpi=600)

def plot_ECD(up, neut, down, targ, canc, savepath=False, filename=False):
    """
    Create the empirical cumulative distribution plot as a function of the ranked feature values.
    """

    if targ == False:
        print("ERROR: Need to input targ set name.")

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    figure, axarr = plt.subplots(5, 2, figsize=(20,20), sharex=True, sharey=True)

    for k in range(0, len(features)):
        tmp_up = up[up['feature'] == features[k]]
        arr_up = tmp_up.sort_values('value', axis=0, ascending=True)

        tmp_neut = neut[neut['feature'] == features[k]]
        arr_neut = tmp_neut.sort_values('value', axis=0, ascending=True)

        tmp_down = down[down['feature'] == features[k]]
        arr_down = tmp_down.sort_values('value', axis=0, ascending=True)

        # Get individual ECDFs
        ecdf_up = np.linspace(0,1, len(tmp_up.index), endpoint=False)
        ecdf_neut = np.linspace(0,1, len(tmp_neut.index), endpoint=False)
        ecdf_down = np.linspace(0,1, len(tmp_down.index), endpoint=False)

        df_up = pd.DataFrame({"Feature value":arr_up['value'], "ECDF":ecdf_up, "Type":arr_up["type"]})
        df_neut = pd.DataFrame({"Feature value":arr_neut['value'], "ECDF":ecdf_neut, "Type":arr_neut["type"]})
        df_down = pd.DataFrame({"Feature value":arr_down['value'], "ECDF":ecdf_down, "Type":arr_down["type"]})

        # Arrange the ten subplots using 5x2 arrangement
        x = 0
        if k <= 4:
            x = k
        else:
            x = k - 5

        y = 0
        if k <= 4:
            y = 0
        else:
            y = 1

        # Plot everything onto a single figure
        sns.lineplot(x="Feature value", y="ECDF", color="#1B5258", data=df_up, ax=axarr[x, y])
        sns.lineplot(x="Feature value", y="ECDF", color="gray", data=df_neut, ax=axarr[x, y])
        sns.lineplot(x="Feature value", y="ECDF", color="#ECC655", data=df_down, ax=axarr[x, y])
        sns.lineplot(x=[0,1], y=[0,1], color='k', ax=axarr[x, y])

        # Figure properties
        axarr[x,y].set_title(features[k])
        axarr[x,y].set_xlabel('')
        axarr[x,y].set_ylabel('')
        figure.text(0.5, 0.04, 'Feature values', ha='center')
        figure.text(0.04, 0.5, 'Percentile', va='center', rotation='vertical')
        figure.suptitle("Predicting "+targ+" in "+canc+" canc", fontsize=20)

        # Legend properties
        import matplotlib.lines as mlines
        up_patch = mlines.Line2D([], [], color='red', label=class_col[0])
        neut_patch = mlines.Line2D([], [], color='gray', label=class_col[1])
        down_patch = mlines.Line2D([], [], color='blue', label=class_col[2])
        handles = [up_patch, neut_patch, down_patch]
        labels = [h.get_label() for h in handles]
        figure.legend(handles=handles, labels=labels, title="Gene Class")

        #plt.savefig(savepath+'/'+filename+'.svg', dpi=600)

def plot_importance(importance, targ, canc, savepath=False, filename=False):
    """
    Create the barplot that shows the feature importance values corresponding to the top 10 features for each canc. The correlation values are also
    """

    if targ == False:
        print("ERROR: Need to input targ set name.")

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    # Capture the correlation value symbols that will be used in the figure.
    pearson = importance['R'].tolist()
    correl = []
    for value in pearson:
        value = float(value)
        if value > 0.75:
            correl.append(str('+'))
        elif value < -0.75:
            correl.append(u"\u2013")
        else:
            correl.append('')

    # Create barplot
    plt.figure(figsize=(5,5))
    ax = sns.barplot(x="Gini", y="Feature", data=importance, color='#4B5E2D')

    # Figure parameters
    sns.set(style="white")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tick_params(axis='x', which='both', top=False)
    plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
    plt.xlabel("Gini Score")
    plt.ylabel('')

    # Create boxes that will be pearson correlation values
    _, xmax = plt.xlim()
    bbox_col = []
    for x in range(0, importance.shape[0]):
        if correl[x] == "+":
            bbox_col.append("#009E73")
        elif correl[x] == u"\u2013":
            bbox_col.append("#E69F00")
        else:
            bbox_col.append("#808080")

    for bar in range(0, importance.shape[0]):
        ax.text(xmax, bar+0.15, correl[bar], size='medium', color='black', bbox=dict(facecolor=bbox_col[bar], edgecolor=None), horizontalalignment='center')

    plt.tight_layout()
    #plt.savefig(savepath+'/'+filename+'.svg', dpi=600)

def plot_dotplot(df, targ, canc, savepath=False, filename=False):
    """
    Create a dotplot, where the large diamond is the mean, the bar is the standard deviation of the dataset, and the transparent diamonds are the datapoints.
    """
    if targ == False:
        print("ERROR: Need to input targ set name.")

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    class_col = {}
    # Map each class to a specific color
    if targ == 'CNV':
        class_col = {"GAIN":"#009E73", "NEUT":"#808080", "LOSS": "#E69F00"}
    else:
        class_col = {"UPREG":"#009E73", "NEUTRAL":"#808080", "DOWNREG": "#E69F00"}

    # Legend parameters
    from matplotlib.lines import Line2D
    legend_elements = [
    Line2D([0], [0], color=class_col.values()[0], lw=3, markerfacecolor=class_col.values()[0], label=class_col.keys()[0]),
    Line2D([0], [0], color=class_col.values()[1], lw=3, markerfacecolor=class_col.values()[1], label=class_col.keys()[1]),
    Line2D([0], [0], color=class_col.values()[2], lw=3, markerfacecolor=class_col.values()[2], label=class_col.keys()[2])]

    # Make the figure
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(15, 5), sharex=True, sharey=True)

    # Show each observation in scatter plot
    sns.stripplot(x="value", y="feature", hue="type", palette = class_col, data=df, dodge=True, jitter=True, alpha=0.3, zorder=1, size=2.75)
    sns.despine(offset=10, trim=True)

    # Show the conditional mean and standard deviation
    sns.pointplot(x="value", y="feature", hue="type", palette = class_col, data=df, dodge=0.532, join=False, markers="D", scale=0.75, ci="sd", errwidth=0.75)

    # Figure parameters

    plt.title("Top 10 feature value distribution for predicting "+targ+" in "+canc+" Cancer")
    plt.xlabel("Feature values")
    plt.ylabel("Top 10 important features")
    plt.tick_params(axis='x', which='both', top=False)
    plt.tick_params(axis='y', which='both', right=False)
    plt.xlim((-0.1, 1.1))
    plt.tight_layout()
    #plt.savefig(savepath+'/'+filename+'.svg', dpi=600)

def make_figure(df1, importance, cm, orig_classes, rfc_pred, cv_acc, pval, zscore, canc, targ, normalize=True, savepath=False, filename=False, title_name=False, analysis='main'):
    """
    This module creates files that replicate the figures in Oruganty et al. It will output the dotplot (left), the gini-importance barplot (middle), and the confusion matrix with normalized values (right).
    """

    import matplotlib.gridspec as gridspec
    from matplotlib import rcParams
    rcParams['text.usetex'] = True
    rcParams['text.latex.unicode'] = True
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica', 'Arial']

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    # Figure titles for main text
    if analysis == ('Main' | "SupplementaryPan"):
        if canc == "Breast Cancer":
            title_name = "$\bf{a. }"+canc
        elif canc == "Lung Cancer":
            title_name = "$\bf{b. }"+canc
        elif canc == "Skin Cancer":
            title_name = "$\bf{b. }"+canc
    elif analysis == ("Pan Cancer" | "SupplementaryCancer"):
        if targ == "SURV":
            title_name = "$\bf{c. )}$" + "Cancer Patient Survival"
        elif targ == "TCGA_annot":
            title_name = "$\bf{a)}$" + " Differential Expression"
        elif targ == "CNV":
            title_name = "$\bf{b)}$" + " Copy Number Variation"
    elif analysis == False:
        print("Default title is used (ie 'Renal Cancer')")
        title_name = canc

    # Define target labels and colors based on the class
    class_col = {}
    if targ == 'CNV':
        targ_labels = ["GAIN","NEUT","LOSS"]
        class_col = {"GAIN":"#9C3848", "NEUT":"#808080", "LOSS": "#1E3888"}
    else:
        targ_labels = ["UPREG","NEUTRAL","DOWNREG"]
        class_col = {"UPREG":"#9C3848", "NEUTRAL":"#808080", "DOWNREG": "#1E3888"}

    # Create confusion matrix and related variables
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    orig_classes = np.array(orig_classes)

    # Capture the correlation value symbols that will be used in the figure.
    pearson = importance['R'].tolist()
    correl = []
    for value in pearson:
        value = float(value)
        if value > 0.75:
            correl.append(str('+'))
        elif value < -0.75:
            correl.append(u"\u2014")
        else:
            correl.append('~')

    # Main figure parameters and arguments
    sns.set_style("whitegrid")
    figure, axarr = plt.subplots(nrows=1, ncols=3, figsize=(7.2,3.6), gridspec_kw={'width_ratios':[1.0,1.0,0.75], 'wspace':0.2}, sharex=False)

    # Legend parameters for the dotplot
    from matplotlib.lines import Line2D
    legend_elements = [
    Line2D([0],[0], color=list(class_col.values())[0], lw=3, markerfacecolor=list(class_col.values())[0], label=list(class_col.keys())[0]),
    Line2D([0],[0], color=list(class_col.values())[1], lw=3, markerfacecolor=list(class_col.values())[1], label=list(class_col.keys())[1]),
    Line2D([0],[0], color=list(class_col.values())[2], lw=3, markerfacecolor=list(class_col.values())[2], label=list(class_col.keys())[2])]

    # Show each observation in scatter plot
    sns.set_style("whitegrid")
    sns.stripplot(x="value", y="feature", hue="type", palette=class_col, data=df1, dodge=True, jitter=True, alpha=0.3, zorder=1, size=2.75, ax=axarr[0])

    # Show the conditional mean and standard deviation
    sns.pointplot(x="value", y="feature", hue="type", palette = class_col, data=df1, dodge=0.532, join=False, markers="D", scale=0.75, ci="sd", errwidth=1.00, ax=axarr[0])

    # Dotplot specific handles
    handles, labels = axarr[0].get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    axarr[0].legend(handles=legend_elements, bbox_to_anchor=(-1.50, -0.1), loc=3, ncol=3, borderaxespad=0.0, frameon=False, fontsize=7)
    axarr[0].yaxis.label.set_visible(False)
    axarr[0].set_xlabel("Feature values", fontsize=7)
    axarr[0].set_ylabel("Top 10 important features", fontsize=7)
    axarr[0].tick_params(axis='x', which='both', top=False, labelsize=7)
    axarr[0].tick_params(axis='y', which='both', right=False, labelsize=7)
    axarr[0].set_xlim((-0.1, 1.1))
    axarr[0].grid(color='gray', axis='y')
    axarr[0].xaxis.grid(False)
    axarr[0].text(-2.25,0, title_name, color='black', fontsize=7)

    # Create barplot
    sns.set(style="whitegrid")
    sns.barplot(x="Gini", y="Feature", data=importance, color='#1E3888', ax=axarr[1])

    # Barplot specific handles
    axarr[1].spines['right'].set_visible(False)
    axarr[1].spines['top'].set_visible(False)
    axarr[1].tick_params(axis='x', which='both', top=False, labelsize=7)
    axarr[1].tick_params(axis='y', which='both', right=False, left=True, labelleft=False, labelsize=7)
    axarr[1].set_xlabel("Gini Score", fontsize=7)
    axarr[1].set_ylabel('')
    axarr[1].grid(color='gray', axis='y')
    axarr[1].xaxis.grid(False)
    axarr[1].get_xaxis().set_major_locator(plt.MaxNLocator(5))

    # Create boxes that will be pearson correlation values
    _, xmax = axarr[1].set_xlim()
    bbox_col = []
    for x in range(0, importance.shape[0]):
        if correl[x] == "+":
            bbox_col.append("#FFAD69")
        elif correl[x] == u"\u2014":
            bbox_col.append("#47A8BD")
        else:
            bbox_col.append("#808080")

    # Create the +/- for the barplot
    for bar in range(0, importance.shape[0]):
        axarr[1].text(xmax+0.005, bar+0.05, correl[bar], color='black', fontsize=7, bbox=dict(facecolor=bbox_col[bar], edgecolor=None), horizontalalignment='center')

    # Confusion matrix
    cmap = ListedColormap(sns.color_palette("Blues", 1000).as_hex())
    im = axarr[2].imshow(cm, interpolation='nearest', cmap=cmap, vmin=0, vmax=1)
    #figure.colorbar(im, fraction=0.046, pad=0.04)
    plt.grid('off')

    tick_marks = np.arange(len(targ_labels))
    plt.setp(axarr[2], xticks=tick_marks, yticks=tick_marks, xticklabels=targ_labels, yticklabels=targ_labels)
    plt.tick_params(axis='x', top=False, labelsize=7)
    plt.tick_params(axis='y', left=False, labelright=True, labelleft=False, labelsize=7)
    axarr[2].set_ylabel('True Class', labelpad=20, fontsize=7).set_rotation(-90)
    axarr[2].yaxis.set_label_position("right")
    axarr[2].set_xlabel('Predicted Class', labelpad=10, fontsize=7)
    plt.xticks(rotation=45)
    axarr[2].xaxis.set_label_position("bottom")

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.0

    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        axarr[2].text(j, i, format(cm[i, j], fmt), fontsize=7,
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black")

    axarr[2].text(-0.4, 4.8, "Cross validation accuracy: "+str(cv_acc)+'%', size=7, ha="left")
    axarr[2].text(-0.4, 5.0, "Z-score of accuracy: "+(str(zscore).lstrip('[').rstrip(']')), size=7, ha="left")

    figure.savefig(savepath+'/'+filename+'.png', format='png', dpi=300, bbox_inches='tight', pad_inches=0.2)

def heatmap_html(df, savepath=False, filename=False):
    """
    Create the heatmap using Bokeh to create html files that will be embedded into the MetOncoFit website.
    """

    if savepath == False:
        savepath = '.'

    # default params == Pan Cancer and differential expression
    mask = (df.loc[(df['Cancer'] == "Pan") & (df["Target"] == "TCGA")])

    if df["Target"] == 'CNV':
        targ_labels = ["GAIN", "NEUT", "LOSS"]
    else:
        targ_labels = ["UPREGULATED", "NEUTRAL", "DOWNREGULATED"]

    # Input controls:
    # Slider to choose number of genes to show
    up_gene_select = Slider(start=0, end=len(df.loc[mask["type"] == targ_labels[0]]), value=5, step=5, title="Number of upregulated genes to display")
    neut_gene_select = Slider(start=0, end=len(df.loc[mask["type"] == targ_labels[1]]), value=5, step=1, title="Number of unregulated genes to display")
    down_gene_select = Slider(start=0, end=len(df.loc[mask["type"] == targ_labels[2]]), value=5, step=1, title="Number of downregulated genes to display")

    # Drop down menu to choose cancer and target
    cancer_type = Select(title="Cancer type:", value="Select cancer type", options=["Select cancer type", "Breast Cancer", "Glioma", "Colorectal Cancer", "Lung Cancer", "Melanoma", "Renal Cancer", "Prostate Cancer", "Ovarian Cancer", "B-cell Lymphoma", "Pan Cancer"])

    target_type = Select(title="MetOncoFit Predictions:", value="Select the cancer prognostic marker", options=["Select the cancer prognostic marker", "Differential Expression", "Copy Number Variation", "Cancer Patient Survival"])

    # Text input box for selecting a custom list of genes from the user
    gene_list = TextInput(value="Gene symbols (ie: IDH2, TDO2, PDGDH, ...)", title="Enter a list of gene symbols to query")

    # Create source that will be used by the plot
    source_up = ColumnDataSource(data=dict(Gene=[], feature=[], value=[]))
    source_neut = ColumnDataSource(data=dict(Gene=[], feature=[], value=[]))
    source_down = ColumnDataSource(data=dict(Gene=[], feature=[], value=[]))

    # Figure toolbar functions for interactions
    tools_in_figure = ["hover, save, pan, box_zoom, reset, wheel_zoom"]
    TOOLTIPS = [('Feature', '@feature'),('Gene', '@Gene'),('Value', '@value')]

    # Set up heat map figure spaces to be filled in
    hm1 = figure(title=targ_labels[0], x_axis_location='above', plot_height=400, plot_width=400, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
    hm2 = figure(title=targ_labels[1], x_axis_location='above', plot_height=400, plot_width=10000, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)
    hm3 = figure(title=targ_labels[2], x_axis_location='above', plot_height=400, plot_width=3000, tools=tools_in_figure, toolbar_location='right', tooltips=TOOLTIPS)

    # The actual figure vessels
    hm1.rect(x="Gene", y="feature", width=1, height=1, source=source_up, line_color=None, fill_color=transform('value', mapper))
    hm2.rect(x="Gene", y="feature", width=1, height=1, source=source_neut, line_color=None, fill_color=transform('value', mapper))
    hm3.rect(x="Gene", y="feature", width=1, height=1, source=source_down, line_color=None, fill_color=transform('value', mapper))

    # Custom subfunctions that will only be used in the scope of heatmap_html
    def usr_changes(df, attr, old, new):
        """
        This module will update the number of columns shown in the plot for a single heatmap based on user input via the button widget. If there are genes in the gene list field, it will take a list of genes from user input to show where on the heatmap they are. If there is no match, it needs to print somewhere, and just not show the gene that isn't in the dataset.
        """

        # Number of genes selected (Default = 5)
        up_gene_val = up_gene_select.value
        neut_gene_val = neut_gene_select.value
        down_gene_val = down_gene_select.value

        # The cancer type and target prediction
        canc_select = cancer_type.value.strip()
        target_select = target_type.value.strip()

        # First get the values corresponding to the cancer and target the user wants to query
        selected = df.loc[(df['Target'] == target_select) & (df['Cancer'] == canc_select)]

        # Second, separate into 3 dataframes that will be used to make the 3 heat maps
        up_select = selected.loc[(selected["type"] == targ_labels[0])]
        neut_select = selected.loc[(selected["type"] == targ_labels[1])]
        down_select = selected.loc[(selected["type"] == targ_labels[2])]

        # Third, get the number of genes they want to display.
        up_select = up_select.sample(up_gene_val)
        neut_select = neut_select.sample(neut_gene_val)
        down_select = down_select.sample(down_gene_val)

        # Fourth, if they specified some genes:
        if gene_list != ("Gene symbols (ie: IDH2, TDO2, PDGDH, ...)" or ""):
            up_select = up_select.loc[up_select["Gene"].str.contains(gene_list)==True]
            neut_select = neut_select.loc[neut_select["Gene"].str.contains(gene_list)==True]
            down_select = down_select.loc[down_select["Gene"].str.contains(gene_list)==True]

        # Finally, load dataframes into source that will be mapped back onto the figure
        source_up.data = dict(Gene=up_select["Gene"], feature=up_select["feature"], value=up_select["value"])
        source_neut.data = dict(Gene=neut_select["Gene"], feature=neut_select["feature"], value=neut_select["value"])
        source_down.data = dict(Gene=down_select["Gene"], feature=down_select["feature"], value=down_select["value"])

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

        # Color bar and params
        colors = brewer["RdBu"][8]
        mapper = LinearColorMapper(palette=colors, low=0, high=1)
        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7pt",  border_line_color=None, location=(0,0))
        p1.add_layout(color_bar, 'left')

    usr_changes()
    plots = gridplot([[hm1, hm2, hm3]], sizing_mode='fixed')
    html = file_html(plots, CDN, filename)
