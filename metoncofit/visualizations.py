#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The Visualizations module contains the code that can generate the confusion matrix, heatmap, clustermap, empirical cumulative distribution (ECD) plots, and dotplot that are in the manuscript

@author: Scott Campit
"""

import sys
import operator
import copy
import itertools
from random import shuffle
from itertools import chain
from math import pi

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from statsmodels.distributions.empirical_distribution import ECDF

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap


def make_figure(df1, importance, cm, orig_classes, rfc_pred, cv_acc, pval, zscore, canc, targ, normalize=True, savepath=False, filename=False, title_name=False):
    """
    This module creates files that replicate the figures in Oruganty et al. It will output the dotplot (left), the gini-importance barplot (middle), and the confusion matrix with normalized values (right).
    """

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    # Define target labels and colors based on the class
    class_col = {}
    if targ == 'CNV':
        targ_labels = ["GAIN", "NEUT", "LOSS"]
        class_col = {"GAIN": "#9C3848", "NEUT": "#808080", "LOSS": "#1E3888"}
    else:
        targ_labels = ["UPREG", "NEUTRAL", "DOWNREG"]
        class_col = {"UPREG": "#9C3848",
                     "NEUTRAL": "#808080", "DOWNREG": "#1E3888"}

    # Create confusion matrix and related variables
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    orig_classes = np.array(orig_classes)

    # Capture the correlation value symbols that will be used in the figure.
    pearson = importance['R'].tolist()
    correl = []
    for value in pearson:
        value = float(value)
        if value >= 0.60:
            correl.append(str('+'))
        elif value <= -0.60:
            correl.append(u"\u2014")
        else:
            correl.append('~')

    # Main figure parameters and arguments
    sns.set_style("whitegrid")
    figure, axarr = plt.subplots(nrows=1, ncols=3, figsize=(7.2, 3.6), gridspec_kw={
                                 'width_ratios': [1.0, 1.0, 0.75], 'wspace': 0.2}, sharex=False)

    # Legend parameters for the dotplot
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=list(class_col.values())[0], lw=3, markerfacecolor=list(
            class_col.values())[0], label=list(class_col.keys())[0]),
        Line2D([0], [0], color=list(class_col.values())[1], lw=3, markerfacecolor=list(
            class_col.values())[1], label=list(class_col.keys())[1]),
        Line2D([0], [0], color=list(class_col.values())[2], lw=3, markerfacecolor=list(class_col.values())[2], label=list(class_col.keys())[2])]

    # Show each observation in scatter plot
    sns.set_style("whitegrid")
    sns.stripplot(x="value", y="feature", hue="type", palette=class_col, data=df1,
                  order=importance['Feature'], hue_order=targ_labels, dodge=True, jitter=True, alpha=0.3, zorder=1, size=2.75, ax=axarr[0])

    from numpy import median
    # Show the conditional mean and standard deviation
    sns.pointplot(x="value", y="feature", hue="type", palette=class_col, data=df1,
                  order=importance['Feature'], hue_order=targ_labels, dodge=0.532, join=False, markers="D", scale=0.75, ci="sd", estimator=median, errwidth=1.00, ax=axarr[0])

    # Dotplot specific handles
    handles, labels = axarr[0].get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    axarr[0].legend(handles=legend_elements, bbox_to_anchor=(-1.50, -0.1),
                    loc=3, ncol=3, borderaxespad=0.0, frameon=False, fontsize=7)
    axarr[0].yaxis.label.set_visible(False)
    axarr[0].set_xlabel("Feature values", fontsize=7)
    axarr[0].set_ylabel("Top 10 important features", fontsize=7)
    axarr[0].tick_params(axis='x', which='both', top=False, labelsize=7)
    axarr[0].tick_params(axis='y', which='both', right=False, labelsize=7)
    axarr[0].set_xlim((-0.1, 1.1))
    axarr[0].grid(color='gray', axis='y')
    axarr[0].xaxis.grid(False)
    axarr[0].text(-2.25, 0, title_name, color='black', fontsize=7)

    # Create barplot
    sns.set(style="whitegrid")
    sns.barplot(x="Gini", y="Feature", data=importance,
                color='#1E3888', ax=axarr[1])

    # Barplot specific handles
    axarr[1].spines['right'].set_visible(False)
    axarr[1].spines['top'].set_visible(False)
    axarr[1].tick_params(axis='x', which='both', top=False, labelsize=7)
    axarr[1].tick_params(axis='y', which='both', right=False,
                         left=True, labelleft=False, labelsize=7)
    axarr[1].set_xlabel("Importance Score", fontsize=7)
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
        axarr[1].text(xmax+0.005, bar+0.05, correl[bar], color='black', fontsize=7,
                      bbox=dict(facecolor=bbox_col[bar], edgecolor=None), horizontalalignment='center')

    # Confusion matrix
    cmap = ListedColormap(sns.color_palette("Blues", 1000).as_hex())
    im = axarr[2].imshow(cm, interpolation='nearest',
                         cmap=cmap, vmin=0, vmax=1)
    #figure.colorbar(im, fraction=0.046, pad=0.04)
    plt.grid('off')

    tick_marks = np.arange(len(targ_labels))
    plt.setp(axarr[2], xticks=tick_marks, yticks=tick_marks,
             xticklabels=targ_labels, yticklabels=targ_labels)
    plt.tick_params(axis='x', top=False, labelsize=7)
    plt.tick_params(axis='y', left=False, labelright=True,
                    labelleft=False, labelsize=7)
    axarr[2].set_ylabel('Experimental group', labelpad=20,
                        fontsize=7).set_rotation(-90)
    axarr[2].yaxis.set_label_position("right")
    axarr[2].set_xlabel('MetOncoFit prediction', labelpad=10, fontsize=7)
    plt.xticks(rotation=45)
    axarr[2].xaxis.set_label_position("bottom")

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.0

    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        axarr[2].text(j, i, format(cm[i, j], fmt), fontsize=7,
                      horizontalalignment="center",
                      color="white" if cm[i, j] > thresh else "black")

    axarr[2].text(-0.4, 4.8, "Cross validation accuracy: "
                  + str(cv_acc)+'%', size=7, ha="left")
    axarr[2].text(-0.4, 5.0, "Z-score of accuracy: "
                  + (str(zscore).lstrip('[').rstrip(']')), size=7, ha="left")

    figure.savefig(savepath+'/'+filename+'.png', format='png',
                   dpi=300, bbox_inches='tight', pad_inches=0.2)


def specific_pathways_heatmap(df, importance, targ, canc, genelist, savepath=False, filename=False):
    """
    Get the genes associated with specific metabolic pathways you're interested in.
    """
    if targ == False:
        print("ERROR: Need to input targ set name.")

    if savepath == False:
        savepath = '.'

    if filename == False:
        filename = str(canc+'_'+targ)

    if targ == 'CNV':
        targ_labels = ["GAIN", "NEUT", "LOSS"]
    else:
        targ_labels = ["UPREG", "NEUTRAL", "DOWNREG"]

    # Read in gene list
    fil = open(genelist, "r")
    genes = fil.read().split('\n')
    del genes[0]
    del genes[0]

    # Get only the values you want:
    df = df.loc[df['Genes'].isin(genes)]
    rank = importance['Feature'].tolist()

    # Divide the df into 3 dataframes separated by target label
    up = df.loc[(df['type'] == targ_labels[0])]
    up = up.pivot(index='feature', columns='Genes', values='value')
    up = up.reindex(rank)
    neut = df.loc[(df['type'] == targ_labels[1])]
    neut = neut.pivot(index='feature', columns='Genes', values='value')
    neut = neut.reindex(rank)
    down = df.loc[(df['type'] == targ_labels[2])]
    down = down.pivot(index='feature', columns='Genes', values='value')
    down = down.reindex(rank)

    # Main figure parameters and arguments
    sns.set_style("whitegrid")
    figure, axarr = plt.subplots(nrows=1, ncols=3, figsize=(7.2, 3.6), gridspec_kw={
                                 'width_ratios': [1, 1, 1], 'wspace': 0.2}, sharex=False)

    # 3 heatmaps
    sns.heatmap(up, cmap='RdBu', robust=True, cbar=False,
                square=False, yticklabels=True, ax=axarr[0])
    sns.heatmap(neut, cmap='RdBu', robust=True, cbar=False,
                square=False, yticklabels=False, ax=axarr[1])
    sns.heatmap(down, cmap='RdBu', robust=True, cbar=True,
                square=False, yticklabels=False, ax=axarr[2])

    # Additional formatting for all 3 axes
    for i in range(0, 3):
        axarr[i].set_xlabel('')
        axarr[i].set_ylabel('')
        axarr[i].set_title(targ_labels[i])

    figure.tight_layout()
    figure.savefig(savepath+'/'+filename+'.png', format='png',
                   dpi=300, bbox_inches='tight', pad_inches=0.2)
