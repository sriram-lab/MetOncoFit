"""
visualization.py
This module contains several functions that construct static figure objects generated in matplotlib and seaborn.

TODO:
  * Add Plotly versions of these functions
  * Add exploratory data analysis functions to visualize data distributions

@author: Scott Campit
"""

import itertools
from itertools import chain

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

def colormapper(targ):
    """
    colormapper specifies the color scheme for different labels depending on the the specific target labels used. For predicting differentially expressed genes or patient survival, the labels are UPREG, NEUTRAL, and DOWNREG. However, for predicting copy number alterations, the target labels are GAIN, NEUT, or LOSS.

    :param targ:         A string denoting the target type
    :return targ_labels: A list of strings for each class
    :return class_col    A dictionary containing the hex codes associated with each color
    """
    # Define target labels and colors based on the class
    if targ == 'CNV':
        targ_labels = ["GAIN", "NEUT", "LOSS"]
        class_col = {"GAIN": "#9C3848", "NEUT": "#808080", "LOSS": "#1E3888"}
    else:
        targ_labels = ["UPREG", "NEUTRAL", "DOWNREG"]
        class_col = {"UPREG": "#9C3848",
                     "NEUTRAL": "#808080", "DOWNREG": "#1E3888"}

    return targ_labels, class_col

def confusionMatrix(cm, targ, stats, normalize=True):
    """
    Create confusion matrix
    :param cm:        A numpy array containing the confusion matrix values.
    :param targ:      A string denoting the specific target for prediction.
    :param normalize: A boolean denoting whether to scale the confusion matrix to percentages
    :return ax:       A matplotlib object containing the confusion matrix plot data.
    """
    _, ax = plt.subplots()
    targ_labels, class_col = colormapper(targ)
    cmap = ListedColormap(sns.color_palette("Blues", 1000).as_hex())

    # Option to scale confusion matrix from 0 to 1 instead of raw true/false positive counts
    if normalize is True:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    fmt = '.2f' if normalize else 'd'

    # Confusion matrix
    _ = ax.imshow(cm, interpolation='nearest',
                         cmap=cmap, vmin=0, vmax=1)
    plt.grid('off')

    # Confusion matrix styling
    tick_marks = np.arange(len(targ_labels))
    plt.setp(ax, xticks=tick_marks, yticks=tick_marks,
                 xticklabels=targ_labels, yticklabels=targ_labels)

    plt.tick_params(axis='x', top=False, labelsize=7)
    plt.tick_params(axis='y', left=False,
                    labelright=True, labelleft=False,
                    labelsize=7)

    ax.set_xlabel('MetOncoFit prediction', labelpad=10, fontsize=7)
    ax.set_ylabel('Experimental group', labelpad=20,
                        fontsize=7).set_rotation(-90)
    ax.yaxis.set_label_position("right")
    ax.xaxis.set_label_position("bottom")
    plt.xticks(rotation=45)

    # Text formatting to either white or black, depending on the cell color
    thresh = cm.max() / 2.0
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        ax.text(j, i, format(cm[i, j], fmt), fontsize=7,
                      horizontalalignment="center",
                      color="white" if cm[i, j] > thresh else "black")

    # Also show the cross validation accuracy and Z-score of accuracy
    ax.text(-0.4, 4.8, "Cross validation accuracy: "
                  + str(stats[0]) + '%', size=7, ha="left")
    ax.text(-0.4, 5.0, "Z-score of accuracy: "
                  + (str(stats[1]).lstrip('[').rstrip(']')), size=7, ha="left")

    return ax

def dotplot(df, importance, targ):
    """
    Create dotplot showing data distribution
    :param df:         A pandas dataframe containing the raw values
    :param importance: A pandas dataframe containing an importance dataframe
    :param targ:       A string denoting the specific target for prediction.
    :return ax:        A matplotlib object containing the dot plot data.
    """
    from matplotlib.lines import Line2D
    from numpy import median

    targ_labels, class_col = colormapper(targ)

    # Legend parameters for the dotplot
    _, ax = plt.subplots()
    legend_elements = [
        Line2D([0], [0],
               color=list(class_col.values())[0],
               lw=3, markerfacecolor=list(class_col.values())[0],
               label=list(class_col.keys())[0]),
        Line2D([0], [0],
               color=list(class_col.values())[1],
               lw=3, markerfacecolor=list(class_col.values())[1],
               label=list(class_col.keys())[1]),
        Line2D([0], [0],
               color=list(class_col.values())[2],
               lw=3, markerfacecolor=list(class_col.values())[2],
               label=list(class_col.keys())[2])
    ]

    # Show each observation in scatter plot
    sns.set_style("whitegrid")
    sns.stripplot(x="value", y="feature", data=df,
                  hue="type", hue_order=targ_labels,
                  palette=class_col,
                  order=importance['Feature'],
                  dodge=True, jitter=True,
                  alpha=0.3, zorder=1, size=2.75, ax=ax)

    # Show the conditional median and standard deviation
    sns.pointplot(x="value", y="feature", data=df1,
                  hue="type", hue_order=targ_labels,
                  palette=class_col,
                  order=importance['Feature'],
                  dodge=0.532, join=False,
                  markers="D", scale=0.75, ci="sd",
                  estimator=median, errwidth=1.00, ax=ax)

    # Dotplot specific handles
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles=legend_elements,
              bbox_to_anchor=(-1.50, -0.1),
              loc=3, ncol=3,
              borderaxespad=0.0,
              frameon=False,
              fontsize=7)
    ax.yaxis.label.set_visible(False)
    ax.set_xlabel("Feature values", fontsize=7)
    ax.set_ylabel("Top 10 important features", fontsize=7)
    ax.tick_params(axis='x', which='both', top=False, labelsize=7)
    ax.tick_params(axis='y', which='both', right=False, labelsize=7)
    ax.set_xlim((-0.1, 1.1))
    ax.grid(color='gray', axis='y')
    ax.xaxis.grid(False)
    ax.text(-2.25, 0, title_name, color='black', fontsize=8)
    return ax

def variableImportance(importance):
    """
    Make variable importance box plots from random forest computation
    :oaram importance: A pandas dataframe containing important features from random forests
    :return ax:        A matplotlib object containing the variable importance plot
    """
    _, ax = plt.subplots()

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

    # Create barplot
    sns.set(style="whitegrid")
    sns.barplot(x="Gini", y="Feature", data=importance,
                color='#1E3888', ax=ax)

    # Barplot specific handles
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='x', which='both', top=False, labelsize=7)
    ax.tick_params(axis='y', which='both', right=False,
                         left=True, labelleft=False, labelsize=7)
    ax.set_xlabel("Importance Score", fontsize=7)
    ax.set_ylabel('')
    ax.grid(color='gray', axis='y')
    ax.xaxis.grid(False)
    ax.get_xaxis().set_major_locator(plt.MaxNLocator(5))

    # Create boxes that will be pearson correlation values
    _, xmax = ax.set_xlim()
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
        ax.text(xmax + 0.005, bar + 0.05, correl[bar], color='black', fontsize=7,
                bbox=dict(facecolor=bbox_col[bar], edgecolor=None), horizontalalignment='center')
    return ax


def concatFigures(cmPlt, dotPlt, varImpPlt):
    """
    This module creates files that replicate the figures in Oruganty et al. It will output the dotplot (left), the gini-importance barplot (middle), and the confusion matrix with normalized values (right).

    :param cmPlt:     A matplotlib object of the confusion matrix.
    :param dotPlt:    A matplotlib object of the dot plot.
    :param varImpPlt: A matplotlib object of the variable importance plot.

    :return axarr:    A matplotlib object containing the concatenated figure

    """

    # Main figure parameters and arguments
    sns.set_style("whitegrid")
    figure, axarr = plt.subplots(nrows=1, ncols=3, figsize=(7.2, 3.6),
                                 gridspec_kw={
                                              'width_ratios': [1.0, 1.0, 0.75],
                                              'wspace': 0.2
                                 }, 
                                 sharex=False)
    axarr[0] = dotPlt
    axarr[1] = varImpPlt
    axarr[2] = cmPlt

    return figure, axarr

def pathwayHeatmaps(df, importance, targ, genelist):
    """
    Create heatmap showing data values corresponding to specific pathways.
    :param df:         A pandas dataframe of the raw values
    :param importance: A pandas dataframe of the variable importances from random forests
    :param targ:       A string denoting the target type
    :param genelist:   A list containing the list of genes to query for the heatmap
    :return figure:    A matplotlib figure object
    :return axarr:     A matplotlib axes object
    """
    targ_labels, class_col = colormapper(targ)

    # Read in gene list
    file = open(genelist, "r")
    genes = file.read().split('\n')
    file.close()
    del genes[0]

    # Get only the values you want:
    df = df.loc[df['Genes'].isin(genes)]
    rank = importance['Feature'].tolist()

    # Main figure parameters and arguments
    sns.set_style("whitegrid")
    figure, axarr = plt.subplots(nrows=1, ncols=3,
                                 figsize=(7.2, 3.6),
                                 gridspec_kw={
                                     'width_ratios': [1, 1, 1],
                                     'wspace': 0.2
                                 }, sharex=False)

    # Divide the df into 3 dataframes separated by target label
    dfs = []
    for i in range(len(targ_labels)):
        dfs[i] = df.loc[(df['type'] == targ_labels[i])]
        dfs[i] = dfs[i].pivot(index='feature', columns='Genes', values='value')
        dfs[i] = dfs[i].reindex(rank)

        sns.heatmap(dfs[i], cmap='RdBu', robust=True, cbar=False,
                    square=False, yticklabels=True, ax=axarr[i])
        axarr[i].set_xlabel('')
        axarr[i].set_ylabel('')
        axarr[i].set_title(targ_labels[i])

    figure.tight_layout()
    return figure, axarr

if __name__ == "__main__":
    cmPlt = confusionMatrix(cm, targ, stats, normalize=True)
    dotPlt = dotplot(df, importance, targ)
    varPlt = variableImportance(importance)
    figure, axarr = concatFigures(cmPlt, dotPlt, varImpPlt)