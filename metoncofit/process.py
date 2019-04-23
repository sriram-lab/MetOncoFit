#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Process module contains functions that will be used as inputs for the random forest classifier and for generating the figures in the manuscript.

@authors: Krishna Dev Oruganty & Scott Campit
"""
import sys
import copy
import operator
import argparse
import re

import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from imblearn.over_sampling import RandomOverSampler
from sklearn.model_selection import train_test_split


def preprocess(datapath='', fil=sys.argv[1], targ=sys.argv[2], exclude=sys.argv[3]):
    """
    preprocess takes in the '*.csv' file and transforms the data that can be
    analyzed or fed into the MetOncoFit classifier.

    STEPS:
        1. Prettifies the column names from reading in the *.csv files
        2. Label encode network associations
        3. Removes TCGA or CNV values based on var_excl argument
        4. Robust scaling, split dataset, and random oversampling

    INPUTS:
        datapath: path to datasets for preprocessing
        file: csv file used for the analysis
        target: the target for random forest prediction
        exclude: specifies if the TCGA patient data will be included or excluded in the dataset

    OUTPUTS:
        df: DataFrame structure without the target classes. Should be used in the random_forest module
        df1: DataFrame structure containing the target class strings. Should be used in generating the visualizations
        canc: The string variable containing the cancer tissue
        targ: The string variable containing the target values
        var_excl: Features that should be removed in the analysis
        data, classes: The NumPy array containing scaled data and classes for training the random forest classifier
        orig_data, orig_classes: The NumPy array containing scaled data and classes for testing the random forest classifier
    """
    classes = []
    data = []
    names = []

    if targ == 'TCGA_annot':
        targ = str("TCGA annotation")

    if datapath is None:
        datapath = './../data/original'

    canc = fil.replace(".csv", "")
    canc_dict = {
        'breast':'Breast Cancer',
        'cns':'Glioma',
        'colon':'Colorectal Cancer',
        'complex':'Pan Cancer',
        'leukemia':'B-cell lymphoma',
        'melanoma':'Melanoma',
        'nsclc':'Lung Cancer',
        'ovarian':'Ovarian Cancer',
        'prostate':'Prostate Cancer',
        'renal':'Renal Cancer'
    }
    canc = canc_dict.get(canc)

    df_names = pd.read_csv("./../labels/real_headers.txt", sep='\t', names=['Original', 'New'])
    names = dict([(i, nam) for i, nam in zip(df_names['Original'], df_names['New'])])
    df = pd.read_csv(datapath+fil, index_col=None)
    df = df.rename(columns=names)
    df = df.set_index(['Genes','Cell Line'])

    # We are label encoding the subsystem and datapath labels
    le = preprocessing.LabelEncoder()
    df["RECON1 subsystem"] = le.fit_transform(df["RECON1 subsystem"])
    df["Biomass subsystem"] = le.fit_transform(df["Biomass subsystem"])

    excl_targ = {'TCGA annotation', 'SURV', 'CNV'}
    tmp = excl_targ.remove(targ)

    # We will drop a few columns in the cases where we have an exclusion list.
    if(len(sys.argv) > 3):
        fil3 = open("./../labels/"+exclude)
        drop_col_names = [i.strip() for i in fil3.readlines()]
        fil3.close()
        df = df.drop(columns=drop_col_names)

    df = df.drop(columns=excl_targ)
    classes = df[targ]
    header = df.columns
    df1 = df.copy(deep=True)  # contains target classes
    df = df.drop(columns=targ)  # doesn't contain target classes

    # Robust scaling the dataset with random oversampling
    data = np.array(df).astype(np.float)
    data = RobustScaler().fit_transform(data)

    new_data, orig_data, new_classes, orig_classes = train_test_split(
        data, classes, test_size=0.3)

    ros = RandomOverSampler()
    data, classes = ros.fit_sample(new_data, new_classes)

    return df, df1, header, canc, targ, data, classes, orig_data, orig_classes, excl_targ

def one_gene_only(df, target):
    """
    one_gene_only will merge the gene target value by majority rules and will take the median values for all numerical values.

    STEPS:
        1. Split gene and cell line index
        2. Capture unique genes, grouping by majority vote and taking the median value for the gene. The data is stored in one_gene_df
        3. Get specific up/neut/down-regulated gene information and store in dataframes / lists

    INPUTS:
        df: DataFrame structure from the preprocess function.
        target: The target we are going to predict (CNV, DE, SURV)

    OUTPUTS:
        one_gene_df, one_gene_class: Dataframe structure will all unique genes data and classes
        up_df, neut_df, down_df: DataFrame structure with data for genes
        up_genes, neut_genes, down_genes: List containing gene names
    """

    global up_df, neut_df, down_df, up_genes, neut_genes, down_genes, one_gene_df, one_gene_class

    if(target == "CNV"):
        targ_labels = ["GAIN", "NEUT", "LOSS"]
        targ_dict = {'NEUT': 0, 'LOSS': 0, 'GAIN': 0}
    else:
        targ_labels = ["UPREG", "NEUTRAL", "DOWNREG"]
        targ_dict = {'NEUTRAL': 0, 'DOWNREG': 0, 'UPREG': 0}

    df = df.reset_index()
    one_gene_df = df.drop(columns="Cell Line").groupby(["Genes", target]).median().reset_index().set_index("Genes")
    one_gene_class = pd.DataFrame(one_gene_df[target])
    one_gene_class = one_gene_class.reset_index()

    # These dataframes contain the df entries with increased, neutral, and decreased values.
    up_df = one_gene_df.loc[(one_gene_df[target] == targ_labels[0])]
    neut_df = one_gene_df.loc[(one_gene_df[target] == targ_labels[1])]
    down_df = one_gene_df.loc[(one_gene_df[target] == targ_labels[2])]

    # To create the figure, we are randomly selecting three genes that are upreg, neutral, or downreg and are storing them in this list.
    up_genes = up_df.index.values.tolist()
    neut_genes = neut_df.index.values.tolist()
    down_genes = down_df.index.values.tolist()

    return up_df, neut_df, down_df, up_genes, neut_genes, down_genes, one_gene_df, one_gene_class


def plotting_preprocess(up_df, neut_df, down_df, up_genes, neut_genes, down_genes, one_gene_df, rfc, header, targ, orig_classes, rfc_pred, one_gene_class, canc):
    """
    plotting_preprocess formats the data so that it can be used in the visualizations module.

    STEPS:
        1. Calculate the Pearson Correlation Coefficient for up/neut/downregulated dataframes
        2. Get the important features determined by the random forest classifier and capture the R value and Gini score in the importance dataframe.
        3. Scale the data using the min/max approach and create the DataFrame that will be used, based on the `plot` argument
    """

    global importance, up, neut, down, df

    # Remove the classes
    _ = one_gene_df.pop(targ)

    # This will calculate the correlation for each feature, if there is one between the biological features.
    column_squigly = {}
    for col in one_gene_df.columns:
        v1 = up_df[col].median()
        v2 = neut_df[col].median()
        v3 = down_df[col].median()

        correl = np.corrcoef([v1, v2, v3], [1.0, 0.0, -1.0])

        if(np.isnan(correl[0][1]) != True):
            column_squigly[col] = correl[0][1]
        else:
            column_squigly[col] = 0.0

    def idx_change(header, to_be_mapped):
        """
        idx_change sorts the feature importances and maps it to the feature name
        """
        temp_dict_feat = {}
        for i, j in zip(header, to_be_mapped):
            temp_dict_feat[i] = j
        sorted_d = sorted(temp_dict_feat.items(),
                          key=operator.itemgetter(1), reverse=True)
        return sorted_d
    sorted_d = idx_change(header, rfc.feature_importances_)

    feat = []
    gini = []
    corr = []

    x = 0
    while(x < 10):  # Get the first 10 features
        tempa = sorted_d[x]
        feat.append(tempa[0])
        gini.append(tempa[1])
        corr.append(str(column_squigly[tempa[0]]))
        x = x+1

    importance = pd.DataFrame({"Feature": feat, "Gini": gini, "R": corr})

    # Map to label
    if targ == 'CNV':
        class_col = ["GAIN", "NEUT", "LOSS"]
    else:
        class_col = ["UPREG", "NEUTRAL", "DOWNREG"]

    # Scale the dataframe from 0 to 1
    idx = one_gene_df.index
    col = one_gene_df.columns
    scaler = MinMaxScaler()
    result = scaler.fit_transform(one_gene_df)
    one_gene_df = pd.DataFrame(result, columns=col, index=idx)

    # Get the genes that are up/neut/downregulated
    tmparr_up = one_gene_df[one_gene_df.index.isin(up_genes)]
    tmparr_neut = one_gene_df[one_gene_df.index.isin(neut_genes)]
    tmparr_down = one_gene_df[one_gene_df.index.isin(down_genes)]

    features = list(importance['Feature'])
    up = tmparr_up[features].T
    up = up.reset_index().rename(columns={'index': "feature"})
    up = pd.melt(up, id_vars=["feature"])
    up["type"] = class_col[0]

    neut = tmparr_neut[features].T
    neut = neut.reset_index().rename(columns={'index': "feature"})
    neut = pd.melt(neut, id_vars=["feature"])
    neut["type"] = class_col[1]

    down = tmparr_down[features].T
    down = down.reset_index().rename(columns={'index': "feature"})
    down = pd.melt(down, id_vars=["feature"])
    down["type"] = class_col[2]

    df = pd.concat([up, neut, down], axis=0)
    df = df.sort_values('Genes')
    df['Cancer'] = canc
    df = df.reset_index().drop('index', axis=1)

    return importance, up, neut, down, df
