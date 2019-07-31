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
from scipy.stats import pearsonr

from sklearn import preprocessing
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from imblearn.over_sampling import RandomOverSampler
from sklearn.model_selection import train_test_split


import prettify_df_labels


def label_encode(df):
    """
    """

    label_encoder = preprocessing.LabelEncoder()
    df["RECON1 subsystem"] = label_encoder.fit_transform(
        df["RECON1 subsystem"])
    df["Metabolic subnetwork"] = label_encoder.fit_transform(
        df["Metabolic subnetwork"])

    return df


def preprocess_input_data(file, prediction_target, exclude='exclusion'):
    """
    """

    if targ == 'TCGA_annot':
        targ = str("TCGA annotation")

    canc = prettify_df_labels.prettifyCancerLabels(file)
    column_names = prettify_df_labels.prettify_DataFrame_colNames()

    df = pd.read_csv(file, index_col=None, names=column_names).set_index(
        ['Genes', 'Cell Line'])
    df = label_encode(df)

    classes = df[targ]
    return df, classes


def remove_targets(df):
    """
    """
    excl_targ = {'TCGA annotation', 'SURV', 'CNV'}.remove(targ)

    if(len(sys.argv) > 3):
        fil3 = open("./../labels/"+exclude.rstrip())
        drop_col_names = [i.strip() for i in fil3.readlines()]
        fil3.close()
        df = df.drop(columns=drop_col_names)

    no_target_df = df.drop([excl_targ, targ])
    return no_target_df


def robust_scaler(df):
    """
    """
    data = np.array(df).astype(np.float)
    robust_scaled_data = RobustScaler().fit_transform(data)
    return robust_scaled_data


def split_data(robust_scaled_data, classes):
    """
    """
    training_data, test_data, training_labels, test_labels = train_test_split(
        robust_scaled_data, classes, test_size=0.3)

    random_overSampler = RandomOverSampler()
    oversampled_training_data, oversampled_training_labels = random_overSampler.fit_sample(
        training_data, training_labels)

    return oversampled_training_data, oversampled_training_labels, test_data, test_labels


def one_gene_only(df):
    """
    """

    prediction_labels, prediction_dict = prettify_df_labels.prettify_prediction_labels(
        prediction)

    df = df.reset_index()
    one_gene_df = df.drop(columns="Cell Line").groupby(
        ["Genes", targ]).median().reset_index().set_index("Genes")
    return one_gene_df


def

# These dataframes contain the df entries with increased, neutral, and decreased values.
up_df = one_gene_df.loc[(one_gene_df[targ] == targ_labels[0])]
neut_df = one_gene_df.loc[(one_gene_df[targ] == targ_labels[1])]
down_df = one_gene_df.loc[(one_gene_df[targ] == targ_labels[2])]

# To create the figure, we are randomly selecting three genes that are upreg, neutral, or downreg and are storing them in this list.
up_genes = up_df.index.values.tolist()
neut_genes = neut_df.index.values.tolist()
down_genes = down_df.index.values.tolist()

# Get rid of the targ column
_ = one_gene_df.pop(targ)


def compute_pearson_correlation():
    # This will calculate the correlation for each feature, if there is one between the biological features.
    pearson_correlation_coef = {}
    for feature in one_gene_df.columns:
        up_df_median = up_df[feature].median()
        neut_df_median = neut_df[feature].median()
        down_df_median = down_df[feature].median()
        perfect_correlation = [1.0, 0.0, -1.0]

        feature_correlation = np.pearsonr(
            [up_df_median, neut_df_median, down_df_median], perfect_correlation)

        if(np.isnan(feature_correlation[0]) != True):
            pearson_correlation_coef[feature] = feature_correlation[0]
        else:
            pearson_correlation_coef[feature] = 0.0
    return pearson_correlation_coef


def feature_importance_mapper(features, feature_importances):
    """
    idx_change sorts the feature importances and maps it to the feature name
    """
    feature_dictionary = {}
    for i, j in zip(features, feature_importances):
        feature_dictionary[i] = j
    sorted_dataframe = sorted(feature_dictionary.items(),
                              key=operator.itemgetter(1), reverse=True)
    return sorted_dataframe

    sorted_d = feature_importance_mapper(header, rfc.feature_importances_)

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
    supp_fig = importance.copy(deep=True)
    importance = importance.head(10)

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
    return importance, df
