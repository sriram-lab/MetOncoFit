#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DataPreparation.py contains the functions that sets up the data to be fit into the MetOncoFit algorithm

@authors: Krishna Dev Oruganty & Scott Campit
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from sklearn import preprocessing
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from imblearn.over_sampling import RandomOverSampler
from sklearn.model_selection import train_test_split

import prettify_df_labels


def preprocess_input_data(model_file, predictionTarget_string, exclusionTargets):
    """
    preprocess_input_data reads in the data (from .csv files) for each cancer tissue model and outputs a pandas dataframe.
    """
    column_names = prettify_labels.long_feature_names()
    cancer = model_file.strip(".")[0]
    model = pd.read_csv(model_file, index_col=None, names=column_names).set_index(
        ['Genes', 'Cell Line'])

    return model, cancer


def label_encode(model):
    """
    label_encode uses the label_encoder from scikit-learn for the RECON1 subsystem and Metabolic subnetwork associations.
    """

    label_encoder = preprocessing.LabelEncoder()
    model["RECON1 subsystem"] = label_encoder.fit_transform(
        model["RECON1 subsystem"])
    model["Metabolic subnetwork"] = label_encoder.fit_transform(
        model["Metabolic subnetwork"])
    label_encoded_model = model.copy(deep=True)

    return label_encoded_model


def remove_target_from_training_data(label_encoded_model, target):
    """
    remove_target_from_training_data takes in an argument "target" and removes it from the label encoded model. The target labels are stored in the pandas Series "classes".
    """
    if target == 'TCGA_annot':
        target = str("TCGA annotation")
    excl_target = {'TCGA annotation', 'SURV', 'CNV'}.remove(target)

    classes = label_encoded_model[target]
    remove_list = open("./../labels/"+sys.argv[3].rstrip())
    drop_col_names = [i.strip() for i in remove_list.readlines()]
    remove_list.close()

    no_target_model = label_encoded_model.drop(columns=drop_col_names)
    no_target_model = no_target_model.drop([excl_target, target])

    return no_target_model, classes


def robust_scaler(no_target_model):
    """
    robust_scaler uses the scikit-learn Robust scaler function on the model data.
    """
    data = np.array(no_target_model).astype(np.float)
    robust_scaled_data = RobustScaler().fit_transform(data)

    return robust_scaled_data


def split_data_with_ROS(robust_scaled_data, classes):
    """
    split_data_with_ROS splits the data into 30% test 70% training with random oversampling to account for class imbalance.
    """
    training_data, test_data, training_labels, test_labels = train_test_split(
        robust_scaled_data, classes, test_size=0.3)

    random_overSampler = RandomOverSampler()
    oversampled_training_data, oversampled_training_labels = random_overSampler.fit_sample(
        training_data, training_labels)

    return oversampled_training_data, oversampled_training_labels, test_data, test_labels


def merge_cellLine_to_tissue_model(label_encoded_model, target):
    """
    merge_cellLine_to_tissue_model returns a tissue model containing single gene entries. The median values are used for the feature values.
    """

    label_encoded_model = label_encoded_model.reset_index()
    tissue_model = label_encoded_model.drop(columns="Cell Line")
    tissue_model = tissue_model.groupby(
        ["Genes", target]).median().reset_index()
    tissue_model = tissue_model.set_index(["Genes"])

    return tissue_model


def get_diffExp_genes(tissue_model, target_labels):
    """
    get_diffExp_genes is a function that returns a set of pandas dataframes that correspond to up / neutral / and down dataframes and gene lists based on the tissue model.
    """

    diffExp_dfs = {}
    diffExp_geneList = {}
    for label in range(1, len(target_labels)):
        diffExp_dfs[label] = tissue_model.loc[tissue_model[target_labels[label]]]
        diffExp_geneList[label] = diffExp_dfs[label].index.values.tolist()

    return diffExp_dfs, diffExp_geneList


def compute_pearson_correlation(diffExp_dfs, target):
    """
    compute_pearson_correlation returns a set of pearson correlation coefficients that correspond to the up / neutral / down dataframes mentioned previously.
    """
    _ = tissue_model.pop(target)
    pearson_correlation_coef = {}
    median = {}

    for label in range(1, len(diffExp_dfs)):
        for feature in diffExp_dfs[0].columns:
            median[label] = diffExp_dfs[label][feature].median()
            perfect_correlation = [1.0, 0.0, -1.0]
            feature_correlation = pearsonr(
                [median[0], median[1], median[2]], perfect_correlation)

            if np.isnan(feature_correlation[0]) != True:
                pearson_correlation_coef[feature] = feature_correlation[0]
            else:
                pearson_correlation_coef[feature] = 0.0

    return pearson_correlation_coef


def feature_importance_mapper(features, feature_importances):
    """
    feature_importance_mapper creates a sorted dataframe of all the features ranked by the importance score.
    """
    feature_dictionary = {}
    for feature, importance in zip(features, feature_importances):
        feature_dictionary[feature] = importance
    sorted_feature_df = sorted(
        feature_dictionary.items(),  key=operator.itemgetter(1), reverse=True)

    return sorted_feature_df


def get_importance_dataframe(sorted_feature_df):
    """
    get_importance_dataframe returns a sorted dataframe containing the feature, importance score, and associated pearson correlation coefficient.
    """
    feature = []
    importance_score = []
    pearson_correlation = []

    feat = 0
    while(feat < len(sorted_feature_df)):
        temp = sorted_feature_df[feat]
        feature.append(temp[0])
        importance_score.append(temp[1])
        pearson_correlation.append(str(pearson_correlation_coef[temp[0]]))
        feat = feat + 1

    importance_dataframe = pd.DataFrame({
        "Feature": feature,
        "Importance Score": importance_score,
        "R-Value": pearson_correlation
        })

    return importance_dataframe


def minMaxScale(tissue_model):
    """
    minMaxScale uses the scikit-learn function to perform min max scaling.
    """
    genes = tissue_model.index
    features = tissue_model.columns
    scaler = MinMaxScaler()
    minMax_tissueModel = scaler.fit_transform(tissue_model)
    minMax_tissueModel = pd.DataFrame(
        minMax_tissueModel, columns=features, index=genes)

    return minMax_tissueModel


def melt_dataframes(df, feature_list):
    """
    melt_dataframe performs a pd.melt to format the data for the figures.
    """
    df = df[feature_list].T
    df = df.reset_index().rename(columns={'Index': "Feature"})
    melted_df = df.melt(df, id_vars=["Feature"])

    return melted_df


def construct_figure_dataframes(minMax_tissueModel, importance, target, cancer):
    """
    construct_figure_dataframes returns a set of dataframes that will be used for the main figures in the manuscript.
    """
    diffExp_dfs, diffExp_geneList = get_diffExp_genes(minMax_tissueModel)
    top_10_feature = importance.head(10)

    prediction_labels, prediction_dict = prettify_df_labels.set_prediction_labels(
        target)
    feature_list = list(importance['Feature'])

    melted_dfs = {}

    for df in range(1, len(diffExp_dfs)):
        melted_dfs[df] = melt_dataframes(diffExp_dfs[df], feature_list)
        melted_dfs[df] = melted_dfs[df].sort_values('Genes')
        melted_dfs[df]["Label"] = prediction_labels[df]
        melted_dfs[df]["Cancer"] = cancer
        melted_dfs[df] = melted_dfs[df].reset_index().drop('index', axis=1)

    return melted_dfs
