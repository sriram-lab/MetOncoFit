#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DataPreparation.py contains the functions that sets up the data to be fit into the MetOncoFit algorithm

@authors: Krishna Dev Oruganty & Scott Edward Campit
"""
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn import preprocessing

def load_data(model_file, labelFileName):
    """
    load_data reads in the cancer model data (.csv file) and outputs a pandas dataframe.

    :params:
        model_file: The path to the .csv file containing the rows as observations and the columns as features.
            Note: there needs to be a corresponding 'Genes' and 'Cell Line' column to set as the index.

    :return:
        model:      A pandas dataframe containing the cancer model data, with observations as rows and
            features as columns.
        cancer:     A string denoting the tissue type from the name of the .csv file.
    """

    import PrettifyLabels

    column_names = PrettifyLabels.long_feature_names(labelFileName)
    cancer = model_file.strip(".")[0]
    model = pd.read_csv(model_file)
    model = model.rename(columns=column_names)
    model = model.set_index(['Genes', 'Cell Line'])

    return model, cancer


def label_encode(model):
    """
    label_encode uses the label_encoder function from scikit-learn for the RECON1 subsystem and Metabolic subnetwork
        features. Note that these features may be removed in future MetOncoFit versions.

    :params:
        model:               A pandas dataframe containing the cancer model data, with observations as rows and
            features as columns.

    :return:
        label_encoded_model: A panda dataframe of the label-encoded model.
    """

    label_encoder = preprocessing.LabelEncoder()
    model['RECON1 subsystem'] = label_encoder.fit_transform(
        model['RECON1 subsystem'])
    model['Metabolic subnetwork'] = label_encoder.fit_transform(
        model['Metabolic subnetwork'])
    label_encoded_model = model.copy(deep=True)

    return label_encoded_model


def prune_targets(model, target="DE", exclude="DE_and_CNV"):
    """
    prune_targets removes values that determined the target labels from the label encoded model.

    target (optional):  A string denoting the specific target variable to train the model on. The default argument
        is 'DE'. There are three arguments:
        - 'DE': Predict differential metabolic enzyme expression.
        - 'CNV': Predict amplification or loss of copy number variation for a given metabolic enzyme.
        - 'SURV' Predict whether a metabolic enzyme up-regulation or down-regulation is associated with patient
           mortality.

    exclude (optional): A string denoting the argument that removes columns corresponding to the raw feature values
        that are also used to determine the target labels. For example, we can make predictions on copy number
        variation while using gene expression values as features for training the cancer model, or remove the gene
        expression values entirely. The default argument is 'DE_and_CNV'.
        There are two arguments:
        - 'DE_and_CNV': Removes fold change values for differential expression and copy number variation.
        - 'CNV_only':   Removes fold change values for copy number variation only.
    """
    label_encoded_model = model.copy(deep=True)
    if exclude is 'DE_and_CNV':
        label_encoded_model = label_encoded_model.drop(['TCGA gene expression fold change',
                                                        'CNV gain/loss ratio'])
    elif exclude is 'CNV_only':
        label_encoded_model = label_encoded_model.drop(['CNV gain/loss ratio'])

    targetVariables = {'DE':'TCGA annotation',
                       'CNV':'CNV',
                       'SURV':'SURV'}
    target = targetVariables.get(target)
    classes = label_encoded_model[target]
    pruned_model = label_encoded_model.drop(
        ["TCGA annotation", "CNV", "SURV"],
        axis=1)

    return pruned_model, classes


def robust_scaler(model):
    """
    robust_scaler uses the scikit-learn RobustScaler function to scale the data using the interquartile ranges.

    :params:
        model: A pandas dataframe containing the model without the target labels

    :return:
        robust_model: A pandas dataframe containing the model that has been standardized by the IQR
    """
    from sklearn.preprocessing import RobustScaler

    data = np.array(model).astype(np.float)
    robust_model = RobustScaler(with_centering=True, with_scaling=True).fit_transform(data)

    return robust_model


def randomOversampling(model, classes, testSize=0.2):
    """
    randomOversampling takes a pandas dataframe and attempts to perform naive random oversampling on classes that are
    naturally under-represented in the dataset.

    :params:
        model:              A pandas dataframe containing the model without the target labels
        classes:            A pandas series containing the labels for a specific target variable
        testSize(optional): A float value corresponding to the size of the test dataset. The default value is 20%.

    :return:
        Xtrain:  A pandas dataframe containing oversampled data used to train the model.
        Xtest:   A pandas dataframe containing the test dataset.
        Ytrain:  A pandas series containing oversampled labels used to train the model.
        Ytest:   A pandas series containing the test labels.
        
    """
    from imblearn.over_sampling import RandomOverSampler
    from sklearn.model_selection import train_test_split
    
    Xtrain, Xtest, Ytrain, Ytest = train_test_split(model, classes,
                                                    test_size=testSize,
                                                    train_size=1-testSize,
                                                    random_state=1,
                                                    shuffle=True)

    over_sampler = RandomOverSampler(sampling_strategy='auto',
                                     random_state=1)

    Xtrain, Ytrain = over_sampler.fit_sample(Xtrain, Ytrain)

    return Xtrain, Xtest, Ytrain, Ytest

def processDataFromFile(filename, target, exclude, labelFileName):
    model, cancer = load_data(filename, labelFileName)
    labelEncodedModel = label_encode(model)
    prunedModels, classes = prune_targets(labelEncodedModel, target, exclude)
    robustModel = robust_scaler(prunedModels)
    Xtrain, Xtest, Ytrain, Ytest = randomOversampling(robustModel, classes, testSize=0.2)
    return Xtrain, Xtest, Ytrain, Ytest

def create_tissue_model(model, target):
    """
    create_tissue_model returns a tissue model containing single gene entries. The median values corresponding to each
    observation are used for the final feature values.

    :params:
        model:  A pandas dataframe containing the tumor dataset, ideally after standardization / scaling / sampling.
        target: A pandas series containing the labels for the target variable.

    :return:
        tissue_model: A pandas dataframe representing the tumor model, which contains the median values for each gene.
    """

    model = model.reset_index()
    tissue_model = model.drop(columns=["Cell Line"])
    tissue_model = tissue_model.groupby(
        ["Genes", target]).median().reset_index()
    tissue_model = tissue_model.set_index(["Genes"])

    return tissue_model


def DE_genes(model, target):
    """
    DE_genes returns pandas dataframes that correspond to up / neutral / and down dataframes and gene lists based on
    the input dataframe.

    :params:
        model:  A pandas dataframe containing n observerations by p predictors. The tumor model is the ideal pandas
            dataframe to input into this function.
        target: A pandas series containing the labels corresponding to a given target variable.

    :return:
        diffExpDFs:   A dictionary containing dataframes corresponding to the unique labels within the target variable.
        diffExpGenes: A dictionary containing gene lists corresponding to the unique labels within the target variable.

    """

    diffExpDFs = {}
    diffExpGenes = {}
    for label in range(0, len(target)):
        diffExpDFs[label] = model.loc[model[target[label]]]
        diffExpGenes[label] = diffExpDFs[label].index.values.tolist()

    return diffExpDFs, diffExpGenes


def feature_importance_map(features, feature_importances):
    """
    feature_importance_map creates a sorted dataframe of all the features ranked by the importance score.

    :params:
        features:            A panadas series containing the features in the model
        feature_importances: A pandas series containing the Gini impurity index corresponding the each feature.

    :return:
        sorted_feature_df:   A pandas dataframe containing the features ranked by the importance score.

    """
    feature_dictionary = {}
    for feature, importance in zip(features, feature_importances):
        feature_dictionary[feature] = importance

    sorted_feature_df = sorted(
        feature_dictionary.items(),  key=operator.itemgetter(1), reverse=True
    )

    return sorted_feature_df


def get_importance_dataframe(sorted_feature_df, pearsonCorrelationDict):
    """
    get_importance_dataframe returns a sorted dataframe containing the feature, importance score, and associated
    pearson correlation coefficient.

    :params:
        sorted_feature_df: A pandas dataframe containing the features ranked by the importance score.

    :return:
        importanceDF: A pandas dataframe containing the feature, importance score, and pearson correlation value.

    """
    feature = []
    importance_score = []
    pearson_correlation = []

    count = 0
    while(count < len(sorted_feature_df)):
        temp = sorted_feature_df[feat]
        feature.append(temp[0])
        importance_score.append(temp[1])
        pearson_correlation.append(str(pearsonCorrelationDict[temp[0]]))
        count += count + 1

    importanceDF = pd.DataFrame({
        "Feature": feature,
        "Importance Score": importance_score,
        "R-Value": pearson_correlation
        })

    return importanceDF


def minMaxScale(model):
    """
    minMaxScale uses the scikit-learn function to perform min max scaling for each feature.

    :params:
        model:       A pandas dataframe. This function is intended on scaling the features of the tissue model.

    :return:
        scaledModel: A pandas dataframe containing the data scaled with values from 0 to 1.
    """
    from sklearn.preprocessing import MinMaxScaler

    genes = model.index
    features = model.columns
    mdl = pd.DataFrame()

    for feat in features:
        scaler = MinMaxScaler(feature_range=(0,1))
        mdl[feat] = scaler.fit_transform(model[feat])

    scaledModel = pd.DataFrame(
        mdl, columns=features, index=genes
    )
    return scaledModel


def melt_dataframe(df, feature_list):
    """
    melt_dataframe performs a pd.melt to format the data for the figures.

    :params:
        df:           A pandas dataframe corresponding to a tumor model.
        feature_list: A list containing the feature names.

    :return:
        melted_df:    A pandas dataframe reformatted for the figures.

    """
    df = df[feature_list].T
    df = df.reset_index().rename(columns={'Index': "Feature"})
    melted_df = df.melt(df, id_vars=["Feature"])

    return melted_df


def constructFigureDF(model, importance, target, cancer):
    """
    constructFigureDF returns a set of dataframes that will be used for the main figures in the manuscript.

    :params:
        model:      A pandas dataframe with the preprocessed data.
        importance: A pandas dataframe containing the features, importance, and pearson correlation coefficients
        target:     A string denoting the target variable of interest
        cancer:     A string denoting the tumor model

    :return:
        meltedDF:   A pandas dataframe consisting of the melted features

    """
    diffExpDFs, diffExpGeneList = DE_genes(model)
    top10Features = importance.head(10)

    predictionLabels, predictionDict = prettify_df_labels.set_prediction_labels(
        target)
    featureList = list(importance['Feature'])

    meltedDFs = {}

    for df in diffExp_dfs:
        meltedDFs[df] = melt_dataframes(diffExpDFs[df], featureList)
        meltedDFs[df] = meltedDFs[df].sort_values('Genes')
        meltedDFs[df]["Label"] = predictionLabels[df]
        meltedDFs[df]["Cancer"] = cancer
        meltedDFs[df] = meltedDFs[df].reset_index().drop('index', axis=1)

    return meltedDFs

if __name__ == '__main__':
    pass
