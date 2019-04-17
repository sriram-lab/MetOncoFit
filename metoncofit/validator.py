#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
`validator.py` contains several tools to assess the performance of the modelself.

Python Dependencies:
    * pandas
    * numpy
    * scipy
    * scikit-learn
    * imbalanced-learn

@authors: Krishna Oruganty & Scott Campit
"""

import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from random import shuffle
import scipy
from scipy import stats, interp

from sklearn.externals import joblib
from sklearn import preprocessing
from sklearn.metrics import cohen_kappa_score as coh_kap
from sklearn.metrics import f1_score, recall_score, precision_score, precision_recall_fscore_support, confusion_matrix, roc_curve, auc, classification_report, roc_auc_score, accuracy_score, average_precision_score, matthews_corrcoef
from sklearn.model_selection import train_test_split, permutation_test_score, cross_val_score
from sklearn.preprocessing import RobustScaler, label_binarize
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from imblearn.over_sampling import RandomOverSampler

def summary_statistics(rfc, rfc_pred, data, classes, orig_classes, orig_data, targ, excl_targ, mean_acc, canc):
    """
    summary_statistics takes in the classes and creates a precision probability distribution. This distribution is used to calculate several values of interest.

    INPUTS:
        clf_pred: the prediction array from the trained model.
        data: the training data.
        classes: the training label.
        orig_classes: the test labels.
        orig_data: the test data.

    OUTPUTS:
        cv_score: 10-fold cross validation score
        cm: confusion matrix (for visualization)
        ave: mean of the random precision score distribution.
        std: the standard deviation of the random precision score distribution.
        kap: Cohen's kappa coefficient
        f1: F1 score or harmonic average of the precision and recall
        prec: the average precision score across all classes
        rec: the average recall score across all classes
        prec_incr, prec_decr: the precision score for the upregulated/gain and downregulated/loss class
        rec_incr, rec_decr: the recall score for the upregulated/gain and downregulated/loss class
        zscore: the Z score of accuracy
        pvalue: the p-value of accuracy using a one-tailed T-test

    """

    # Create a distribution to measure the p-value associated with the accuracies obtained from MetOncoFit
    acc_distri = []
    x=0
    while(x <= 1000):
        temp_classes = list(classes)
        shuffle(temp_classes)
        acc_distri.append(precision_score(temp_classes, classes, average='micro'))
        x=x+1
    dist = np.array(acc_distri)

    # Basic statistical measures from the model
    ave = np.mean(np.array(acc_distri))
    std = np.std(dist)
    kap = coh_kap(orig_classes, rfc_pred)
    f1 = f1_score(orig_classes, rfc_pred, average='micro')
    cv_score = np.round(np.mean(cross_val_score(rfc, data, classes, scoring='accuracy', cv=10))*100, 2)
    mean_acc = np.round(mean_acc*100, 2)

    # Calculate the confusion matrix after 10 fold cross validation with new random forest classifier
    feat = 130
    while(feat < 140):
        trees = 5
        while(trees <= 500):
            x=0
            while(x < 10):
                X_train, X_test, y_train, y_test = train_test_split(data, classes, test_size=0.3)
                new_rfc = RandomForestClassifier(n_estimators=trees, max_features=feat)
                new_rfc.fit(X_train, y_train)
                cm_pred = new_rfc.predict(X_test)
                if x == 0:
                    cm = confusion_matrix(y_test, cm_pred)
                elif x > 1:
                    cm = np.add(cm, confusion_matrix(y_test, cm_pred))
                x=x+1
            trees = trees + 1500
        feat = feat + 20

    mcc = matthews_corrcoef(y_test, cm_pred)
    report = classification_report(y_test, cm_pred, output_dict=True)
    report = pd.DataFrame.from_dict(report).round(2)

    # Calculating the upper tailed p-value from the cumulative distribution function
    average_precision = report.loc[['precision'],['micro avg']].values[0]
    zscore = np.round((average_precision - ave)/std, 2)
    cdf = scipy.stats.norm.cdf(average_precision, ave, std)
    pvalue = 1 - cdf
    if(pvalue < 1e-50):
        pvalue = 1e-50

    if(targ == "CNV"):
        targ_labels = ["GAIN","NEUT","LOSS"]
    else:
        targ_labels = ["UPREG","NEUTRAL","DOWNREG"]
    dat = {"Average Hold-Out Accuracy (%)":mean_acc, "10-fold CV Accuracy (%)":cv_score, "Average Precision":(report.loc[['precision'],['micro avg']].values[0]), "UPREG/GAIN Precision":(report.loc[['precision'],[targ_labels[0]]].values[0]), "DOWNREG/LOSS Precision":(report.loc[['precision'],[targ_labels[2]]].values[0]), "Average Recall":(report.loc[['recall'],['micro avg']].values[0]), "UPREG/GAIN Recall":(report.loc[['recall'],[targ_labels[0]]].values[0]), "DOWNREG/LOSS Recall":(report.loc[['recall'],[targ_labels[2]]].values[0]), "P-value":pvalue, "Z-score":zscore, "F1-Score":(report.loc[['f1-score'],['micro avg']].values[0]), "MCC":mcc}

    summary = pd.DataFrame.from_dict(dict([(k, pd.Series(v)) for k,v in dat.items()])).T
    summary = summary.rename(columns={0:(canc+' Cancer')})
    summary = summary.T
    return cm, pvalue, zscore, cv_score, summary

def area_under_curve_calc(df1, canc, targ):
    """
    Calculates the score for the Area Under the Curve for MetOncoFit. This is to be compared to Auslander et al., 2016, specifically for

    INPUTS:
        data: the training data
        classes: the training targets
        orig_data: the test data
        orig_classes: the test targets

    OUTPUTS:
        auroc: the area under the receiever operating characteristic curve comparing the test dataset with the model's predictions.
    """

    if(targ == "CNV"):
        targ_labels = ["GAIN","NEUT","LOSS"]
    else:
        targ_labels = ["UPREG","NEUTRAL","DOWNREG"]

    # Pop out the classification labels & binarize the array. This is necessary for calculating the auroc using sklearn.
    cls = df1.pop(targ)
    cls = label_binarize(cls, classes=[targ_labels[0], targ_labels[1], targ_labels[2]])

    data = np.array(df1).astype(np.float)
    #data = RobustScaler().fit_transform(data)

    # MetOncoFit Random Forest Classifier with specifications to the method used by Auslander et al
    new_data, orig_data, new_classes, orig_classes = train_test_split(data, cls, test_size=0.2)
    feat = (data.shape[1]-10)
    while(feat < data.shape[1]-1):
        trees = 5
        while(trees <= 500):
            rfc = RandomForestClassifier()
            rfc.fit(new_data, new_classes)
            trees = trees + 1500
        feat = feat + 20
    rfc_pred = rfc.predict(orig_data)
    #rfc_score = np.asarray(rfc.predict_proba(orig_data)[:,1])
    #rfc_score = np.transpose(rfc_score[:,:,0])

    # Calculate AUROC, averaged by taking into consideration label imbalance
    fpr, tpr, _ = roc_curve(orig_classes.ravel(), rfc_pred.ravel())
    auroc = auc(fpr, tpr)
    dat = {"Cancer":[canc], "Target":[targ], "AUROC":[auroc]}
    df = pd.DataFrame.from_dict(dat)
    return df

def leave_one_feat_out(df, canc, targ):
    """
    Leave one feature out reports the accuracy obtained from removing the following features:

    1. Topological features only
    2. Dynamic features only
    3. Expression only
    4. Expression and kcat
    """

    # create 4 dataframes that will be used for each of the features after robust scaler
    df1 = df.copy(deep=True)
    classes = df1.pop(targ)
    num = RobustScaler().fit_transform(np.array(df1).astype(np.float))
    df1 = pd.DataFrame(num, columns=df1.columns, index=df1.index)

    dynm = df1.drop(df1.columns[0:52], axis=1)
    topo = df1.drop(df1.columns[53:132], axis=1)
    kexp = df1.drop(df1.columns[132:], axis=1)
    genexp = df1.drop(df1.columns[133:], axis=1)
    subsys = df1.drop(df1.columns[134], axis=1)

    # concatenate the dataframes to a single structure and get accuracies
    dfs = [topo, dynm, kexp, genexp]

    output = []
    for df in dfs:
        if df is topo:
            lofo = "Toplogical Features"
        elif df is dynm:
            lofo = "Dynamic Features"
        elif df is kexp:
            lofo = "Gene expression and kcat"
        elif df is genexp:
            lofo = "Gene expression only"
        elif df is subsys:
            lofo = "RECON1 Subsystem only"
        else:
            return("ERROR: Not suitable df input")

        new_data, orig_data, new_classes, orig_classes = train_test_split(df, classes, test_size=0.3)

        feat = (new_data.shape[1]-10)
        while(feat < new_data.shape[1]-1):
            trees = 5
            while(trees <= 500):
                rfc = RandomForestClassifier(n_estimators=trees, max_features=feat)
                rfc.fit(new_data, new_classes)
                trees = trees + 1500
            feat = feat + 20
            rfc_pred = rfc.predict(orig_data)
            mean_acc = rfc.score(orig_data, orig_classes)
            output.append([canc, targ, lofo, mean_acc])

    # Return data frame to be saved
    lofo_df = pd.DataFrame(output, columns=["Cancer", "Target", "Held-out feature set", "Mean class accuracy"])
    return lofo_df

def leave_one_cell_out(df2, canc, targ):
    """
    Leave one cell out outputs the mean accuracy obtained after holding out a single NCI-60 cancer cell line from the dataset.

    """

    # Split the index into Gene symbol and Cell Line
    df2 = df2.reset_index()
    _ = df2.pop('Genes') # don't need gene labels
    cell_line = df2.pop("Cell Line")
    cell_line = pd.DataFrame(cell_line)

    # Temporarily remove the classes to do robust scaling
    classes = df2.pop(targ)
    classes = pd.DataFrame(classes, index=df2.index)

    # Robust scaling (since this is not done in the process script)
    num = RobustScaler().fit_transform(np.array(df2).astype(np.float))
    df2 = pd.DataFrame(num, columns=df2.columns, index=df2.index)

    # Unite the datasets together again
    df2 = pd.concat([cell_line, df2, classes], axis=1)
    groups = df2["Cell Line"].unique()

    # Leave one cell line out
    output = []
    for cell in groups:
        new_df = df2[~df2["Cell Line"].str.contains(str(cell))]

        # This time pop out for sure
        classes = new_df.pop(targ)
        _ = new_df.pop("Cell Line")

        # Train test and split as usual
        new_data, orig_data, new_classes, orig_classes = train_test_split(new_df, classes, test_size=0.3)

        # Now we can do random oversampling
        ros = RandomOverSampler()
        data, classes = ros.fit_sample(new_data, new_classes)

        feat = (data.shape[1]-10)
        while(feat < data.shape[1]-1):
            trees = 5
            while(trees <= 500):
                rfc = RandomForestClassifier(n_estimators=trees, max_features=feat)
                rfc.fit(data, classes)
                trees = trees + 1500
            feat = feat + 20
            rfc_pred = rfc.predict(orig_data)
            mean_acc = rfc.score(orig_data, orig_classes)
            output.append([canc, targ, cell, mean_acc])

    # Return data frame to be saved
    loco = pd.DataFrame(output, columns=["Cancer", "Target", "Held-out cell line", "Mean class accuracy"])
    return loco
