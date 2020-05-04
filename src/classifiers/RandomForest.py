#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Random Forest module for training the random forest classifier and outputting accuracy scores. There are 5 functions:
    1. random_forest: train the random forest classifier
    2. save_model: pickle the random forest classifier
    3. load_model: load the pickled random forest classifier
@authors: Krishna Dev Oruganty & Scott Campit
"""
import os, sys

import numpy as np
import pandas as pd

from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

def random_forest(canc, targ, data, classes, orig_data, orig_classes):
    """
    random_forest will train the random forest classifier using the pre-processed data.
    """

    global rfc, rfc_pred, mean_acc

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

    return rfc, rfc_pred, mean_acc

def save_model(canc, targ, var_excl, clf):
    """
    cPickle the random forest classification model for future use to ensure robustness.
    """
    # Pickle model
    save_path = './../models/'
    filename = os.path.join(save_path, (canc+'_'+targ+"_"+var_excl+'_model.pkl'))
    joblib.dump(clf, filename)

def load_model(model, orig_data):
    """
    Load the pickled random forest model to use.
    """
    clf = joblib.load(model)
    prediction = clf.predict(orig_data)

    return(prediction)
