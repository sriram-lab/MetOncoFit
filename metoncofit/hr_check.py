"""
HR Check
"""
import sys
import numpy as np
import scipy
import pandas as pd

import process
import random_forest
import validator
import visualizations
import save

# Get the frequency for each label in each dataset
df, df1, header, canc, targ, data, classes, orig_data, orig_classes, excl_targ, freq = process.preprocess(datapath=sys.argv[1], fil=sys.argv[2], targ=sys.argv[3], exclude=sys.argv[4])

# Random Forest
rfc, rfc_pred, mean_acc = random_forest.random_forest(canc, targ, data, classes, orig_data, orig_classes)

# Summary statistics
cm, pvalue, zscore, cv_score, summary = validator.summary_statistics(rfc, rfc_pred, data, classes, orig_classes, orig_data, targ, excl_targ, mean_acc, canc)

freq["10-fold CV Accuracy"] = cv_score
print(freq)

book = load_workbook('./../output/Tables/')
