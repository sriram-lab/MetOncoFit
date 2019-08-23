#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
`metoncofit.py` generates the following:
  1) The top 10 biological features from the 9 tumor models using Random Forests
  2) The heatmap showing the relationship between the figures and the increased/netural/decreased gene-associated values
  3) The correlation value associated with each feature and differential gene expression, copy number variation, or cancer patient survival
  4) And the confusion matrix for model validation. The precision and recall values for the model are saved in separated excel files.

@author: Scott Campit
"""

import sys
import process
import random_forest
import validator
import visualizations
import save

import pandas as pd
from openpyxl import load_workbook

# Create data structures that will be used in the analysis
df, df1, header, canc, targ, data, classes, orig_data, orig_classes, excl_targ, freq = process.preprocess(
    datapath='./../data/geneko/', fil=sys.argv[1], targ=sys.argv[2], exclude=sys.argv[3])

# Random Forest Classifier, prediction, and hold out accuracy
rfc, rfc_pred, mean_acc = random_forest.random_forest(
    canc, targ, data, classes, orig_data, orig_classes)

# Model performance and statistical measures. THIS FUNCTION IS ALSO NECESSARY TO GENERATE THE FIGURES.
cm, pvalue, zscore, cv_score, summary = validator.summary_statistics(
    rfc, rfc_pred, data, classes, orig_classes, orig_data, targ, excl_targ, mean_acc, canc)

# Model comparison with Auslander et al., 2016. Use only gene expression in these predictions for a true comparison.
#df2 = df1.copy(deep=True)
#compare_models = validator.area_under_curve_calc(df2, canc, targ)

# Leave-One-Cell Out Model Validation
#df3 = df1.copy(deep=True)
#loco = validator.leave_one_cell_out(df3, canc, targ)

# Leave-One-Feature-Set Out Model Validation
#df4 = df1.copy(deep=True)
#lofo = validator.leave_one_feat_out(df4, canc, targ)

# Save the data into supplementary figures:
#save.make_excel(summary, compare_models, loco, lofo, filename='SI.xlsx')

# Create data structures that will only be used while making the figures
importance, final_df = process.one_gene_only(df1, targ, header, rfc, canc)

# Make the figures
visualizations.make_figure(final_df, importance, cm, orig_classes, rfc_pred, cv_score,
                           pvalue, zscore, canc, targ, normalize=True, savepath=False, filename=False, title_name=canc)

# Make supplementary heatmaps
#fils = [
#  "alaaspglu.txt",
#  "argpro.txt",
#  "purine.txt",
#  "pyrimidine.txt",
#  "pyruvate.txt",
#  "tyrosine.txt"
#]

#genelist = r"/mnt/c/Users/scampit/Desktop/pyruvate.txt"
#visualizations.specific_pathways_heatmap(final_df, importance, targ, canc, genelist, savepath=False, filename=False)
