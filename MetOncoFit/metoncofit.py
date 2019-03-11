#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
`metoncofit.py` generates the main figures presented in Oruganty et al., containing:
  1) The top 10 biological features predicted by random forests
  2) The heatmap showing the relationship between the figures and the increased/netural/decreased gene-associated values
  3) The correlation value associated with each feature and differential gene expression, copy number variation, or cancer patient survival
  4) And the confusion matrix for model validation. The precision and recall values for the model are saved in separated excel files.
@authors: Krishna Oruganty & Scott Campit
"""

import sys
import process
import random_forest
import validator
import visualizations

import pandas as pd
from openpyxl import load_workbook

# Create data structures that will be used in the analysis
df, df1, header, canc, targ, data, classes, orig_data, orig_classes, excl_targ = process.preprocess('./../data/', sys.argv[1], sys.argv[2], sys.argv[3])

# Random Forest and Leave Out Validations
rfc, rfc_pred, mean_acc = random_forest.random_forest(canc, targ, data, classes, orig_data, orig_classes)

# Model performance and statistical measures
cm, pvalue, zscore, cv_score, summary = validator.summary_statistics(rfc, rfc_pred, data, classes, orig_classes, orig_data, targ, excl_targ, mean_acc, canc)

# Model comparison with Auslander et al., 2016. Use only gene expression in these predictions for a true comparison.
#genexp = df1[['NCI-60 gene expression',targ]]
#compare_models = validator.area_under_curve_calc(genexp, canc, targ)
#print(compare_models)

#df2 = df1.copy(deep=True)
# Leave-One-Cell Out Model Validation
#loco = validator.leave_one_cell_out(df2, canc, targ)

#df3 = df1.copy(deep=True)
# Leave-One-Feature-Set Out Model Validation
#lofo = validator.leave_one_feat_out(df3, canc, targ)

# Create data structures that will only be used while making the figures
up_df, neut_df, down_df, up_genes, neut_genes, down_genes, one_gene_df, one_gene_class = process.one_gene_only(df1, targ)

importance, up, neut, down, final_df = process.plotting_preprocess(up_df, neut_df, down_df, up_genes, neut_genes, down_genes, one_gene_df, rfc, header, targ, orig_classes, rfc_pred, one_gene_class)

# Make the figures
visualizations.make_figure(final_df, importance, cm, orig_classes, rfc_pred, cv_score, pvalue, zscore, canc, targ, normalize=True, savepath=False, filename=False)
"""
# Save the data into excel files
book = load_workbook('./../output/Tables/SI.xlsx')
writer = pd.ExcelWriter('./../output/Tables/SI.xlsx', engine='openpyxl')
writer.book = book

if targ == "CNV":
    summary.to_excel(writer, sheet_name="Table S4 CNV Summary Statistics")
elif targ == "SURV":
    summary.to_excel(writer, sheet_name="Table S5 SURV Summary Statistics")
elif targ == "DE":
    summary.to_excel(writer, sheet_name="Table S6 DE Summary Statistics")
elif targ == "TCGA_annot":
    summary.to_excel(writer, sheet_name="Table S7 TCGA Summary Statistics")
elif targ == "TCGA_annot_CNV":
    summary.to_excel(writer, sheet_name="Table S8 TCGA CNV Summary Statistics")

# Model validation, both internal and against Auslander et al
lofo.to_excel(writer, sheet_name="Table S9 Leave One Feature Set Out")
loco.to_excel(writer, sheet_name="Table S10 Leave One Cell Line Out")
compare_model.to_excel(writer, sheet_name="Table S11 AUC Performance")

writer.save()
"""
