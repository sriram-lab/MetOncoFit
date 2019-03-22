#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
db.py creates the metoncofit dataframe that will be used to create html figures 
"""

import sys, copy, operator, argparse, re, os

import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from imblearn.over_sampling import RandomOverSampler
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

global fil, t, var_excl
datapath = None
all_dfs = []
targ = ["TCGA_annot", "CNV", "SURV"]
var_excl = ["TCGA gene expression fold change", "CNV gain/loss ratio"]
for fil in os.listdir('./../data'):
    for t in targ:
        if datapath is None:
            datapath = './../data/'
        canc = fil.replace(".train.csv","")
        if canc == "breast":
            canc = "Breast"
        elif canc == "cns":
            canc = "CNS"
        elif canc == "colon":
            canc = "Colorectal"
        elif canc == "complex":
            canc = "Pan"
        elif canc == "leukemia":
            canc = "Leukemia"
        elif canc == "melanoma":
            canc = "Melanoma"
        elif canc == "nsclc":
            canc = "Lung"
        elif canc == "ovarian":
            canc = "Ovarian"
        elif canc == "prostate":
            canc = "Prostate"
        elif canc == "renal":
            canc = "Renal"

        classes = []
        data = []
        names = []

        if t == 'TCGA_annot':
            t = str("TCGA annotation")

        df_names = pd.read_csv("./../labels/real_headers.txt", sep='\t')
        df = pd.read_csv(datapath+fil, names=list(df_names.iloc[:,1]), index_col=0, skiprows=1)

        # We are label encoding the subsystem and datapath labels
        le = preprocessing.LabelEncoder()
        df["RECON1 subsystem"] = le.fit_transform(df["RECON1 subsystem"])
        df["Biomass subsystem"] = le.fit_transform(df["Biomass subsystem"])

        excl_targ = {'TCGA annotation', 'SURV', 'CNV'}
        tmp = excl_targ.remove(t)

        # We will drop a few columns in the cases where we have an exclusion list.
        if(len(sys.argv) > 3):
            fil2=open("./../labels/"+var_excl)
            drop_col_names = [i.strip() for i in fil2.readlines()]
            fil2.close()
            df = df.drop(columns=drop_col_names)

        df = df.drop(columns=excl_targ)
        classes = df[t]
        header = df.columns
        df1 = df.copy(deep=True) # contains target classes
        df = df.drop(columns=t) # doesn't contain target classes

        data = np.array(df).astype(np.float)
        data = RobustScaler().fit_transform(data)

        new_data, orig_data, new_classes, orig_classes = train_test_split(data, classes, test_size=0.3)

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

        if(t == "CNV"):
            targ_labels = ["GAIN","NEUT","LOSS"]
            targ_dict = {'NEUT': 0, 'LOSS': 0, 'GAIN': 0}
        else:
            targ_labels = ["UPREG","NEUTRAL","DOWNREG"]
            targ_dict = {'NEUTRAL': 0, 'DOWNREG': 0, 'UPREG': 0}

        df1 = df1.reset_index()
        df1["Gene"], df1["Cell Line"] = df1["index"].str.split("_", 1).str
        one_gene_df = df1.drop(columns=["index", "Cell Line"]).groupby(["Gene", t]).median().reset_index().set_index("Gene")
        one_gene_class = pd.DataFrame(one_gene_df[t])
        one_gene_class = one_gene_class.reset_index()

        up_df = one_gene_df.loc[(one_gene_df[t] == targ_labels[0])]
        neut_df = one_gene_df.loc[(one_gene_df[t] == targ_labels[1])]
        down_df = one_gene_df.loc[(one_gene_df[t] == targ_labels[2])]

        up_genes = up_df.index.values.tolist()
        neut_genes = neut_df.index.values.tolist()
        down_genes = down_df.index.values.tolist()

        _ = one_gene_df.pop(t)

        column_squigly = {}
        for col in one_gene_df.columns:
            v1 = up_df[col].median()
            v2 = neut_df[col].median()
            v3 = down_df[col].median()

            correl = np.corrcoef([v1,v2,v3],[1.0,0.0,-1.0])

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
            sorted_d = sorted(temp_dict_feat.items(), key=operator.itemgetter(1), reverse=True)
            return sorted_d
        sorted_d = idx_change(header, rfc.feature_importances_)

        feat = []
        gini = []
        corr = []

        x=0
        while(x<10): # Get the first 10 features
            tempa = sorted_d[x]
            feat.append(tempa[0])
            gini.append(tempa[1])
            corr.append(str(column_squigly[tempa[0]]))
            x = x+1

        importance = pd.DataFrame({"Feature":feat, "Gini":gini, "R":corr})

        # Map to label
        if t == 'CNV':
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
        up = up.reset_index().rename(columns={'index':"feature"})
        up = pd.melt(up, id_vars=["feature"])
        up["type"] = class_col[0]

        neut = tmparr_neut[features].T
        neut = neut.reset_index().rename(columns={'index':"feature"})
        neut = pd.melt(neut, id_vars=["feature"])
        neut["type"] = class_col[1]

        down = tmparr_down[features].T
        down = down.reset_index().rename(columns={'index':"feature"})
        down = pd.melt(down, id_vars=["feature"])
        down["type"] = class_col[2]

        final_df = pd.concat([up, neut, down], axis=0)
        final_df['Cancer'] = canc
        final_df = final_df.reset_index().drop('index', axis=1)

        if t == "TCGA annotation":
            t = "TCGA"

        final_df["Target"] = t
        all_dfs.append(final_df)

big_df = pd.concat(all_dfs, axis=0, ignore_index=True)
#big_df.to_csv("metoncofit.csv")
big_df.to_json("metoncofit.json")
