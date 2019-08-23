#!/usr/bin/env python3


import sys
import pandas as pd
import csv
from pandas import DataFrame


df_ko = pd.read_csv("./Cobra/Gene KO Tables/Melanoma Gene Knockout Table.csv")
df_breast = pd.read_csv("../data/median/melanoma.csv")
#print(list(df_ko.columns))

df_new=df_ko[['grRatio','Cell Line','Gene']]
#print(df_new)



df_final = pd.merge(df_breast, df_new, on=['Gene','Cell Line'])
export_csv = df_final.to_csv (r'/Users/kirksmith/Documents/GitHub.nosync/MetOncoFit/data/geneko/melanoma.csv', index = None, header=True) #Don't forget to add '.csv' at the end of the path
