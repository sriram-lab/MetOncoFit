# -*- coding: utf-8 -*-
"""
"""
import os
import numpy as np
import pandas as pd

datapath = None
for fil in os.listdir('./../data'):
    if datapath is None:
        datapath = './../data/'
    canc = fil.replace(".train.csv", "")
    df = pd.read_csv(datapath+fil)
    df['kcat'] = df['kcat'].replace(0, np.NaN)
    print(df['kcat'])
    #median_kcat = df['kcat'].median()
    #df['kcat'] = df['kcat'].replace(np.NaN, median_kcat)
    #df.to_json(datapath+canc+'.json', orient='columns')
