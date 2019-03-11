# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 11:56:57 2018
Filter all csv files for duplicate entries
@author: scampit
"""

import pandas as pd
import os

for fil in os.listdir():
    if fil.endswith(".csv"):
        nam = fil.split('.')[0].str
        df = pd.read_csv(fil)
        df = df[~df.index.duplicated(keep='first')]
        df.to_csv(nam+'train.csv')
