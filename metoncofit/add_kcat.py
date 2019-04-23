#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os, sys

datapath = ''
for fil in os.listdir('./../data'):
    if datapath is None:
        datapath = './../data/'
    df = pd.read_csv(datapath+fil)
    print(df)
    df = df.replace(0, np.NaN)
    median_kcat = df['kcat'].median()
    df = df.replace(np.NaN, median_kcat)
    print(df['kcat'].head(100))
