#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
add_kcat.py turns the 0 values in the cancer models to the median kcat value.
    * As a result, some of the top 10 most important predictions were questionable
    * Now, I will add the option of improementing the mean, mode, min, and maximum as well
"""
import pandas as pd
import numpy as np
import os, sys, re

def add_kcat(datapath='str', savepath='str', type='str'):
    """
    The type argument will be either the mean, median, mode, min, or maximum kcat value. The default argument for the type will be 'median'.
    """

    # Read in all files
    if datapath is None:
        datapath = './../data/'
    if savepath is None:
        savepath = './../data/'

    for fil in os.listdir(datapath):
        filname = fil.split('.')[0]

        # Read in file and open as dataframe
        df = pd.read_csv(datapath+fil)
        df['kcat'] = df['kcat'].replace(0, np.NaN)

        # Accessory functions
        if type == 'median':
            median_kcat = df['kcat'].median()
            df['kcat'] = df['kcat'].replace(np.NaN, median_kcat)
        elif type == 'mean':
            mean_kcat = df['kcat'].mean()
            df['kcat'] = df['kcat'].replace(np.NaN, mean_kcat)
        elif type == 'mode':
            mode_kcat = df['kcat'].mode('columns')
            df['kcat'] = df['kcat'].replace(np.NaN, [mode_kcat])
        elif type == 'max':
            max_kcat = df['kcat'].max()
            df['kcat'] = df['kcat'].replace(np.NaN, max_kcat)
        elif type == 'min':
            min_kcat = df['kcat'].min()
            df['kcat'] = df['kcat'].replace(np.NaN, min_kcat)

        # Save dataframe as file. Will save as same filename unless extension is specified
        df.to_csv(savepath+filname+'.csv', index=False)
        #return df

def remove_kcat(datapath='str', savepath='str', ext='str'):
    """
    In the event that this procedure fails, I wrote a way to reverse it. Since the median value will essentially be the mode, this function finds the model and replaces it with a np.NaN. The file can be fed back into the add_kcat function above.
    """

    # Read in all files
    if datapath is None:
        datapath = './../data/'
    if savepath is None:
        savepath = './../data/'

    for fil in os.listdir(datapath):
        filname = fil.split('.')[0]

        df = pd.read_csv(datapath+fil)

        # Function to remove most abundant kcat value in column
        val = df['kcat'].mode('columns')
        df['kcat'] = df['kcat'].replace([val], 0)

        # Save dataframe as original file name.
        df.to_csv(savepath+filname+'.csv', index=False)

# Remove crap:
#remove_kcat(datapath='./../data/median/', savepath='./../data/original/', ext='')

# Get Mean
#add_kcat(datapath='./../data/okdev/', savepath='./../data/mean/', type='mean')

# Get Median
#add_kcat(datapath='./../data/okdev/', savepath='./../data/median/', type='median')

# Get Mode
add_kcat(datapath='./../data/okdev/', savepath='./../data/mode/', type='mode')

# Get Max
#add_kcat(datapath='./../data/okdev/', savepath='./../data/max/', type='max')

# Get Min
#add_kcat(datapath='./../data/okdev/', savepath='./../data/min/', type='min')
