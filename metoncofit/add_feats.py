# -*- coding: utf-8 -*-
"""
Add features

@author: scott
"""

import sys
import os
import re
import difflib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def sep_obj(df, col='', regex=' '):
    """
    sep_obj is used to separate elements in a single cell within a dataframe to have its own unique entry. This will rapidly expand the dataset
    """
    if type(regex) is str:
        s = df[col].str.split(regex).apply(pd.Series, 1).stack()
        s.index = s.index.droplevel(-1)
        s.name = col
        del df[col]
        df[col] = s
        df = df.drop_duplicates(keep='first')
    elif type(regex) is list:
        for character in regex:
            s = df[col].str.split(character).apply(pd.Series, 1).stack()
            s.index = s.index.droplevel(-1)
            s.name = col
            del df[col]
            df[col] = s
            df = df.drop_duplicates(keep='first')
    return df


# Add median metabolomics values from the NCI-60 cancer cell line database
path = ""
fil = ""
df = pd.read_table(path+fil)
