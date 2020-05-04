# -*- coding: utf-8 -*-
"""
addFeatures.py
This script contains accessory functions that help expand the feature set for the MetOncoFit databases.

@author: Scott Campit
"""

import sys
import os
import re
import difflib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def explode(input, col=None, sep=None):
    """
    explode separates elements in a single pandas dataframe cell that are separated by a single delimiter or a list
    of delimiters to have its own unique row entry, while duplicating the cell contents.

    Note that this function has the potential to rapidly expand the dataset.

    :param input: A pandas dataframe containing the dataset.
    :param col:   A string denoting the column to search for a specific separator.
    :param regex: A string denoting the specific regex to search for.

    :return df:   A pandas dataframe containing the expanded dataset.
    """
    if type(sep) is str:
        s = input[col].str.split(sep).apply(pd.Series, 1).stack()
        s.index = s.index.droplevel(-1)
        s.name = col
        del input[col]
        input[col] = s
        df = input.drop_duplicates(keep='first')
    elif type(sep) is list:
        for character in sep:
            s = input[col].str.split(character).apply(pd.Series, 1).stack()
            s.index = s.index.droplevel(-1)
            s.name = col
            del input[col]
            input[col] = s
            df = input.drop_duplicates(keep='first')
    return df

def concatFeatures(originalData, newData):
    """
    concatFeatures will take the features from the original dataset (pandas dataframe) and the new dataset (pandas
    dataframe), and concat them based on similar indices.

    :param originalData:
    :param newData:
    :return output:
    """

    return pd.merge(originalData, newData, how='inner', left_index=True, right_index=True)



if __name__ == "__main__":
    None
