"""
Data cleanup function
@author: Scott Campit
"""
import os
import pandas as pd
import numpy as np

# Clean up original dataset
path = r"./../data/okdev/"
folder = os.listdir(path)

def clean(path, fil, canc):

    if path is None:
        path = r"./../data/okdev/"

    # Read in the existing model and format it for our analysis
    model = pd.read_csv(path+fil)

    # Main cleaning functions:
    model = model.drop_duplicates(subset='GENE', keep='first')
    model["Gene"], model["Cell Line"] = model["GENE"].str.split('_', 1).str
    model = model.drop(columns='GENE', axis=1)
    model = model.set_index(['Gene', 'Cell Line'])
    model = model.reset_index()
    print(canc+' is done!')
    return model

for fil in folder:
    canc, _ = os.path.splitext(fil)
    canc = canc.split('.')[0]
    model = clean(path, fil, canc)
    model.to_csv(path+canc+'.csv', index=False)
