"""
Misc data clean up
@author: Scott Campit
"""
import os
import pandas as pd
import numpy as np

# Clean up original dataset
path = r"./data/original/"
folder = os.listdir(path)

def clean(path, fil):

    if path is None:
        path = r"./data/original/"

    # Read in the existing model and format it for our analysis
    model = pd.read_csv(path+fil)

    # Main cleaning function
    model["Gene"], model["Cell Line"] = model["GENE"].str.split('_', 1).str
    model = model.drop(columns='GENE', axis=1)
    model = model.set_index(['Gene', 'Cell Line'])
    model = model.reset_index()
    return model

    print(canc+' is done!')

for fil in folder:
    canc, _ = os.path.splitext(fil)
    model = clean(path, fil)
    model.to_csv(path+canc+'.csv', index=False)
