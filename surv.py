"""
New survival analysis
@author: Scott Campit
"""

import pandas as pd
import numpy as np

# Filters
remove_col = ["TYPE", "ID_DESCRIPTION", "DATA_POSTPROCESSING", "DATASET", "SUBTYPE", "ENDPOINT", "COHORT", "CONTRIBUTOR", "PROBE ID", "ARRAY TYPE", "N", "CUTPOINT", "MINIMUM P-VALUE", "CORRECTED P-VALUE", "ln(HR-high / HR-low)", "ln(HR)"]
cancers = ["Breast cancer", "Ovarian cancer", "Colorectal cancer", "Lung cancer", "Prostate cancer", "Skin cancer", "Brain cancer", "Renal cell carcinoma", "Blood cancer"]

# Process data and only get the COX P-value and Hazard ratio
df = pd.read_excel("./raw/prognoscan/prognoscan.xlsx")
df = df.drop(columns=remove_col, axis=1)
df["HR [95% CI-low CI-upp]"] = df["HR [95% CI-low CI-upp]"].str.replace('\[(.*?)\]', '', regex=True)
df["CANCER TYPE"] = df["CANCER TYPE"].str.replace(' cancer', '', regex=True)
df["CANCER TYPE"] = df["CANCER TYPE"].str.lower()+
df = df[df["CANCER TYPE"].isin(cancers)]
df["HR [95% CI-low CI-upp]"] = df["HR [95% CI-low CI-upp]"].apply(pd.to_numeric)
df["SURV"] = ""

# Make the actual labels
df["SURV"].loc[df["HR [95% CI-low CI-upp]"] >= 2] = "UPREG"
df["SURV"].loc[df["HR [95% CI-low CI-upp]"] <= 0.5] = "DOWNREG"
df["SURV"].loc[(df["HR [95% CI-low CI-upp]"] < 2) & (df["HR [95% CI-low CI-upp]"] > 0.5)] = "NEUTRAL"

# Majority vote on the labels if there are multiple genes and they each have different labels
df = df.groupby(['ID_NAME', 'CANCER TYPE'])['SURV'].agg(pd.Series.mode).to_frame()
df = df.reset_index()

# I physically edited the csv file
#df.to_csv('surv.csv')

df = pd.read_csv('surv.csv')
print(df)

# Read in the existing model and format it for our analysis
fil = r"./data/new/breast.csv"
model = pd.read_csv(fil)
model["Gene"], model["Cell Line"] = model["GENE"].str.split('_', 1).str
model = model.drop(columns='GENE', axis=1)

# Drop existing labels
model = model.drop(columns="SURV", axis=1)
