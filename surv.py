"""
Labeling survival data:

To label the survival data, I will do the following:
  * An HR > 2 and Cox p-value < 0.05 will be UPREG
  * An HR < 0.5 and Cox p-value < 0.05 will be DOWNREG
  * Everything else: NEUTRAL

@author: Scott Campit
"""
import os
import pandas as pd
import numpy as np

def make_surv(input, cox, hr_up, hr_low):
    """
    make_surv will make a csv file containing the annotations by the Cox p-value and Hazard Ratio thresholds specified by the user. This file needs to be manually edited for multiple modes.
    """
    # Filters
    remove_col = ["TYPE", "ID_DESCRIPTION", "DATA_POSTPROCESSING", "DATASET", "SUBTYPE", "ENDPOINT", "COHORT", "CONTRIBUTOR", "PROBE ID", "ARRAY TYPE", "N", "CUTPOINT", "MINIMUM P-VALUE", "CORRECTED P-VALUE", "ln(HR-high / HR-low)", "ln(HR)"]
    cancers = ["Breast cancer", "Ovarian cancer", "Colorectal cancer", "Lung cancer", "Prostate cancer", "Skin cancer", "Brain cancer", "Renal cell carcinoma", "Blood cancer"]

    # Process data and only get the COX P-value and Hazard ratio
    df = pd.read_excel(input)
    df = df.drop(columns=remove_col, axis=1)
    df["HR [95% CI-low CI-upp]"] = df["HR [95% CI-low CI-upp]"].str.replace('\[(.*?)\]', '', regex=True)
    df = df[df["CANCER TYPE"].isin(cancers)]
    df["HR [95% CI-low CI-upp]"] = df["HR [95% CI-low CI-upp]"].apply(pd.to_numeric)
    df["SURV"] = ""

    # Make the actual labels
    df["SURV"].loc[(df["HR [95% CI-low CI-upp]"] >= hr_up) & (df["COX P-VALUE"] <= cox)] = "UPREG"
    df["SURV"].loc[(df["HR [95% CI-low CI-upp]"] <= hr_low) & (df["COX P-VALUE"] <= cox)] = "DOWNREG"
    df["SURV"].loc[(df["HR [95% CI-low CI-upp]"] >= hr_up) & (df["COX P-VALUE"] > cox)] = "NEUTRAL"
    df["SURV"].loc[(df["HR [95% CI-low CI-upp]"] <= hr_low) & (df["COX P-VALUE"] > cox)] = "NEUTRAL"
    df["SURV"].loc[(df["HR [95% CI-low CI-upp]"] < hr_up) & (df["HR [95% CI-low CI-upp]"] > hr_low)] = "NEUTRAL"

    # Majority vote on the labels if there are multiple genes and they each have different labels
    df = df.groupby(['ID_NAME', 'CANCER TYPE'])['SURV'].agg(pd.Series.mode).to_frame()
    df = df.reset_index()

    # I physically edited the xlsx file. Need to devise conditional rule set to automatically determine labels for multiple modes
    df.to_excel('surv.xlsx', index=False)

#make_surv("./raw/prognoscan/prognoscan.xlsx", cox=0.05, hr_up=2, hr_low=0.05)

canc_dict = {
    'Breast cancer':'breast',
    'Brain cancer':'cns',
    'Colorectal cancer':'colon',
    'Blood cancer':'leukemia',
    'Skin cancer':'melanoma',
    'Lung cancer':'nsclc',
    'Ovarian cancer':'ovarian',
    'Prostate cancer':'prostate',
    'Renal cell carcinoma':'renal'
}

df = pd.read_excel('surv.xlsx')
df = df.replace({'CANCER TYPE':canc_dict})

# Read in the existing model and format it for our analysis
fil = r"./data/original/breast.csv"
model = pd.read_csv(fil)
filname = os.path.basename(fil)
canc, _ = os.path.splitext(filname)

# Split gene name and cell line.
model["Gene"], model["Cell Line"] = model["GENE"].str.split('_', 1).str
model = model.drop(columns='GENE', axis=1)

# Drop existing survival labels
model = model.drop(columns="SURV", axis=1)

# Create new survival label as empty and populate later
tmp = df[df["CANCER TYPE"] == canc]
tmp = tmp.drop(columns='CANCER TYPE')

model = pd.merge(model, tmp, how='left', left_on='Gene', right_on='ID_NAME').drop(columns='ID_NAME')
model = model.fillna('NEUTRAL')

model = model.set_index(['Gene', 'Cell Line'])
model = model.reset_index()
