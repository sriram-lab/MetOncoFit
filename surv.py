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
df = df[df["CANCER TYPE"].isin(cancers)]

df["HR [95% CI-low CI-upp]"] = df["HR [95% CI-low CI-upp]"].apply(pd.to_numeric)
df["SURV"] = ""

df["SURV"].loc[df["HR [95% CI-low CI-upp]"] >= 2] = "UPREG"
df["SURV"].loc[df["HR [95% CI-low CI-upp]"] <= 0.5] = "DOWNREG"
df["SURV"].loc[(df["HR [95% CI-low CI-upp]"] < 2) & (df["HR [95% CI-low CI-upp]"] > 0.5)] = "NEUTRAL"

"""
# Make labels
for _, row in df.iterrows():
    if row["COX P-VALUE"] < 0.05:
        if row["HR [95% CI-low CI-upp]"] >= 2:
            df["SURV"] = "UPREG"
        elif row["HR [95% CI-low CI-upp]"] <= 0.5:
            df["SURV"] = "DOWNREG"
        else:
            df["SURV"] = "NEUTRAL"
    else:
        df["SURV"] = "NEUTRAL"
"""
print(df)
