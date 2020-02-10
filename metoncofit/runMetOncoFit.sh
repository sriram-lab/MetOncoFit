#!/bin/sh
################################################################################
### Main Figures 2-5 (With additional supplementary figures --> the first three for every cancer)
path="~/Data/MetOncoFit/" # path where data is
file="var_excl"           # file type for exclusion
fullpath="${path}${file}"

# Figure 2: Differential expression
python3 metoncofit.py breast.csv TCGA_annot fullpath
python3 metoncofit.py nsclc.csv TCGA_annot fullpath
python3 metoncofit.py melanoma.csv TCGA_annot fullpath

# Figure 3: Predicing copy number variation
python3 metoncofit.py breast.csv CNV fullpath
python3 metoncofit.py nsclc.csv CNV fullpath
python3 metoncofit.py melanoma.csv CNV fullpath

# Figure 4: Predicting cancer patient survival
python3 metoncofit.py breast.csv SURV fullpath
python3 metoncofit.py nsclc.csv SURV fullpath
python3 metoncofit.py melanoma.csv SURV fullpath

# Figure 5: Pan cancer predictions
python3 metoncofit.py complex.csv CNV fullpath
python3 metoncofit.py complex.csv SURV fullpath
python3 metoncofit.py complex.csv TCGA_annot fullpath

