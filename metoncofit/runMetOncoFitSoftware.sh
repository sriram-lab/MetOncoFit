#!/bin/sh

# Differential expression prediction
#python3 metoncofit.py './../data/median/breast.csv' TCGA_annot './../labels/var_excl'
#python3 metoncofit.py nsclc.csv TCGA_annot var_excl
#python3 metoncofit.py melanoma.csv TCGA_annot var_excl

# Copy number variation prediction
#python3 metoncofit.py breast.csv CNV var_excl
#python3 metoncofit.py nsclc.csv CNV var_excl
#python3 metoncofit.py melanoma.csv CNV var_excl

# Prognoscan patient survival prediction
#python3 metoncofit.py breast.csv SURV var_excl
#python3 metoncofit.py nsclc.csv SURV var_excl
#python3 metoncofit.py melanoma.csv SURV var_excl

# Pan cancer prediction
#python3 metoncofit.py complex.csv CNV var_excl
#python3 metoncofit.py complex.csv SURV var_excl
#python3 metoncofit.py complex.csv TCGA_annot var_excl
#
# Predicting copy number variation using all data values (TCGA gene expression fold change and CNV ratio)
# python3 metoncofit.py breast.csv CNV var_excl_cnv
# python3 metoncofit.py nsclc.csv CNV var_excl_cnv
# python3 metoncofit.py melanoma.csv CNV var_excl_cnv
#
# Predicting patient survival including all values (TCGA gene expression fold change and CNV ratio)
# python3 metoncofit.py breast.csv SURV no_excl
# python3 metoncofit.py nsclc.csv SURV no_excl
# python3 metoncofit.py melanoma.csv SURV no_excl

# Leukemia cancer model
#python3 metoncofit.py leukemia.csv SURV var_excl
#python3 metoncofit.py leukemia.csv TCGA_annot var_excl
#python3 metoncofit.py leukemia.csv CNV var_excl

# Colon cancer model
#python3 metoncofit.py colon.csv SURV var_excl
#python3 metoncofit.py colon.csv TCGA_annot var_excl
#python3 metoncofit.py colon.csv CNV var_excl

# Prostate cancer model
#python3 metoncofit.py prostate.csv SURV var_excl
#python3 metoncofit.py prostate.csv TCGA_annot var_excl
#python3 metoncofit.py prostate.csv CNV var_excl

# Ovarian cancer model
#python3 metoncofit.py ovarian.csv SURV var_excl
#python3 metoncofit.py ovarian.csv TCGA_annot var_excl
#python3 metoncofit.py ovarian.csv CNV var_excl

# Central Nervous System cancer model
#python3 metoncofit.py cns.csv SURV var_excl
#python3 metoncofit.py cns.csv TCGA_annot var_excl
#python3 metoncofit.py cns.csv CNV var_excl

# Renal cancer model
#python3 metoncofit.py renal.csv SURV var_excl
#python3 metoncofit.py renal.csv TCGA_annot var_excl
#python3 metoncofit.py renal.csv CNV var_excl
