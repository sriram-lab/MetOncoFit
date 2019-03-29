#!/bin/sh
################################################################################
### Main Figures 2-5 (With additional supplementary figures --> the first three for every cancer)

# Figure 2: Differential expression
python3 metoncofit.py breast.train.csv TCGA_annot var_excl
python3 metoncofit.py nsclc.train.csv TCGA_annot var_excl
python3 metoncofit.py melanoma.train.csv TCGA_annot var_excl

# Figure 3: Predicing copy number variation
python3 metoncofit.py breast.train.csv CNV var_excl
python3 metoncofit.py nsclc.train.csv CNV var_excl
python3 metoncofit.py melanoma.train.csv CNV var_excl

# Figure 4: Predicting cancer patient survival
python3 metoncofit.py breast.train.csv SURV var_excl
python3 metoncofit.py nsclc.train.csv SURV var_excl
python3 metoncofit.py melanoma.train.csv SURV var_excl

# Figure 5: Pan cancer predictions
python3 metoncofit.py complex.train.csv CNV var_excl
python3 metoncofit.py complex.train.csv SURV var_excl
python3 metoncofit.py complex.train.csv TCGA_annot var_excl

################################################################################
### S. Figures 1 and 8: Including more expression features for differential expression and cancer patient survival

# S. Figure X: Predicting copy number variation using values
#python3 metoncofit.py breast.train.csv TCGA_annot var_excl_cnv
#python3 metoncofit.py nsclc.train.csv TCGA_annot var_excl_cnv
#python3 metoncofit.py melanoma.train.csv TCGA_annot var_excl_cnv

# S. Figure 1: Predicting copy number variation using values
#python3 metoncofit.py breast.train.csv CNV var_excl_cnv
#python3 metoncofit.py nsclc.train.csv CNV var_excl_cnv
#python3 metoncofit.py melanoma.train.csv CNV var_excl_cnv

# S. Figure 8: Predicting patient survival including all values (TCGA gene expression fold change and CNV ratio)
#python3 metoncofit.py breast.train.csv SURV no_excl
#python3 metoncofit.py nsclc.train.csv SURV no_excl
#python3 metoncofit.py melanoma.train.csv SURV no_excl

################################################################################
### S. Figures 2-7: The other cancer models not described in the main text

# Leukemia
python3 metoncofit.py leukemia.train.csv SURV var_excl
python3 metoncofit.py leukemia.train.csv TCGA_annot var_excl
python3 metoncofit.py leukemia.train.csv CNV var_excl

# Colon Cancer
python3 metoncofit.py colon.train.csv SURV var_excl
python3 metoncofit.py colon.train.csv TCGA_annot var_excl
python3 metoncofit.py colon.train.csv CNV var_excl

# Prostate Cancer
python3 metoncofit.py prostate.train.csv SURV var_excl
python3 metoncofit.py prostate.train.csv TCGA_annot var_excl
python3 metoncofit.py prostate.train.csv CNV var_excl

# Ovarian Cancer
python3 metoncofit.py ovarian.train.csv SURV var_excl
python3 metoncofit.py ovarian.train.csv TCGA_annot var_excl
python3 metoncofit.py ovarian.train.csv CNV var_excl

# Central Nervous System Cancer
python3 metoncofit.py cns.train.csv SURV var_excl
python3 metoncofit.py cns.train.csv TCGA_annot var_excl
python3 metoncofit.py cns.train.csv CNV var_excl

# Renal Cancer
python3 metoncofit.py renal.train.csv SURV var_excl
python3 metoncofit.py renal.train.csv TCGA_annot var_excl
python3 metoncofit.py renal.train.csv CNV var_excl
