#!/bin/sh
################################################################################
### Main Figures 2-5 (With additional supplementary figures --> the first three for every cancer)

# Figure 2: Differential expression
#python3 metoncofit.py breast.csv TCGA_annot var_excl
#python3 metoncofit.py nsclc.csv TCGA_annot var_excl
#python3 metoncofit.py melanoma.csv TCGA_annot var_excl

# Figure 3: Predicing copy number variation
#python3 metoncofit.py breast.csv CNV var_excl
#python3 metoncofit.py nsclc.csv CNV var_excl
#python3 metoncofit.py melanoma.csv CNV var_excl

# Figure 4: Predicting cancer patient survival
#python3 metoncofit.py breast.csv SURV var_excl
#python3 metoncofit.py nsclc.csv SURV var_excl
#python3 metoncofit.py melanoma.csv SURV var_excl

# Figure 5: Pan cancer predictions
#python3 metoncofit.py complex.csv CNV var_excl
#python3 metoncofit.py complex.csv SURV var_excl
#python3 metoncofit.py complex.csv TCGA_annot var_excl

################################################################################
### S. Figures 1 and 8: Including more expression features for differential expression and cancer patient survival

# S. Figure X: Predicting differential expression using values
#python3 metoncofit.py breast.csv TCGA_annot var_excl_cnv
#python3 metoncofit.py nsclc.csv TCGA_annot var_excl_cnv
#python3 metoncofit.py melanoma.csv TCGA_annot var_excl_cnv

# S. Figure 1: Predicting copy number variation using values
#python3 metoncofit.py breast.csv CNV var_excl_cnv
#python3 metoncofit.py nsclc.csv CNV var_excl_cnv
#python3 metoncofit.py melanoma.csv CNV var_excl_cnv

# S. Figure 8: Predicting patient survival including all values (TCGA gene expression fold change and CNV ratio)
#python3 metoncofit.py breast.csv SURV no_excl
#python3 metoncofit.py nsclc.csv SURV no_excl
#python3 metoncofit.py melanoma.csv SURV no_excl

################################################################################
### S. Figures 2-7: The other cancer models not described in the main text

# Leukemia
python3 metoncofit.py leukemia.csv SURV var_excl
python3 metoncofit.py leukemia.csv TCGA_annot var_excl
python3 metoncofit.py leukemia.csv CNV var_excl

# Colon Cancer
python3 metoncofit.py colon.csv SURV var_excl
python3 metoncofit.py colon.csv TCGA_annot var_excl
python3 metoncofit.py colon.csv CNV var_excl

# Prostate Cancer
python3 metoncofit.py prostate.csv SURV var_excl
python3 metoncofit.py prostate.csv TCGA_annot var_excl
python3 metoncofit.py prostate.csv CNV var_excl

# Ovarian Cancer
python3 metoncofit.py ovarian.csv SURV var_excl
python3 metoncofit.py ovarian.csv TCGA_annot var_excl
python3 metoncofit.py ovarian.csv CNV var_excl

# Central Nervous System Cancer
python3 metoncofit.py cns.csv SURV var_excl
python3 metoncofit.py cns.csv TCGA_annot var_excl
python3 metoncofit.py cns.csv CNV var_excl

# Renal Cancer
python3 metoncofit.py renal.csv SURV var_excl
python3 metoncofit.py renal.csv TCGA_annot var_excl
python3 metoncofit.py renal.csv CNV var_excl
