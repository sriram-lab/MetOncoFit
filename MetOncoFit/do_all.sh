#!/bin/sh
################################################################################
### Main Figures 2-5 (With additional supplementary figures --> the first three for every cancer)

# Figure 2: Differential expression
#python metoncofit.py breast.train.csv TCGA_annot var_excl
#python metoncofit.py nsclc.train.csv TCGA_annot var_excl
#python metoncofit.py melanoma.train.csv TCGA_annot var_excl

# Figure 3: Predicing copy number variation
#python metoncofit.py breast.train.csv CNV var_excl
#python metoncofit.py nsclc.train.csv CNV var_excl
#python metoncofit.py melanoma.train.csv CNV var_excl

# Figure 4: Predicting cancer patient survival
#python metoncofit.py breast.train.csv SURV var_excl
#python metoncofit.py nsclc.train.csv SURV var_excl
#python metoncofit.py melanoma.train.csv SURV var_excl

# Figure 5: Pan cancer predictions
#python metoncofit.py complex.train.csv CNV var_excl
#python metoncofit.py complex.train.csv SURV var_excl
#python metoncofit.py complex.train.csv TCGA_annot var_excl

################################################################################
### S. Figures 1 and 8: Including more expression features for differential expression and cancer patient survival

# S. Figure 1: Predicting copy number variation using values
#python metoncofit.py breast.train.csv CNV var_excl_cnv
#python metoncofit.py nsclc.train.csv CNV var_excl_cnv
#python metoncofit.py melanoma.train.csv CNV var_excl_cnv

# S. Figure 8: Predicting patient survival including all values (TCGA gene expression fold change and CNV ratio)
python metoncofit.py breast.train.csv SURV no_excl
python metoncofit.py nsclc.train.csv SURV no_excl
python metoncofit.py melanoma.train.csv SURV no_excl

################################################################################
### S. Figures 2-7: The other cancer models not described in the main text

# Leukemia
#python metoncofit.py leukemia.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py leukemia.train.csv SURV var_excl
#python metoncofit.py leukemia.train.csv TCGA_annot var_excl
#python metoncofit.py leukemia.train.csv CNV var_excl
#python metoncofit.py leukemia.train.csv SURV_wGeneCNV

# Colon Cancer
#python metoncofit.py colon.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py colon.train.csv SURV var_excl
#python metoncofit.py colon.train.csv TCGA_annot var_excl
#python metoncofit.py colon.train.csv CNV var_excl
#python metoncofit.py colon.train.csv SURV_wGeneCNV

# Prostate Cancer
#python metoncofit.py prostate.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py prostate.train.csv SURV var_excl
#python metoncofit.py prostate.train.csv TCGA_annot var_excl
#python metoncofit.py prostate.train.csv CNV var_excl
#python metoncofit.py prostate.train.csv SURV_wGeneCNV

# Ovarian Cancer
#python metoncofit.py ovarian.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py ovarian.train.csv SURV var_excl
#python metoncofit.py ovarian.train.csv TCGA_annot var_excl
#python metoncofit.py ovarian.train.csv CNV var_excl
#python metoncofit.py ovarian.train.csv SURV_wGeneCNV

# Central Nervous System Cancer
#python metoncofit.py cns.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py cns.train.csv SURV var_excl
#python metoncofit.py cns.train.csv TCGA_annot var_excl
#python metoncofit.py cns.train.csv CNV var_excl
#python metoncofit.py cns.train.csv SURV_wGeneCNV

# Renal Cancer
#python metoncofit.py renal.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py renal.train.csv SURV var_excl
#python metoncofit.py renal.train.csv TCGA_annot var_excl
#python metoncofit.py renal.train.csv CNV var_excl
#python metoncofit.py renal.train.csv SURV_wGeneCNV
