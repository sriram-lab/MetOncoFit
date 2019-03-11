#!/bin/sh
### Main Figures 2-5 (With additional supplementary figures --> the first three for every cancer)

# Breast Cancer
#python metoncofit.py breast.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py breast.train.csv SURV var_excl
#python metoncofit.py breast.train.csv TCGA_annot var_excl
#python metoncofit.py breast.train.csv CNV var_excl
#python metoncofit.py breast.train.csv SURV_wGeneCNV

# NSCLC
#python metoncofit.py nsclc.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py nsclc.train.csv SURV var_excl
#python metoncofit.py nsclc.train.csv TCGA_annot var_excl
#python metoncofit.py nsclc.train.csv CNV var_excl
#python metoncofit.py nsclc.train.csv SURV_wGeneCNV

# Melanoma
#python metoncofit.py melanoma.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py melanoma.train.csv SURV var_excl
#python metoncofit.py melanoma.train.csv TCGA_annot var_excl
#python metoncofit.py melanoma.train.csv CNV var_excl
#python metoncofit.py melanoma.train.csv SURV_wGeneCNV

# Pan Cancer
#python metoncofit.py complex.train.csv CNV var_excl
#python metoncofit.py complex.train.csv SURV var_excl
#python metoncofit.py complex.train.csv TCGA_annot var_excl
#python metoncofit.py complex.train.csv CNV_wGene var_excl_cnv
#python metoncofit.py complex.train.csv SURV_wGeneCNV


### Other Cancer Models
# Leukemia
#python metoncofit.py leukemia.train.csv CNV_wGene var_excl_cnv
python metoncofit.py leukemia.train.csv SURV var_excl
python metoncofit.py leukemia.train.csv TCGA_annot var_excl
python metoncofit.py leukemia.train.csv CNV var_excl
#python metoncofit.py leukemia.train.csv SURV_wGeneCNV

# Colon Cancer
#python metoncofit.py colon.train.csv CNV_wGene var_excl_cnv
python metoncofit.py colon.train.csv SURV var_excl
python metoncofit.py colon.train.csv TCGA_annot var_excl
python metoncofit.py colon.train.csv CNV var_excl
#python metoncofit.py colon.train.csv SURV_wGeneCNV

# Prostate Cancer
#python metoncofit.py prostate.train.csv CNV_wGene var_excl_cnv
python metoncofit.py prostate.train.csv SURV var_excl
python metoncofit.py prostate.train.csv TCGA_annot var_excl
python metoncofit.py prostate.train.csv CNV var_excl
#python metoncofit.py prostate.train.csv SURV_wGeneCNV

# Ovarian Cancer
#python metoncofit.py ovarian.train.csv CNV_wGene var_excl_cnv
python metoncofit.py ovarian.train.csv SURV var_excl
python metoncofit.py ovarian.train.csv TCGA_annot var_excl
python metoncofit.py ovarian.train.csv CNV var_excl
#python metoncofit.py ovarian.train.csv SURV_wGeneCNV

# Central Nervous System Cancer
#python metoncofit.py cns.train.csv CNV_wGene var_excl_cnv
python metoncofit.py cns.train.csv SURV var_excl
python metoncofit.py cns.train.csv TCGA_annot var_excl
python metoncofit.py cns.train.csv CNV var_excl
#python metoncofit.py cns.train.csv SURV_wGeneCNV

# Renal Cancer
#python metoncofit.py renal.train.csv CNV_wGene var_excl_cnv
python metoncofit.py renal.train.csv SURV var_excl
python metoncofit.py renal.train.csv TCGA_annot var_excl
python metoncofit.py renal.train.csv CNV var_excl
#python metoncofit.py renal.train.csv SURV_wGeneCNV
