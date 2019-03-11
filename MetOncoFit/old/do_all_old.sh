#!/bin/sh
### Main Figures 2-5 (With additional supplementary figures --> the first three for every cancer)
# Breast Cancer
python metoncofit_expand_heatmap.py breast.train.csv CNV var_excl > ./../output/breast_cnv.html
python metoncofit_expand_heatmap.py breast.train.csv SURV var_excl > ./../output/breast_surv.html
python metoncofit_expand_heatmap.py breast.train.csv TCGA_annot var_excl > ./../output/breast_diff.html
python metoncofit_expand_heatmap.py breast.train.csv CNV var_excl_cnv > ./../output/breast_notcga_cnv.html
python metoncofit_expand_heatmap.py breast.train.csv SURV > ./../output/breast_notcga_surv.html
# NSCLC
python metoncofit_expand_heatmap.py nsclc.train.csv CNV var_excl > ./../output/lung_cnv.html
python metoncofit_expand_heatmap.py nsclc.train.csv SURV var_excl > ./../output/lung_surv.html
python metoncofit_expand_heatmap.py nsclc.train.csv TCGA_annot var_excl > ./../output/lung_diff.html
python metoncofit_expand_heatmap.py nsclc.train.csv CNV var_excl_cnv > ./../output/lung_notcga_cnv.html
python metoncofit_expand_heatmap.py nsclc.train.csv SURV  > ./../output/lung_notcga_surv.html
# Melanoma
python metoncofit_expand_heatmap.py melanoma.train.csv CNV var_excl > ./../output/skin_cnv.html
python metoncofit_expand_heatmap.py melanoma.train.csv SURV var_excl > ./../output/skin_surv.html
python metoncofit_expand_heatmap.py melanoma.train.csv TCGA_annot var_excl > ./../output/skin_diff.html
python metoncofit_expand_heatmap.py melanoma.train.csv CNV var_excl_cnv > ./../output/skin_notcga_cnv.html
python metoncofit_expand_heatmap.py melanoma.train.csv SURV  > ./../output/skin_notcga_surv.html
# Pan Cancer
python metoncofit_expand_heatmap.py complex.train.csv CNV var_excl_cnv > ./../output/complex_notcga_cnv.html
python metoncofit_expand_heatmap.py complex.train.csv SURV  > ./../output/complex_notcga_surv.html
python metoncofit_expand_heatmap.py complex.train.csv CNV var_excl > ./../output/complex_cnv.html
python metoncofit_expand_heatmap.py complex.train.csv SURV var_excl  > ./../output/complex_surv.html
python metoncofit_expand_heatmap.py complex.train.csv TCGA_annot var_excl  > ./../output/complex_diff.html
### Supplementary Figures
# Leukemia
python metoncofit_expand_heatmap.py leukemia.train.csv CNV var_excl_cnv > ./../output/leukemia_notcga_cnv.html
python metoncofit_expand_heatmap.py leukemia.train.csv SURV  > ./../output/leukemia_notcga_surv.html
python metoncofit_expand_heatmap.py leukemia.train.csv CNV var_excl > ./../output/leukemia_cnv.html
python metoncofit_expand_heatmap.py leukemia.train.csv SURV var_excl  > ./../output/leukemia_surv.html
python metoncofit_expand_heatmap.py leukemia.train.csv TCGA_annot var_excl  > ./../output/leukemia_diff.html
# Colon Cancer
python metoncofit_expand_heatmap.py colon.train.csv CNV var_excl_cnv > ./../output/colon_notcga_cnv.html
python metoncofit_expand_heatmap.py colon.train.csv SURV  > ./../output/colon_notcga_surv.html
python metoncofit_expand_heatmap.py colon.train.csv CNV var_excl > ./../output/colon_cnv.html
python metoncofit_expand_heatmap.py colon.train.csv SURV var_excl  > ./../output/colon_surv.html
python metoncofit_expand_heatmap.py colon.train.csv TCGA_annot var_excl  > ./../output/colon_diff.html
# Prostate Cancer
python metoncofit_expand_heatmap.py prostate.train.csv CNV var_excl_cnv > ./../output/prostate_notcga_cnv.html
python metoncofit_expand_heatmap.py prostate.train.csv SURV  > ./../output/prostate_notcga_surv.html
python metoncofit_expand_heatmap.py prostate.train.csv CNV var_excl > ./../output/prostate_cnv.html
python metoncofit_expand_heatmap.py prostate.train.csv SURV var_excl  > ./../output/prostate_surv.html
python metoncofit_expand_heatmap.py prostate.train.csv TCGA_annot var_excl  > ./../output/prostate_diff.html
# Ovarian Cancer
python metoncofit_expand_heatmap.py ovarian.train.csv CNV var_excl_cnv > ./../output/ovarian_notcga_cnv.html
python metoncofit_expand_heatmap.py ovarian.train.csv SURV  > ./../output/ovarian_notcga_surv.html
python metoncofit_expand_heatmap.py ovarian.train.csv CNV var_excl > ./../output/ovarian_cnv.html
python metoncofit_expand_heatmap.py ovarian.train.csv SURV var_excl  > ./../output/ovarian_surv.html
python metoncofit_expand_heatmap.py ovarian.train.csv TCGA_annot var_excl  > ./../output/ovarian_diff.html
# Central Nervous System Cancer
python metoncofit_expand_heatmap.py cns.train.csv CNV var_excl_cnv > ./../output/cns_notcga_cnv.html
python metoncofit_expand_heatmap.py cns.train.csv SURV  > ./../output/cns_notcga_surv.html
python metoncofit_expand_heatmap.py cns.train.csv CNV var_excl > ./../output/cns_cnv.html
python metoncofit_expand_heatmap.py cns.train.csv SURV var_excl  > ./../output/cns_surv.html
python metoncofit_expand_heatmap.py cns.train.csv TCGA_annot var_excl  > ./../output/cns_diff.html
# Renal Cancer
python metoncofit_expand_heatmap.py renal.train.csv CNV var_excl_cnv > ./../output/renal_notcga_cnv.html
python metoncofit_expand_heatmap.py renal.train.csv SURV  > ./../output/renal_notcga_surv.html
python metoncofit_expand_heatmap.py renal.train.csv CNV var_excl > ./../output/renal_cnv.html
python metoncofit_expand_heatmap.py renal.train.csv SURV var_excl  > ./../output/renal_surv.html
python metoncofit_expand_heatmap.py renal.train.csv TCGA_annot var_excl  > ./../output/renal_diff.html
