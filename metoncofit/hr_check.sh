#!/bin/sh

# Relaxed model

## Differential expression
python3 hr_check.py ./../data/lax/ breast.csv TCGA_annot var_excl
python3 hr_check.py ./../data/lax/ nsclc.csv TCGA_annot var_excl
python3 hr_check.py ./../data/lax/ melanoma.csv TCGA_annot var_excl

## Copy number variation
python3 hr_check.py ./../data/lax/ breast.csv CNV var_excl
python3 hr_check.py ./../data/lax/ nsclc.csv CNV var_excl
python3 hr_check.py ./../data/lax/ melanoma.csv CNV var_excl

## Cancer patient survival
python3 hr_check.py ./../data/lax/ breast.csv SURV var_excl
python3 hr_check.py ./../data/lax/ nsclc.csv SURV var_excl
python3 hr_check.py ./../data/lax/ melanoma.csv SURV var_excl

# Normal model
python3 hr_check.py ./../data/median/ breast.csv TCGA_annot var_excl
python3 hr_check.py ./../data/median/ nsclc.csv TCGA_annot var_excl
python3 hr_check.py ./../data/median/ melanoma.csv TCGA_annot var_excl

## Copy number variation
python3 hr_check.py ./../data/median/ breast.csv CNV var_excl
python3 hr_check.py ./../data/median/ nsclc.csv CNV var_excl
python3 hr_check.py ./../data/median/ melanoma.csv CNV var_excl

## Cancer patient survival
python3 hr_check.py ./../data/median/ breast.csv SURV var_excl
python3 hr_check.py ./../data/median/ nsclc.csv SURV var_excl
python3 hr_check.py ./../data/median/ melanoma.csv SURV var_excl

# Stringent model
python3 hr_check.py ./../data/stringent/ breast.csv TCGA_annot var_excl
python3 hr_check.py ./../data/stringent/ nsclc.csv TCGA_annot var_excl
python3 hr_check.py ./../data/stringent/ melanoma.csv TCGA_annot var_excl

## Copy number variation
python3 hr_check.py ./../data/stringent/ breast.csv CNV var_excl
python3 hr_check.py ./../data/stringent/ nsclc.csv CNV var_excl
python3 hr_check.py ./../data/stringent/ melanoma.csv CNV var_excl

## Cancer patient survival
python3 hr_check.py ./../data/stringent/ breast.csv SURV var_excl
python3 hr_check.py ./../data/stringent/ nsclc.csv SURV var_excl
python3 hr_check.py ./../data/stringent/ melanoma.csv SURV var_excl
