import pandas as pd
import numpy as np
import sys

#fil=open("good_present.txt")
#cols = [i.strip() for i in fil.readlines()]

#df = pd.read_csv("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct",sep="\t",index_col=1)
#df[cols].to_csv("small_gtex.tsv",sep="\t")

df = pd.read_csv(sys.argv[1],sep="\t",index_col=0)
df['Median'] = df.median(axis=1)
df['Mean'] = df.mean(axis=1)
cols = ['Median','Mean']
df[cols].to_csv("proteome_mean_median.tsv",sep="\t")

