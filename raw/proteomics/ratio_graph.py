import sys
import itertools
import numpy as np
from scipy.stats import ks_2samp
from matplotlib import pyplot
import seaborn as sns
import statsmodels
import pandas as pd

#col_names = ['TotalCNV', 'LossCNV', 'GainCNV', 'RatioGain/Loss']
#col_names = ['HighTotal Substitution - coding silent',	'HighTotal Substitution - Missense',	'HighTotal Substitution - Nonsense',	'MedTotal Substitution - coding silent',	'MedTotal Substitution - Missense',	'MedTotal Substitution - Nonsense',	'LowTotal Substitution - coding silent',	'LowTotal Substitution - Missense',	'LowTotal Substitution - Nonsense', 'Total Missense',	'Total Silent',	'Total Nonsense']
#col_names = ['HighTotal Substitution - coding silent','HighTotal Substitution - Missense','HighTotal Substitution - Nonsense','HighRecur3 Substitution - coding silent','HighRecur3 Substitution - Missense','HighRecur3 Substitution - Nonsense','MedTotal Substitution - coding silent','MedTotal Substitution - Missense','MedTotal Substitution - Nonsense','MedRecur3 Substitution - coding silent','MedRecur3 Substitution - Missense','MedRecur3 Substitution - Nonsense','LowTotal Substitution - coding silent','LowTotal Substitution - Missense','LowTotal Substitution - Nonsense','LowRecur3 Substitution - coding silent','LowRecur3 Substitution - Missense','LowRecur3 Substitution - Nonsense','NoneTotal Substitution - coding silent','NoneTotal Substitution - Missense','NoneTotal Substitution - Nonsense','NoneRecur3 Substitution - coding silent','NoneRecur3 Substitution - Missense','NoneRecur3 Substitution - Nonsense','Total Missense','Total Recur3 Missense','Total Nonsense','Total Recur3 Nonsense','Total Mut','Total Recur3']
col_names = ['Mean']


df = pd.read_table(sys.argv[1],sep="\t",engine="python")
#df2 = df[df.columns.difference(['Grpno',   'Gene',    'HighKcat',    'MedianKcat', 'Ratio Gain/Loss'])]
central_CEM =  df.loc[df['Module'] == 'Primary - Carbohydrate & Energy Metabolism']
central_AFM =  df.loc[df['Module'] == 'Primary - amino acids, fatty acids and nucleotides']
intermediate =  df.loc[df['Module'] == 'Intermediate']
secondary =  df.loc[df['Module'] == 'Secondary']
sns.set_style("dark")
for col in col_names:
    print central_CEM[col].median(axis=1)
    g = sns.distplot(central_CEM[col],kde_kws={'cumulative':True,'lw':3},hist=False,label="Central (carbohydrate)")
    g = sns.distplot(central_AFM[col],kde_kws={'cumulative':True,'lw':3},hist=False,label="Central (amino,fatty,nucl)")
    g = sns.distplot(intermediate[col],kde_kws={'cumulative':True,'lw':3},hist=False,label="Intermediate)")
    g = sns.distplot(secondary[col],kde_kws={'cumulative':True,'lw':3},hist=False,label="Secondary")
    pyplot.title(col)
    #pyplot.show()
    pyplot.savefig(col.replace(" ","_").replace("/","_")+".png",dpi=300)
    pyplot.clf()
    
#g = sns.boxplot(x="c", hue="Type", y="a", data=df_long,fliersize=1.0,orient="h")
#g = sns.boxplot(x="Type", hue="Type", y="c", data=df2,fliersize=0.0)
pyplot.xticks(rotation=45)



