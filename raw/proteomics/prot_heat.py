import sys
import itertools
import numpy as np
from scipy.stats import ks_2samp
from matplotlib import pyplot
import seaborn as sns
import statsmodels
import pandas as pd


sns.set_context("talk")
sns.set_style("darkgrid", {"axes.facecolor": ".9"})

df = pd.read_csv(sys.argv[1],sep="\t",engine="c",index_col=0)
#df = np.log10(df)
#ind_df = df.set_index('Gene')
xlab = df.columns.tolist()
ylab = df.T.columns.tolist()
#sns.heatmap(data=df, vmax=3.0,vmin=-3.0,cmap="RdBu",yticklabels=ylab,xticklabels=xlab) 
sns.heatmap(data=df,vmax=9.0,cmap="RdBu",yticklabels=ylab,xticklabels=xlab) 
pyplot.yticks(rotation=0) 
pyplot.xticks(rotation=90) 
pyplot.title(sys.argv[1])
pyplot.show()
pyplot.close()
#df.to_csv("log_scale.tsv",sep="\t")
