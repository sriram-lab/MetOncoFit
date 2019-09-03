#!/usr/bin/env python3


import sys
import pandas as pd
import csv
import mygene
from pandas import DataFrame

mg = mygene.MyGeneInfo()

genelist = pd.read_csv("CNS_GeneList.csv")
#print(genelist['genes'])
all = genelist['genelistrounded'][:]

#genes = ['6541','6542']

#results = mg.querymany(all, scopes=["entrezgene"], fields=["symbol"], species="human", verbose=False)
#for res in results:
    #sym = res['query']
    #print(sym)
#print(results)


out=mg.getgenes(all, as_dataframe=True)
#print(out)


df = DataFrame(out['symbol'])
export_csv = df.to_csv (r'/Users/kirksmith/Documents/GitHub.nosync/MetOncoFit/CNS_Symbols.csv', index = None, header=True) #Don't forget to add '.csv' at the end of the path
print(df)
