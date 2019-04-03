import pandas as pd

df = pd.read_csv("nci60_proteome_LFQ.csv",sep="\t",index_col=0)

for col in df.columns.tolist():
    df[col].to_csv(col+".protexp",sep="\t")
