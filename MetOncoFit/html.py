# -*- coding: utf-8 -*-
"""
html_build contains the functions that outputs the data in html-friendly format required to save the svg files.

@author: Scott Campit
"""
import pandas as pd
import scipy

def header():
    """
    The header file contains information that will construct the html file used to obtain the svg figure in the end.
    """
    fil=open("./../html_builder/header")
    for lin in fil.readlines():
        return(lin.strip())
    fil.close()

def tailer():
    """
    The tailer file contains useful parameters that can be adjusted to change the heatmap
    """
    fil=open("../html_builder/tailer_heatmap")
    for lin in fil.readlines():
        return(lin.strip())
    fil.close()

def model_validation_scores():
    """
    """
    return("var result_vars = [")
    return("{\"naam\":\"Cross validation accuracy (%)\",\"value\":"+str(acc)+",\"value2\":"+str(stddev)+"},")
    return("{\"naam\":\"p-value of accuracy\",\"value\":"+str(new_pvalue)+"},")
    return("{\"naam\":\"Z score of accuracy\",\"value\":"+str(zscore)+"},")
    return("];")

def one_gene_only(df, target):
    """
    one_gene_only will merge the gene target value by majority rules and will take the median values for all numerical values.
    """

    global up_df, neut_df, down_df, up_genes, neut_genes, down_genes

    if(target == "CNV"):
        targ_labels = ["GAIN","NEUT","LOSS"]
        targ_dict = {'NEUT': 0, 'LOSS': 0, 'GAIN': 0}
    else:
        targ_labels = ["UPREG","NEUTRAL","DOWNREG"]
        targ_dict = {'NEUTRAL': 0, 'DOWNREG': 0, 'UPREG': 0}

    df = df.reset_index()
    df["Gene"], df["Cell Line"] = df["index"].str.split("_", 1).str
    one_gene_df = df.drop(columns=["index", "Cell Line"]).groupby(["Gene", target]).median().reset_index().set_index("Gene")
    one_gene_class = pd.DataFrame(one_gene_df[target])
    one_gene_class = one_gene_class.reset_index()

    # These dataframes contain the df entries with increased, neutral, and decreased values.
    up_df = one_gene_df.loc[(one_gene_df[target] == targ_labels[0])]
    neut_df = one_gene_df.loc[(one_gene_df[target] == targ_labels[1])]
    down_df = one_gene_df.loc[(one_gene_df[target] == targ_labels[2])]

    # To create the figure, we are randomly selecting three genes that are upreg, neutral, or downreg and are storing them in this list.
    up_genes = up_df.index.values.tolist()
    neut_genes = neut_df.index.values.tolist()
    down_genes = down_df.index.values.tolist()

    return up_df, neut_df, down_df, up_genes, neut_genes, down_genes, one_gene_df, one_gene_class

def correlation(up_df, neut_df, down_df):
    """
    """

    # This will calculate the correlation, if there is one between the biological features.
    column_squigly = {}
    for col in df.columns:
        v1 = up_df[col].median()
        v2 = neut_df[col].median()
        v3 = down_df[col].median()

        correl = np.corrcoef([v1,v2,v3],[1.0,0.0,-1.0])

        if(np.isnan(correl[0][1]) != True):
            column_squigly[col] = correl[0][1]
        else:
            column_squigly[col] = 0.0

    return(column_squigly)
"""
def build_heatmap(one_gene_df, up_df, neut_df, down_df, size=3, clf, norm='standard'):

    build_heatmap will construct the table that will be used to visualize the heatmap.

    INPUTS:
        `size` refers to the number of genes that will be outputted in the heatmap.
        The default size is 3.

        `feat_imp` is the feature importances that will be used on the y axis of the heatmap.

        `fir_df, sec_df, trd_df` refer to the dataframes specific to the three individual classes. These are made in the `correlation` function.

    OUTPUTS:
        feat_imp table
        heatmap_arr table
        gene_names associated with each column of the heatmap


    # To create the figure, we are randomly selecting three genes that are upreg, neutral, or downreg and are storing them in this list.
    example_genes = []

    example_genes.append(up_df.sample(n=3).index.tolist())
    example_genes.append(neut_df.sample(n=3).index.tolist())
    example_genes.append(down_df.sample(n=3).index.tolist())

    # This code will output the list of feature importances that will be incorporated into the heatmap with the associated example genes.
    temp_dict_feat = {}
    for i, j in zip(header, clf.feature_importances_):
        temp_dict_feat[i] = j
    sorted_d = sorted(temp_dict_feat.items(),key=operator.itemgetter(1),reverse=True)
    x=0
    ex_heat = []

    print("var feat_imp=[")
    while(x<70):
        tempa = sorted_d[x]
        print("\t{\"xval\":"+str(x)+",\"yval\":"+str(tempa[1])+",\"name\":\""+str(real_names[tempa[0]])+"\",\"correl\":"+str(column_squigly[tempa[0]])+"},")
        absmax = one_gene_df[tempa[0]].quantile(0.33)
        real_min = one_gene_df[tempa[0]].min()
        absmin  = real_min
        real_max = one_gene_df[tempa[0]].max()
        real_diff = real_max - real_min
        tmparr = []
        # First set of genes
        for i in range(0,3):
            tmparr.append((one_gene_df[tempa[0]][example_genes[0][i]]-real_min)/real_diff)
        # Second set of genes
        for j in range(0,3):
            tmparr.append((one_gene_df[tempa[0]][example_genes[1][j]]-real_min)/real_diff)
        # Third set of genes
        for k in range(0,3):
            tmparr.append((one_gene_df[tempa[0]][example_genes[2][k]]-real_min)/real_diff)
        ex_heat.append(tmparr)
        x=x+1
    print("];")

    # This is the code to actually construct the heatmap array that will be outputted.
    print("var heatmap_arr = [")
    x=0
    for row in ex_heat:
        y=0
        offset = 0.0
        for col in row:
            print("{\"xval\":"+str(x)+",\"yval\":"+str(float(y) + offset)+",\"value\":"+str(col)+"},")
            if(y == 2):
                offset = 0.5
            elif(y == 5):
                offset = 1.0
            y=y+1
        x=x+1
    print("];")

    # This will generate the upregulated, neutral, or downregulated genes associated with each biological feature in the figure.
    x=0
    print("var gene_names = [")
    for i in range (0,3):
        print("{\"name\":\""+example_genes[0][i].split("_")[0]+" ("+targ_labels[0]+")\",\"xval\":"+str(x)+"},")
        x=x+1
    for j in range(0,3):
        print("{\"name\":\""+example_genes[1][j].split("_")[0]+" ("+targ_labels[1]+")\",\"xval\":"+str(float(x)+0.5)+"},")
        x=x+1
    for k in range(0,3):
        print("{\"name\":\""+example_genes[2][k].split("_")[0]+" ("+targ_labels[2]+")\",\"xval\":"+str(float(x)+1.0)+"},")
        x=x+1
    print("];")
"""
def build_confusion():
    """
    """
    # The labels that are used (increased, neutral, decreased)
    return("var labelsData = ["+",".join(["\""+i+"\"" for i in targ_labels])+"];")

    # The confusion matrix
    return("var conf_matr = ")
    return(conf_matr.tolist())
    return(";")

    class_max = str(np.amax(conf_matr))

def color_map():
    """
    """
    return("var colorMap = d3.scaleLinear()\n\t.domain([0.0,"+class_max+"])\n\t.range([\"#f7fbff\", \"#08306b\"]);\n")
