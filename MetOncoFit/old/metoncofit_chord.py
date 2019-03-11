import numpy as np
from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler, RobustScaler, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import cohen_kappa_score as coh_kap
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn.metrics import precision_score
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import permutation_test_score
import pandas as pd
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
import sys
from sklearn.externals import joblib
from sklearn import preprocessing
import operator
import copy
from random import shuffle
import scipy

classes = []
data = []
names = []

# pre-processing of data starts here. All of it is done in pandas.

real_names = {}
fil=open("real_headers.tsv")
for lin in fil.readlines():
    flds = lin.strip().split("\t")
    real_names[flds[0]] = flds[1]
fil.close()


fil=open(sys.argv[1])
canc = sys.argv[1].replace(".train.csv","")
df = pd.read_csv(sys.argv[1],index_col=0)

# You can change this to one-hot encoder. But label encoder seems to work better.

le = preprocessing.LabelEncoder()
onhot = preprocessing.OneHotEncoder(sparse=False)
lab_bin = preprocessing.LabelBinarizer()
df["subsys"] = lab_bin.fit_transform(df["subsys"])
df["path_label"] = lab_bin.fit_transform(df["path_label"])


# The names etc for making the HTML are set based on the target given as argv[2]

targets = set(['TCGA_annot', 'SURV','CNV'])
targ = set([sys.argv[2]])
excl_targ = list(targets - targ)
targ = list(targ)
targ = sys.argv[2]

if(sys.argv[2] == "CNV"):
    targ_labels = ["GAIN","NEUT","LOSS"]
    targ_dict = {'NEUT': 0, 'LOSS': 0, 'GAIN': 0}
else:
    targ_labels = ["UPREG","NEUTRAL","DOWNREG"]
    targ_dict = {'NEUTRAL': 0, 'DOWNREG': 0, 'UPREG': 0}

# Getting the classes and all that. We can also drop a few columns in case we have an exclusion list.

if(len(sys.argv) > 3):
    fil=open(sys.argv[3])
    drop_col_names = [i.strip() for i in fil.readlines()]
    fil.close()
    df = df.drop(columns=drop_col_names)
    #print "Excluding the following variables: "+", ".join(drop_col_names)

df = df.drop(columns=excl_targ)
classes = df[sys.argv[2]]
header=df.columns
df1 = df.copy(deep=True)
df = df.drop(columns=targ)

fir_df = df1.loc[df1[targ] == targ_labels[0]]
sec_df = df1.loc[df1[targ] == targ_labels[1]]
trd_df = df1.loc[df1[targ] == targ_labels[2]]

column_squigly = {}
for col in df.columns:
    v1 = fir_df[col].median()
    v2 = sec_df[col].median()
    v3 = trd_df[col].median()
    correl = np.corrcoef([v1,v2,v3],[1.0,0.0,-1.0])
    if(np.isnan(correl[0][1]) != True):
        column_squigly[col] = correl[0][1]
    else:
        column_squigly[col] = 0.0

data = np.array(df).astype(np.float)
#classes = np.array(classes).astype(np.float)
#data = MinMaxScaler().fit_transform(data)
#data = StandardScaler().fit_transform(data)
data = RobustScaler().fit_transform(data)

X_train, X_test, y_train, y_test = train_test_split(data, classes, test_size=0.3)

#X_train, X_disc, y_train, y_disc = train_test_split(X_train, y_train, test_size=0.5)

data = X_train
classes = y_train

orig_data = X_test
orig_classes = y_test

ros = RandomOverSampler()
X_resampled, y_resampled = ros.fit_sample(data, classes)
data = X_resampled
classes = y_resampled

fil=open("header")
for lin in fil.readlines():
    print lin.strip()
fil.close()

# Random Forest
#print "Random Forest predictions"
#print "Type\t# of features\t# of trees\tKappa score\tF1 (micro)\tPrecision (micro)\tRecall (micro)\tCV average accuracy\tCV stddev accuracy"
feat = 130
while(feat < 140):
    trees = 500
    while(trees <= 500):
        clf = RandomForestClassifier(n_estimators=trees, max_features=feat)
        clf1 = RandomForestClassifier(n_estimators=10, max_features=feat)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        clf.fit(data, classes)
        orig_pred = clf.predict(orig_data)
        new_pred = clf.predict(data)
        conf_matr = confusion_matrix(classes,new_pred,labels=targ_labels)
        conf_matr_norm = conf_matr.astype('float') / conf_matr.sum(axis=1)[:, np.newaxis]

        #score, permutation_scores, pvalue = permutation_test_score(clf1, orig_data, orig_classes, scoring="accuracy", cv=5, n_permutations=50, n_jobs=1)
        pvalue= 0.0

        acc_distri = []
        x=0
        while(x<= 10000):
            temp_classes = list(classes)
            shuffle(temp_classes)
            acc_distri.append(precision_score(temp_classes,classes,average='micro'))
            x=x+1
        ave = np.mean(np.array(acc_distri))
        std = np.std(np.array(acc_distri))
            
        kap = coh_kap(orig_classes,orig_pred)
        f1 = f1_score(orig_classes,orig_pred,average='micro')
        prec = precision_score(orig_classes,orig_pred,average='micro')
        rec = recall_score(orig_classes,orig_pred,average='micro')

        new_pvalue = scipy.stats.norm.sf(abs((prec - ave)/std))
        pvalue = (prec - ave)/std
        #print new_pvalue
        #print (prec - ave)/std
        
        print "var result_vars = ["

        print "{\"naam\":\"CV F1 (mean)\",\"value\":"+str(np.mean(cv_score))+"},"
        print "{\"naam\":\"CV F1 (std dev)\",\"value\":"+str(np.std(cv_score))+"},"
        print "{\"naam\":\"Z-score of accuracy\",\"value\":"+str(pvalue)+"},"
        print "{\"naam\":\"Holdout precision\",\"value\":"+str(prec)+"},"
        print "{\"naam\":\"Holdout recall\",\"value\":"+str(rec)+"}"

        #print "\t{\"cv_f1_mean\":"+str(np.mean(cv_score))+",\"cv_f1_std\":"+str(np.std(cv_score))+",\"holdout_prec\":"+str(prec)+",\"holdout_rec\":"+str(rec)+",\"p-value\":"+str(pvalue)+"}"
            
        print "];"
        temp_dict_feat = {}
        for i, j in zip(header, clf.feature_importances_):
            temp_dict_feat[i] = j
        sorted_d = sorted(temp_dict_feat.items(),key=operator.itemgetter(1),reverse=True)
        print "var feat_imp=["
        x=0
        matrix = np.zeros((12,12))
        z=3
        chord_names = copy.deepcopy(targ_labels)
        while(x<70):
            tempa = sorted_d[x]
            print "\t{\"xval\":"+str(x)+",\"yval\":"+str(tempa[1])+",\"name\":\""+str(real_names[tempa[0]])+"\",\"correl\":"+str(column_squigly[tempa[0]])+"},"
            absmax = df1[tempa[0]].quantile(0.33)
            absmin = df1[tempa[0]].min()
            if(absmin < 0.0):
                fir_df = df1.loc[df1[tempa[0]] <= -1.0*absmax ]
                chord_names.append(real_names[tempa[0]]+" (Low)")
                sec_df = df1.loc[(df1[tempa[0]] > -1.0*absmax) & (df1[tempa[0]] < absmax)]
                chord_names.append(real_names[tempa[0]]+" (Med)")
                trd_df = df1.loc[df1[tempa[0]] >= absmax]
                chord_names.append(real_names[tempa[0]]+" (High)")
            else:
                fir_df = df1.loc[df1[tempa[0]] <= absmax ]
                chord_names.append(real_names[tempa[0]]+" (Low)")
                sec_df = df1.loc[(df1[tempa[0]] > absmax) & (df1[tempa[0]] < 2.0*absmax)]
                chord_names.append(real_names[tempa[0]]+" (Med)")
                trd_df = df1.loc[df1[tempa[0]] > absmax]
                chord_names.append(real_names[tempa[0]]+"(High)")
            if(x < 3):
                fir_dc = fir_df.groupby(targ)[targ].count().to_dict()
                sec_dc = sec_df.groupby(targ)[targ].count().to_dict()
                trd_dc = trd_df.groupby(targ)[targ].count().to_dict()
                for trg in targ_labels:
                    if(trg not in fir_dc):
                        fir_dc = targ_dict
                    if(trg not in sec_dc):
                        sec_dc = targ_dict
                    if(trg not in trd_dc):
                        trd_dc = targ_dict
                matrix[z+0][0] = fir_dc[targ_labels[0]]
                matrix[0][z+0] = fir_dc[targ_labels[0]]
                matrix[z+0][1] = fir_dc[targ_labels[1]]
                matrix[1][z+0] = fir_dc[targ_labels[1]]
                matrix[z+0][2] = fir_dc[targ_labels[2]]
                matrix[2][z+0] = fir_dc[targ_labels[2]]

                matrix[z+1][0] = sec_dc[targ_labels[0]]
                matrix[0][z+1] = sec_dc[targ_labels[0]]
                matrix[z+1][1] = sec_dc[targ_labels[1]]
                matrix[1][z+1] = sec_dc[targ_labels[1]]
                matrix[z+1][2] = sec_dc[targ_labels[2]]
                matrix[2][z+1] = sec_dc[targ_labels[2]]

                matrix[z+2][0] = trd_dc[targ_labels[0]]
                matrix[0][z+2] = trd_dc[targ_labels[0]]
                matrix[z+2][1] = trd_dc[targ_labels[1]]
                matrix[1][z+2] = trd_dc[targ_labels[1]]
                matrix[z+2][2] = trd_dc[targ_labels[2]]
                matrix[2][z+2] = trd_dc[targ_labels[2]]
            x=x+1
            z=z+3
        print "];"

        print "var chord_matrix = "
        print matrix.tolist()
        print ";"

        print "var Names = ["+",".join(["\""+i+"\"" for i in chord_names])+"];"
        print "var labelsData = ["+",".join(["\""+i+"\"" for i in targ_labels])+"];"

        print "var conf_matr = "
        print conf_matr.tolist()
        print ";"
        
        class_max = str(np.amax(conf_matr))

        print "var colorMap = d3.scaleLinear()\n\t.domain([0.0,"+class_max+"])\n\t.range([\"#f7fbff\", \"#08306b\"]);\n"
 

        # We need to print a matrix of xy values including zeros. First three can be GAIN/NEUT/LOSS
        # Take each of the top features, and get 0.33 and 0.66 quantile (df.quantile([0.33,0.66])
        # Take all values < 0.33 and find the number of GAIN/NEUT/LOSS.
        # Repeat for values 0.33 to 0.66 and >0.66. All feature-feature values should be zero.
        # Same colors for the features? Lets see. 
        trees = trees + 1500
    feat = feat + 20

fil=open("tailer_chord")
for lin in fil.readlines():
    print lin.strip()
fil.close()

