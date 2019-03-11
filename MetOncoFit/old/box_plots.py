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
import pandas as pd
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
import sys
from sklearn.externals import joblib
from sklearn import preprocessing
import operator
import seaborn as sns
from matplotlib import pyplot

classes = []
data = []
names = []
fil=open(sys.argv[1])
canc = sys.argv[1].replace(".train.csv","")
#x=0
#for lin in fil.readlines():
#    if(x==0):
#        x=1
#        header = lin.strip().split(",")[1:-2]
#    else:
#        flds = lin.strip().split(",")
#        classes.append(flds[-2])
#        data.append(flds[1:-3])
#        names.append(flds[0])
#fil.close()
sns.set_style("whitegrid")
df = pd.read_csv(sys.argv[1],index_col=0)
sns.boxplot(y="tot_biom",x="path_label",data=df,palette="Set2")
pyplot.show()

le = preprocessing.LabelEncoder()
le.fit(df["subsys"])
df["subsys"] = le.transform(df["subsys"])
le.fit(df["path_label"])
df["path_label"] = le.transform(df["path_label"])

df.to_csv("temp.csv")
df = df.drop(columns=['TCGA_annot', 'SURV'])
classes = df['CNV']
header=df.columns
df1 = df.copy(deep=True)
df = df.drop(columns=['CNV'])

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

# Random Forest
#print "Random Forest predictions"
#print "Type\t# of features\t# of trees\tKappa score\tF1 (micro)\tPrecision (micro)\tRecall (micro)\tCV average accuracy\tCV stddev accuracy"
feat = 130
while(feat < 140):
    trees = 10
    while(trees <= 10):
        clf = RandomForestClassifier(n_estimators=trees, max_features=feat)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        clf.fit(data, classes)
        orig_pred = clf.predict(orig_data)
        kap = coh_kap(orig_classes,orig_pred)
        f1 = f1_score(orig_classes,orig_pred,average='micro')
        prec = precision_score(orig_classes,orig_pred,average='micro')
        rec = recall_score(orig_classes,orig_pred,average='micro')
        #print "RF\t"+str(feat)+"\t" +str(trees)+"\t"+str(kap)+"\t"+str(f1)+"\t"+str(prec)+"\t"+str(rec)+"\t\t"+str(np.mean(cv_score))+"\t"+str(np.std(cv_score))
        print "var result_vars = ["
        print "\t{\"cv_f1_mean\":"+str(np.mean(cv_score))+",\"cv_f1_std\":"+str(np.std(cv_score))+",\"holdout_prec\":"+str(prec)+",\"holdout_rec\":"+str(rec)+"}"
        print "];"
        temp_dict_feat = {}
        for i, j in zip(header, clf.feature_importances_):
            temp_dict_feat[i] = j
        sorted_d = sorted(temp_dict_feat.items(),key=operator.itemgetter(1),reverse=True)
        plot_vars = [tempa[0] for tempa in sorted_d[:5]]
        print plot_vars
        melted_df = df1.melt(id_vars=plot_vars,var_name="Important features")
        print melted_df.columns
        trees = trees + 1500
    feat = feat + 40
