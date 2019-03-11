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

classes = []
data = []
names = []
fil=open(sys.argv[1])
canc = sys.argv[1].replace(".train.csv","")
x=0
for lin in fil.readlines():
    if(x==0):
        x=1
        header = lin.strip().split(",")[1:-2]
    else:
        flds = lin.strip().split(",")
        classes.append(flds[-3])
        data.append(flds[1:-5])
        names.append(flds[0])
fil.close()

data = np.array(data).astype(np.float)
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

wrf=open(canc+"_diffexp_feature_importances.txt","w")
# Random Forest
print "Random Forest predictions"
print "Type\t# of features\t# of trees\tKappa score\tF1 (micro)\tPrecision (micro)\tRecall (micro)\tCV average accuracy\tCV stddev accuracy"
feat = 100
while(feat <= 100):
    trees = 500
    while(trees <= 500):
        clf = RandomForestClassifier(n_estimators=trees, max_features=feat)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        clf.fit(data, classes)
        orig_pred = clf.predict(orig_data)
        kap = coh_kap(orig_classes,orig_pred)
        f1 = f1_score(orig_classes,orig_pred,average='micro')
        prec = precision_score(orig_classes,orig_pred,average='micro')
        rec = recall_score(orig_classes,orig_pred,average='micro')
        print "RF\t"+str(feat)+"\t" +str(trees)+"\t"+str(kap)+"\t"+str(f1)+"\t"+str(prec)+"\t"+str(rec)+"\t\t"+str(np.mean(cv_score))+"\t"+str(np.std(cv_score))
        wrf.write(str(feat)+"\t"+str(trees)+"\t")
        for i, j in zip(header, clf.feature_importances_): wrf.write(i +";;"+ str(j)+"\t")
        wrf.write("\n")
        trees = trees + 100
    feat = feat + 10
wrf.close()

sys.exit(1)
# Random Forest
print "SVM predictions"
print "Type\tgamma\tC-value\tKappa score\tF1 (micro)\tPrecision (micro)\tRecall (micro)\tCV average accuracy\tCV stddev accuracy"
gam = 0.001
while(gam <= 10.0):
    cval = 1.0
    while(cval <= 10000.0):
        clf = clf = svm.SVC(kernel='rbf',gamma=gam,C=cval,max_iter=5000)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        clf.fit(data, classes)
        orig_pred = clf.predict(orig_data)
        kap = coh_kap(orig_classes,orig_pred)
        f1 = f1_score(orig_classes,orig_pred,average='micro')
        prec = precision_score(orig_classes,orig_pred,average='micro')
        rec = recall_score(orig_classes,orig_pred,average='micro')
        print "SVM\t"+str(gam)+"\t" +str(cval)+"\t"+str(kap)+"\t"+str(f1)+"\t"+str(prec)+"\t"+str(rec)+"\t\t"+str(np.mean(cv_score))+"\t"+str(np.std(cv_score))
        cval = cval *10.0
    gam = gam*10.0


# Neural network
print "Neural net predictions"
print "Type\tAlpha\t# of iterations\tKappa score\tF1 (micro)\tPrecision (micro)\tRecall (micro)\tCV average accuracy\tCV stddev accuracy"
hidden_layer_sizes = [(2048,16),(1024,16),(512,16),(2048,64),(1024,64),(512,64),(2048,128),(1024,128),(512,128)]
alp = 0.01
while(alp < 0.3):
    itr = 50
    while(itr < 500):
        clf = MLPClassifier(solver='adam', alpha=alp, hidden_layer_sizes=(1024,128),max_iter=itr)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        clf.fit(data, classes)
        orig_pred = clf.predict(orig_data)
        kap = coh_kap(orig_classes,orig_pred)
        f1 = f1_score(orig_classes,orig_pred,average='micro')
        prec = precision_score(orig_classes,orig_pred,average='micro')
        rec = recall_score(orig_classes,orig_pred,average='micro')
        print "NN\t"+str(alp)+"\t" +str(itr)+"\t"+str(kap)+"\t"+str(f1)+"\t"+str(prec)+"\t"+str(rec)+"\t\t"+str(np.mean(cv_score))+"\t"+str(np.std(cv_score))
        itr=itr+50
    alp=alp*2.0


