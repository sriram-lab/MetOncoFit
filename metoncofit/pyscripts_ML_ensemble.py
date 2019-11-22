# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 20:43:57 2018
Machine Learning Ensemble
@author: Krishna Oruganty
"""

import numpy as np
from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler, RobustScaler
from sklearn.model_selection import train_test_split
import pandas as pd
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from sklearn.preprocessing import LabelEncoder,OneHotEncoder,LabelBinarizer
import sys, os

def 
for fname in os.listdir('.'):
    if fname.endswith('.csv'):
        dataset= pd.read_csv(fname)
        le = LabelEncoder()
        onhot = OneHotEncoder(sparse=False)
        lab_bin = LabelBinarizer()
        dataset["subsys"] = lab_bin.fit_transform(dataset["subsys"])
        dataset["path_label"] = lab_bin.fit_transform(dataset["path_label"])
        X=dataset.iloc[:,1:(len(dataset.columns)-1)].values
        y=dataset['SURV'].values

        X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.2)

        data = X_train
        classes = y_train

        orig_data = X_test
        orig_classes = y_test

        ros = RandomOverSampler()
        X_resampled, y_resampled = ros.fit_sample(data, classes)

        #print classes
        data = X_resampled
        classes = y_resampled

        #for c in classes:
        #    print c
        # SVM params
        print ("SVM predictions")
        clf = svm.SVC(kernel='rbf',gamma=1000.0, C=100.0,max_iter=500)
        clf.fit(data, classes)
        svm_pred = clf.predict(orig_data)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        print (cv_score)
        print (str(np.mean(cv_score))+"\t"+str(np.std(cv_score)))
        grid = GridSearchCV(svm.SVC(),param_grid={'C':[10, 100, 1000,10000],'gamma':[10.0,100.0,1000.0],'kernel':['linear','rbf','poly'],'max_iter':[100,500,1000],'degree':[2,3,4,5]},scoring='accuracy',cv=5)
        grid.fit(data,classes)
        print(grid.best_params_)
        means = grid.cv_results_['mean_test_score']
        stds = grid.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, grid.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))


        # Neural network
        print ("Neural net predictions")
        clf = MLPClassifier(solver='adam', alpha=0.1, hidden_layer_sizes=(125, 5),max_iter=1000)
        clf.fit(data, classes)
        cnn_pred = clf.predict(orig_data)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        print (cv_score)
        print (str(np.mean(cv_score))+"\t"+str(np.std(cv_score)))
        grid = GridSearchCV(MLPClassifier(),param_grid={'solver':['adam'],'alpha':[0.1,0.2,0.5],'max_iter':[500,1000,5000],'hidden_layer_sizes':[(100,10),(125,10),(200,10),(100,5),(125,5),(200,5)]},scoring='accuracy',cv=5)
        grid.fit(data,classes)
        print(grid.best_params_)
        means = grid.cv_results_['mean_test_score']
        stds = grid.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, grid.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))

        # Random Forest
        print ("Random Forest predictions")
        clf = RandomForestClassifier(n_estimators=100, max_features=8)
        clf.fit(data, classes)
        rf_pred = clf.predict(orig_data)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        print (cv_score)
        print (str(np.mean(cv_score))+"\t"+str(np.std(cv_score)))
        grid = GridSearchCV(RandomForestClassifier(),param_grid={'max_features':[8,9, 10,15,20,25,30,40],'n_estimators':[100,500,1000]},scoring='accuracy',cv=5)
        grid.fit(data,classes)
        print(grid.best_params_)
        means = grid.cv_results_['mean_test_score']
        stds = grid.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, grid.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))
        #print('RF feature_importances:')
        #for i, j in zip(header, clf.feature_importances_): print(i,j)

        # Gradient boosting tree
        print ("Gradient Boosting Tree predictions")
        clf = GradientBoostingClassifier(n_estimators=50, learning_rate=1.0,max_features=40)
        clf.fit(data, classes)
        gbt_pred = clf.predict(orig_data)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        print (cv_score)
        print (str(np.mean(cv_score))+"\t"+str(np.std(cv_score)))
        grid = GridSearchCV(GradientBoostingClassifier(),param_grid={'max_features':[8,9,10,20,25,30,40],'n_estimators':[50,100,500],'learning_rate':[0.1,1.0,2.0,5.0]},scoring='accuracy',cv=5)
        grid.fit(data,classes)
        print(grid.best_params_)
        means = grid.cv_results_['mean_test_score']
        stds = grid.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, grid.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))


        # AdaBoost
        print ("AdaBoost Tree predictions")
        clf = AdaBoostClassifier(n_estimators=500,learning_rate=1.0)
        clf.fit(data, classes)
        ada_pred = clf.predict(orig_data)
        cv_score = cross_val_score(clf, data, classes, scoring='accuracy',cv=5)
        print (cv_score)
        print (str(np.mean(cv_score))+"\t"+str(np.std(cv_score)))
        grid = GridSearchCV(AdaBoostClassifier(),param_grid={'n_estimators':[50,100,500,1000],'learning_rate':[0.1,1.0,2.0,3.0,5.0]},scoring='accuracy',cv=5)
        grid.fit(data,classes)
        print(grid.best_params_)
        means = grid.cv_results_['mean_test_score']
        stds = grid.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, grid.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))

        x=0
        for gen in svm_pred:
            print (names[x]+"\t\t"+str(svm_pred[x])+"\t"+str(cnn_pred[x])+"\t"+str(rf_pred[x])+"\t"+str(gbt_pred[x])+"\t"+str(ada_pred[x])+"\t"+orig_classes[x],x=x+1)
        break
