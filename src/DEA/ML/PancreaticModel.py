def buildModel(trainingDataFilePath, testingDataFilePath, hyperparamSaveFile):
    import pandas as pd
    import numpy as np
    from sklearn import preprocessing
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score, confusion_matrix, matthews_corrcoef, \
        plot_confusion_matrix, classification_report
    from sklearn.model_selection import cross_val_score
    from imblearn.over_sampling import RandomOverSampler
    from hyperopt import hp, STATUS_OK, tpe, fmin, Trials
    import matplotlib.pyplot as plt
    import csv

    # # Hyper-parameter tuning
    all_data = pd.read_excel(trainingDataFilePath, sheet_name='all_data')
    X = all_data.drop(['reg', 'gene', 'FC'], axis=1)

    LE = preprocessing.LabelEncoder()
    LE.fit(all_data.reg)
    y = LE.transform(all_data.reg)

    ROS = RandomOverSampler(random_state=1)
    X_OS, y_OS = ROS.fit_resample(X, y)


    def objective(params, n_folds=5, out_file=hyperparamSaveFile):

        # maximize recall and accuracy
        params['max_leaf_nodes'] = int(params['max_leaf_nodes'])
        params['min_samples_leaf'] = int(params['min_samples_leaf'])
        model = RandomForestClassifier(**params)
        accuracies = cross_val_score(model, X_OS, y_OS, cv=n_folds)
        recall = cross_val_score(model, X_OS, y_OS, cv=n_folds, scoring='recall_micro')
        loss = 2-(max(accuracies)+max(recall))

        file = open(out_file, 'a')
        writer = csv.writer(file)
        writer.writerow([loss, params])
        file.close()

        return {'loss': loss, 'params': params, 'status': STATUS_OK}


    space = {
        'n_estimators': hp.qloguniform('n_estimators', 1.5, 4, 1),
        'criterion': hp.choice('criterion', ['gini', 'entropy']),
        'max_depth': hp.quniform('max_depth', 1, 15, 1),
        'min_samples_leaf': hp.quniform('min_samples_leaf', 10, 100, 1),
        'max_features': hp.uniform('max_features', 0, 1)
    }

    trials = Trials()

    best = fmin(fn=objective, space=space, algo=tpe.suggest, max_evals=50, trials=trials)
    print('Hyper-parameters: ', best)

    # Test model
    train_data = pd.read_excel(trainingDataFilePath, sheet_name='all_data')
    test_data = pd.read_excel(testingDataFilePath, sheet_name='all_data')

    clf = RandomForestClassifier(**best)

    y_train = LE.transform(train_data.reg)
    y_test = LE.transform(test_data.reg)
    reg_classes = LE.classes_

    x_train = train_data.drop(['reg', 'gene', 'FC'], axis=1)
    x_test = test_data.drop(['reg', 'gene', 'FC'], axis=1)

    ROS = RandomOverSampler(random_state=1)
    x_train_OS, y_train_OS = ROS.fit_resample(x_train, y_train)

    clf.fit(x_train_OS, y_train_OS)
    y_pred = clf.predict(x_test)
    print('The accuracy is: ', accuracy_score(y_pred, y_test))
    print('The confusion matrix is: \n', confusion_matrix(y_test, y_pred))
    print('The Matthews correlation coefficient is: ', matthews_corrcoef(y_test, y_pred))

    print(classification_report(y_test, y_pred, digits=3, target_names=reg_classes))

    plot_confusion_matrix(clf, x_test, y_test, display_labels=reg_classes, normalize='true', cmap=plt.cm.Blues)
    plt.show()

    # test for over fitting
    y_pred_train = clf.predict(x_train_OS)
    print('The accuracy for the training data is: ', accuracy_score(y_pred_train, y_train_OS))

    # # Leave one subset out cross validation
    sheet_names = ['FC', 'KO_flux', 'redist_flux', 'kinetics', 'distances']
    accuracies = np.empty([len(sheet_names)])
    count = 0
    for sheet in sheet_names:
        rmv_data = pd.read_excel(trainingDataFilePath, sheet_name=sheet)
        x_train = train_data.drop(rmv_data.columns, axis=1)
        x_train = x_train.drop(['reg', 'gene', 'FC'], axis=1)
        x_test = test_data.drop(rmv_data.columns, axis=1)
        x_test = x_test.drop(['reg', 'gene', 'FC'], axis=1)

        y_train = train_data.reg
        y_test = test_data.reg

        x_train_OS, y_train_OS = ROS.fit_resample(x_train, y_train)

        clf.fit(x_train_OS, y_train_OS)
        y_pred = clf.predict(x_test)
        accuracies[count] = accuracy_score(y_pred, y_test)
        count += 1

    plt.figure()
    axisRange = np.arange(len(sheet_names))
    plt.bar(axisRange, accuracies, align='center')
    plt.xticks(axisRange, sheet_names, rotation=45, fontsize=8)
    plt.title('Leave one feature set out')
    plt.xlabel('Feature left out')
    plt.ylabel('Accuracy')
    plt.show()
