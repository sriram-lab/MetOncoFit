def buildModel(trainingDataFilePath, testingDataFilePath, hyperparamSaveFile):
    # The purpose of this function is to build a random forest model to predict differentially
    # expressed genes given kinetics data, flux data, and model topological features

    # This function performs automated hyper-parameter tuning and stores the results from this
    # tuning in a file as given in the input

    # The function will also display the accuracy, confusion matrix, Matthews correlation
    # coefficient, and a classification report from sci-kit learn.

    # Two figures will be created, a normalized confusion matrix and a bar plot showing the accuracy
    # when removing each subset of data from the model as a feature importance plot

    # Inputs:
    # trainingDataFilePath: the file path for the training data
    # testingDataFilePath: the file path for the testing data
    # hyperparamSaveFile: a csv file name for saving the results from each
    # step of the automated hyper-parameter tuning

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

    # Hyper-parameter tuning

    # This step imports the training data and removes columns that are not wanted for training
    # (reg - the regulation label ie. up, down, neutral, gene - the gene name, FC - fold
    # change value)
    all_data = pd.read_excel(trainingDataFilePath, sheet_name='all_data')
    X = all_data.drop(['reg', 'gene', 'FC'], axis=1)

    # This step creates a label encoder to assign labels to up, down, and neutral
    LE = preprocessing.LabelEncoder()
    LE.fit(all_data.reg)
    y = LE.transform(all_data.reg)

    # This step create a random over sampler to balance the classes for training
    ROS = RandomOverSampler(random_state=1)
    X_OS, y_OS = ROS.fit_resample(X, y)

    # This function defines the objective function for hyper-parameter tuning
    # This function aims to minimize accuracy and recall during 5 fold cross validation

    # Inputs:
    # params - a dictionary containing the parameters to tune
    # n_folds - the number of folds for cross validation
    # out_file - the file path for saving each step of the hyper-parameter
    # tuning process

    # Outputs:
    # loss - the loss calculated by the loss function
    # params - the dictionary of parameters
    # status - returns the status of the tuning (OK, error, etc)

    def objective(params, n_folds=5, out_file=hyperparamSaveFile):

        # maximize recall and accuracy
        # sets certain parameters as integers as necessary for use in RandomForestClassifier
        params['max_leaf_nodes'] = int(params['max_leaf_nodes'])
        params['min_samples_leaf'] = int(params['min_samples_leaf'])
        model = RandomForestClassifier(**params)
        accuracies = cross_val_score(model, X_OS, y_OS, cv=n_folds)
        recall = cross_val_score(model, X_OS, y_OS, cv=n_folds, scoring='recall_micro')
        loss = 2-(max(accuracies)+max(recall))

        # adds row to the file containing the loss and parameters
        file = open(out_file, 'a')
        writer = csv.writer(file)
        writer.writerow([loss, params])
        file.close()

        return {'loss': loss, 'params': params, 'status': STATUS_OK}

    # This defines the search space for each hyper-parameter
    # qloguniform outputs an integer from a log distribution between the values given (ie. e^1.5 - e^4)
    # choice is used for hyper-parameters that are non-numeric
    # quniform outputs an integer from a uniform distribution between the values given
    # uniform outputs a decimal value from a uniform distribution between the values given
    space = {
        'n_estimators': hp.qloguniform('n_estimators', 1.5, 4, 1),
        'criterion': hp.choice('criterion', ['gini', 'entropy']),
        'max_depth': hp.quniform('max_depth', 1, 15, 1),
        'min_samples_leaf': hp.quniform('min_samples_leaf', 10, 100, 1),
        'max_features': hp.uniform('max_features', 0, 1)
    }

    trials = Trials()

    # The best set of hyper-parameters is outputted from this function. The number of
    # max evals can be adjusted.
    best = fmin(fn=objective, space=space, algo=tpe.suggest, max_evals=50, trials=trials)
    print('Hyper-parameters: ', best)

    # Test model
    train_data = pd.read_excel(trainingDataFilePath, sheet_name='all_data')
    test_data = pd.read_excel(testingDataFilePath, sheet_name='all_data')

    # The classifier is defined using the best set of hyper-parameters determined from tuning
    clf = RandomForestClassifier(**best)

    y_train = LE.transform(train_data.reg)
    y_test = LE.transform(test_data.reg)
    reg_classes = LE.classes_

    # removing unwanted predictors
    x_train = train_data.drop(['reg', 'gene', 'FC'], axis=1)
    x_test = test_data.drop(['reg', 'gene', 'FC'], axis=1)

    x_train_OS, y_train_OS = ROS.fit_resample(x_train, y_train)

    # Training the data and predicting based on the testing data
    clf.fit(x_train_OS, y_train_OS)
    y_pred = clf.predict(x_test)

    # This portion is used to output desired metrics and plot the confusion matrix
    print('The accuracy is: ', accuracy_score(y_pred, y_test))
    print('The confusion matrix is: \n', confusion_matrix(y_test, y_pred))
    print('The Matthews correlation coefficient is: ', matthews_corrcoef(y_test, y_pred))
    print(classification_report(y_test, y_pred, digits=3, target_names=reg_classes))
    plot_confusion_matrix(clf, x_test, y_test, display_labels=reg_classes, normalize='true', cmap=plt.cm.Blues)
    plt.show()

    # test for over fitting
    y_pred_train = clf.predict(x_train_OS)
    print('The accuracy for the training data is: ', accuracy_score(y_pred_train, y_train_OS))

    # Leave one subset out cross validation

    # The list of sheet names is used for removing each subset of data from the training data
    # An empty accuracy vector is created to store the result from each round
    sheet_names = ['KO_flux', 'redist_flux', 'kinetics', 'distances']
    accuracies = np.empty([len(sheet_names)])
    count = 0
    for sheet in sheet_names:

        # In addition to the subset being removed, reg, gene, and FC will be removed as in
        # the model prior
        rmv_data = pd.read_excel(trainingDataFilePath, sheet_name=sheet)
        x_train = train_data.drop(rmv_data.columns, axis=1)
        x_train = x_train.drop(['reg', 'gene', 'FC'], axis=1)
        x_test = test_data.drop(rmv_data.columns, axis=1)
        x_test = x_test.drop(['reg', 'gene', 'FC'], axis=1)

        y_train = train_data.reg
        y_test = test_data.reg

        # The data is oversampled in the same manor as the prior model
        x_train_OS, y_train_OS = ROS.fit_resample(x_train, y_train)

        clf.fit(x_train_OS, y_train_OS)
        y_pred = clf.predict(x_test)

        # The accuracy of the model is stored here
        accuracies[count] = accuracy_score(y_pred, y_test)
        count += 1

    # This portion is creating and formatting the bar chart
    plt.figure()
    axisRange = np.arange(len(sheet_names))
    plt.bar(axisRange, accuracies, align='center')
    plt.xticks(axisRange, sheet_names, rotation=45, fontsize=8)
    plt.title('Leave one feature set out')
    plt.xlabel('Feature left out')
    plt.ylabel('Accuracy')
    plt.show()
