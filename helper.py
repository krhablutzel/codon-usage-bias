# helper.py
# all the helper functions for the analysis in report.ipynb

import pandas as pd
import numpy as np

import math
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.utils import shuffle
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier


def reduce_dim(data_pd, dim=2, drop_columns=[], vertebrate=False):
    # separate out numeric data
    if vertebrate:
        drop_columns.append("Vertebrate")

    data_num = data_pd.drop(columns=drop_columns)
    data_std = data_num.to_numpy()

    # run PCA
    pca_alg = PCA(n_components=dim)
    component_names = [f'Dim{i+1}' for i in range(dim)]  # name PCA columns
    data_2d = pd.DataFrame(pca_alg.fit_transform(data_std), columns=component_names)

    # add back metadata
    metadata = data_pd[drop_columns]
    metadata.reset_index(drop=True, inplace=True)
    data_2d_meta = pd.concat([data_2d, metadata], axis=1)

    return data_2d_meta


def standardize(data):
    """standardize each column of a numpy dataframe"""
    mean_vec = np.mean(data, axis=0)
    sd_vec = np.std(data, axis=0)

    data_std = data.copy()
    for i in range(data.shape[1]):  # for every column
        if math.isclose(sd_vec[i], 0):  # account for data w/ standard dev of 0
            data_std[:, i] = data[:, i]
        else:
            data_std[:, i] = (data[:, i] - mean_vec[i] * np.ones(data.shape[0])) / sd_vec[i]

    return data_std


def normalize(data):
    """normalize each column of a numpy dataframe:
    set range between 0 and 1"""
    columns = []
    for i in range(data.shape[1]):  # for every column
        # find max and min of column
        col = data[:, 1]
        mx = np.max(col)
        mn = np.min(col)

        # scale column between 0 and 1
        if math.isclose(mx - mn, 0):
            col_norm = col
        else:
            col_norm = (col - mn) / (mx - mn)

        columns.append(col_norm)

    data_norm = np.stack(columns, axis=-1)

    return data_norm


def kmeans_sklearn(data, k):
    """returns labels for k clusters in the data"""
    km_alg = KMeans(n_clusters=k)
    fit1 = km_alg.fit(data)
    labels = fit1.labels_
    centers = fit1.cluster_centers_
    return labels, centers


# classifier methods from CSC294 final portfolio
def divide_data(data):
    """divide dataset into two sets: 90% test/train and 10% validation"""
    n = data.shape[0]

    # take out 10% of the data for validation
    # https://numpy.org/doc/stable/reference/random/generated/numpy.random.choice.html
    ind_valid = np.random.choice(n, size=n // 10, replace=False)
    data_valid = data.iloc[ind_valid]

    # take the other 90% for building the model
    # https://stackoverflow.com/questions/27824075/accessing-numpy-array-elements-not-in-a-given-index-list
    ind_tt = [x for x in range(n) if x not in ind_valid]  # not in index
    data_tt = data.iloc[ind_tt]

    # shuffle data for test/train so no patterns in folds
    # https://stackoverflow.com/questions/29576430/shuffle-dataframe-rows
    data_tt = shuffle(data_tt)

    return data_valid, data_tt


def classification_mse(class_truth, pred_class):
    """compute classification mse"""
    return np.mean(class_truth != pred_class)


def cross_validation(data, method, k):
    """k-fold cross-validation"""
    # calculate fold divisions
    n = data.shape[0]
    n_predictors = data.shape[1] - 1
    foldSize = n // k  # int divide
    foldDivisions = [foldSize * x for x in range(k + 1)]

    # adjust for uneven fold size
    if n % k != 0:
        r = n % k  # remainder
        for i in range(1, k + 1):
            # add 1 + previous size increase to each group until r
            # then just shift by r to account for previous size increases
            foldDivisions[i] += min(i, r)

    # divide into folds
    folds = []
    for i in range(k):
        folds.append(data.iloc[foldDivisions[i]:foldDivisions[i + 1], :])

    # linear model w/ each fold as test once
    test_errors = []

    for i in range(k):
        # get test fold
        test = folds[i]

        # combine other folds into training set
        train_folds = folds.copy()
        train_folds.pop(i)
        # https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html
        train = pd.concat(train_folds)  # concatenate folds

        if method == 'SVM':
            # fit SVM to training data
            mod = SVC(kernel="sigmoid")  # or "rbf" kernel
        elif method == 'neighbor':
            # fit kNN to training data
            mod = KNeighborsClassifier(n_neighbors=9)
        elif method == 'tree':
            # fit decision tree to training data
            mod = DecisionTreeClassifier()
        elif method == 'forest':
            # build random forest classifier for training data
            mod = RandomForestClassifier(n_estimators=20, max_features=min(3, n_predictors), max_depth=3,
                                         random_state=0)

        mod.fit(train.iloc[:, :-1], train.iloc[:, -1])  # class var in last column

        # compute testing error
        test_preds = mod.predict(test.iloc[:, :-1])
        test_error = classification_mse(test_preds, test.iloc[:, -1])
        test_errors.append(test_error)

    # cross validation error - avg of test errors
    cross_val_error = np.mean(test_errors)

    return cross_val_error


def all_cv_errors(weather_tt, methods):
    """get cross-validation error for all possible models"""
    cv_errors = []
    # test each possible model type
    for method in methods:
        # compute cross-validation error
        cv_err = cross_validation(weather_tt, method, 10)

        # store errors
        cv_errors.append([method, cv_err])

    # sort cv errors w/ lowest in first row
    cv_err_np = np.array(cv_errors, dtype=object)
    cv_err_np = cv_err_np[np.argsort(cv_err_np[:, 1])]
    return cv_err_np
