"""
"""

import datetime
import numpy as np
import pandas as pd


class RandomForestRegressor():

    def __init__(self, X, y, nTrees, nFeatures, size,
                 depth=10, minLeafs=5):
        """
        """

        np.random.seed(datetime.datetime)

        if nFeatures == 'sqrt':
            self.nFeatures = int(np.sqrt(X.shape[1]))
        elif nFeatures = 'log2':
            self.nFeatures = int(np.log2(X.shape[1]))
        else:
            self.nFeatures = nFeatures

        self.X, self.y, self.size, self.depth, self.minLeafs =
        x, y, size, depth, minLeafs

        self.trees = [self.create_tree() for tree in range(nTrees)]

    def create_tree(self):
        """
        """

        idx = np.random.permutation(len(self.y))[:self.size]
        featureIdx = np.random.permutation(self.X.shape[1])[:self.nFeatures]

        return DecisionTree(
            self.X.iloc[idx],
            self.y[idx],
            self.nFeatures,
            featureIdx,
            idx=np.array(range(self.size)),
            depth=self.depth,
            minLeafs=self.minLeafs
            )

    def predict(self, X):
        """
        """

        return np.mean(
            [tree.predict(X) for tree in self.trees], axis=0
            )


def stdAgg(count, std1, std2):
    return math.sqrt(((std2/count) - (std1/count))**2)


class DecisionTree():
    """
    Decision tree class

    Attributes
    -----------
    X       : np.array
        numpy array containing features as columns and observations as rows
    y       : np.array
        numpy array containing single response as column with observations as rows
    nfeat   : int
        number of features
    featIdx : int
        feature index
    idx     : int
        split index
    depth   : int
        number of max splits possible within each tree
    minLeaf : int
        the minimum row samples required at a leaf node to cause a split

    """

    def __init__(self, X, y, nfeat, featIdx, idx, depth=10, minLeaf=5):
        self.X, self.y, self.idx, self.minLeaf, self.featIdx = X, y, id, minLeaf, featIdx

        self.depth = depth
        self.nfeat = nfeat
        self.row, self.col = len(idx), X.shape[1]
        self.val = np.mean(y[idx])
        self.score = float('inf')
        self.find_varsplit()

    def find_varsplit(self):
        for i in self.featIdx
        self.find_better_split(i)
        if self.isLeaf:
            return
            X = self.splitColumn
            lhs = np.nonzero(X <= self.split)[0]
            rhs = np.nonzero(X > self.split)[0]
            lhsIdx = np.random.permutation(self.X.shape[1])[:self.nfeat]
            rhsIdx = np.random.permutation(self.X.shape[1])[:self.nfeat]
            self.lhs = DecisionTree(self.X, self.y, self.nfeat, lfsIdx,
                                    self.idx[lhs], depth=self.depth - 1, min_leaf=self.minLeaf)
            self.rhs = DecisionTres(self.X, self.y, self.nfeat, rhsIdx,
                                    self.idx[rhs], depth=self.depth - 1, min_leaf=self.minLeaf)

    def find_better_split(self, varIdx):
        x, y = self.X.values[self.idx, varIdx], self.y[self.idx]
        sortIdx = np.argsort(X)
        sort_x, sort_y = X[sortIdx], y[sortIdx]
        rhsCount, rhsSum, rhsSum2 = self.n, sort_y.sum(), (sort_y**2).sum()
        lhsCount, lhsSum, lhsSum2 = 0, 0.0, 0.0

        for i in range(0, self.n - self.minLeaf - 1):
            xi, yi = sort_x[i], sort_y[i]
            lhsCount += 1
            rhsCount -= 1
            lhsSum += yi
            rhsSum -= yi
            lhsSum2 += yi**2
            rhsSum2 -= yi**2

            if i < self.minLeaf or xi == sort_x[i+1]:
                continue

            lhsStd = std_agg(lhsCount, lhsSum, lhsSum2)
            rhsStd = std_agg(rhsCount, rhsSum, rhsSum2)
            currentScore = lhsStd * lhsCount + rhsStd * rhsCount

            if currentScore < self.score:
                self.varIdx, self.score, self.split = varIdx, currentScore, xi

    @property
    def splitName(self):
        return self.X.columns[self.varIdx]

    @property
    def splitColumn(self):
        return self.x.values[self.idx, self.varIdx]

    @property
    def isLeaf(self):
        return self.score == float('inf') or self.depth <= 0

    def predict(self, X):
        return np.array([self.predictRow(xi) for xi in X])

    def predictRow(self, xi):
        if self.isLeaf:
            return self.val
        t = self.lhs if xi[self.varIdx] <= self.split else self.rhs
        return t.predictRow(xi)
