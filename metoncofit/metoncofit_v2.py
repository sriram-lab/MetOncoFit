"""

"""

import datetime
import numpy as np
import pandas as pd


class MetOncoFit:
    def __init__(self, X, y, nTrees, nFeatures, size,
                 depth=10, minLeafs=5):

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
