"""
"""


def long_cancer_nnames(fileName):
    canc = fileName.replace(".csv", "")
    canc_dict = {
        'breast': 'Breast Cancer',
        'cns': 'Glioma',
        'colon': 'Colorectal Cancer',
        'complex': 'Pan Cancer',
        'leukemia': 'B-cell lymphoma',
        'melanoma': 'Melanoma',
        'nsclc': 'Lung Cancer',
        'ovarian': 'Ovarian Cancer',
        'prostate': 'Prostate Cancer',
        'renal': 'Renal Cancer'
    }
    canc = canc_dict.get(canc)
    return canc


def long_feature_names():
    labelFile = "./../labels/real_headers.txt"
    columnName_map = pd.read_csv(
        labelFile, sep='\t', names=['Original', 'New'])
    bestNames = dict([(index, name)
                      for index, name in zip(columnName_map['Original'], columnName_map['New'])])
    return bestNames


def set_prediction_labels(target):
    """
    """
    if(target == "CNV"):
        prediction_labels = ["GAIN", "NEUT", "LOSS"]
        prediction_dict = {'NEUT': 0, 'LOSS': 0, 'GAIN': 0}
    else:
        prediction_labels = ["UPREG", "NEUTRAL", "DOWNREG"]
        prediction_dict = {'NEUTRAL': 0, 'DOWNREG': 0, 'UPREG': 0}

    return prediction_labels, prediction_dict
