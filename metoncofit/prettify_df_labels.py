"""
"""


def prettifyCancerLabels(file=sys.argv[1]):
    canc = fil3.replace(".csv", "")
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


def prettify_DataFrame_colNames():
    columnName_map = pd.read_csv("./../labels/real_headers.txt",
                                 sep='\t', names=['Original', 'New'])
    bestNames = dict([(i, name)
                      for i, name in zip(columnName_map['Original'], columnName_map['New'])])
    return bestNames


def prettify_prediction_labels(prediction):
    """
    """
    if(prediction == "CNV"):
        prediction_labels = ["GAIN", "NEUT", "LOSS"]
        prediction_dict = {'NEUT': 0, 'LOSS': 0, 'GAIN': 0}
    else:
        prediction_labels = ["UPREG", "NEUTRAL", "DOWNREG"]
        prediction_dict = {'NEUTRAL': 0, 'DOWNREG': 0, 'UPREG': 0}

    return prediction_labels, prediction_dict
