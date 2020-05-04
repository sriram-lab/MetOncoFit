

def make_hr_statistics_table(file, targ):
    """
    """
    hr_dataset_type = file
    cancer_tissue = file.split('.')[0].split('/')[-1]
    if './../data/lax/' in hr_dataset_type:
        hr_threshold = "[0.90 - 1.10]"
    elif './../data/median/' in hr_dataset_type:
        hr_threshold = "[0.75 - 1.33]"
    elif './../data/stringent/' in hr_dataset_type:
        hr_threshold = "[0.50 - 2.00]"
    else:
        print("Error: no file input")

    df = pd.read_csv(file)
    target_label_frequency = df[targ].value_counts()
    hr_statistics_table = pd.DataFrame(target_label_frequency).reset_index()
    hr_statistics_table = hr_statistics_table.rename(
        {targ: "Label Frequency", "index": "Label"}, axis='columns')
    hr_statistics_table["Cancer Type"] = cancer_tissue
    hr_statistics_table["Prediction Target"] = targ
    hr_statistics_table["HR Thresholds"] = hr_threshold

    return freq
