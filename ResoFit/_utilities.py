import pandas as pd


def load_txt_csv(path_to_file):
    _sep = ','
    df = pd.read_csv(path_to_file, sep=_sep, header=None)
    if type(df[0][0]) is str:
        if df[0][0].count('\t') != 0:
            _sep = '\t'
            df = pd.read_csv(path_to_file, sep=_sep, header=None)
            #
            # if type(df[0][0]) is str:
            #     if df[0][0].islower() or df[0][0].isupper() is True:
            #         df = pd.read_csv(path_to_file, sep=_sep, header=None, skiprows=1)

    # for _i in range(6):
    if type(df[0][0]) is str:
            # if df[0][0].islower() or df[0][0].isupper() is True:
        df = pd.read_csv(path_to_file, sep=_sep, header=None, skiprows=1)
    return df


