import pandas as pd
import numpy as np
import os
import pprint
import re


def load_txt_csv(path_to_file):
    """
    Load and format data from .txt or .csv files
    :param path_to_file:
    :return: pd.Dataframe
    """
    # Error for file format and existence
    _format = path_to_file[-4:]
    if _format not in ['.txt', '.csv']:
        raise ValueError("File must be in the format of '.txt' or '.csv'")
    if os.path.exists(path_to_file) is False:
        raise ValueError(
            "Can not locate file '{}' in '{}' ".format(os.path.basename(path_to_file), os.path.dirname(path_to_file)))

    _sep = ','
    df = pd.read_csv(path_to_file, sep=_sep, header=None)

    if type(df[0][0]) is str:
        # if the first element is still a str, use ',' to pd.read_csv
        if df[0][0].count('\t') != 0:
            _sep = '\t'
            df = pd.read_csv(path_to_file, sep=_sep, header=None)

    if type(df[0][0]) is str:
        # if the first element is still a str, skip the row of the 'X' 'Y' axis labels
        df = pd.read_csv(path_to_file, sep=_sep, header=None, skiprows=1)

    if list(df[0][:4]) == [1, 2, 3, 4]:
        df[0] = df[1]
        df.drop(df.columns[1], axis=1, inplace=True)
    return df


def get_foil_density_gcm3(length_mm, width_mm, thickness_mm, mass_g):
    """
    Get density from mass/(L*W*H)
    :param length_mm:
    :param width_mm:
    :param thickness_mm:
    :param mass_g:
    :return: density in g/cm^3
    """
    _mm3_to_cm3 = 0.001
    density_gcm3 = mass_g / (length_mm * width_mm * thickness_mm * _mm3_to_cm3)
    return density_gcm3


def restructure_input(name):
    if type(name) is not str:
        raise ValueError("'{}' entered is not a string.".format(name))
    if len(name) == 0:
        raise ValueError("'{}' entered has no length.".format(name))

    _path_of_input = []
    if any(str.isdigit(i) for i in name) is True:
        if str.isdigit(name[0]) is False:
            raise ValueError("Please format you isotope name {} in the form of '238-U'.".format(name))
        _parsed = re.findall(r'([A-Z][a-z]*)(\d*)', name)
        _element_name = _parsed[0][0]
        _path_of_input.append(_element_name)
        _path_of_input.append(_element_name)
        _path_of_input.append(name)
    else:
        if len(name) > 2:
            raise ValueError("'{}' entered is not a single element symbol.".format(name))
        if len(name) == 1:
            if name.isupper() is False:
                name = name.upper()
        if len(name) == 2:
            if name[0].isupper() and name[1].islower() is True:
                _path_of_input.append(name)
                _path_of_input.append(name)
            else:
                raise ValueError("'{}' entered is not a valid element symbol.".format(name))
    return _path_of_input


class Layer(object):
    def __init__(self):
        self.info = {}

    def add_layer(self, layer, thickness_mm, density_gcm3=None):
        # raise error if input is not string.
        if type(layer) is not str:
            raise ValueError("Please enter layer as string. Example: 'Gd' or 'U'")
        _formula = re.findall(r'([A-Z][a-z]*)(\d*)', layer)
        _elements = []
        for _element in _formula:
            _single_element = list(_element)[0]
            _elements.append(_single_element)
        # raise error if input is contains more than one element for single layer.
        if len(_elements) > 1:
            raise ValueError("Please enter single element as layer in string. Example: 'Gd' or 'U'")

        if density_gcm3 is not None:
            self.info[layer] = {'layer': layer,
                                'thickness': {'value': thickness_mm,
                                              'units': 'mm'},
                                'density': {'value': density_gcm3,
                                            'units': 'g/cm3'},
                                'molar_mass': {'value': None,
                                               'units': None},
                                'molar_conc': {'value': None,
                                               'units': None}
                                }
        else:
            self.info[layer] = {'layer': layer,
                                'thickness': {'value': thickness_mm,
                                              'units': 'mm'},
                                'density': {'value': np.NaN,
                                            'units': 'g/cm3'},
                                'molar_mass': {'value': None,
                                               'units': None},
                                'molar_conc': {'value': None,
                                               'units': None}
                                }

    def show(self):
        pprint.pprint(self.info)
