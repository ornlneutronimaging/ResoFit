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


class Layer(object):
    def __init__(self):
        self.info = {}

    def add_layer(self, layer, thickness_mm, density_gcm3=None):
        if type(layer) is not str:
            raise ValueError("Please enter layer as string. Example: 'Gd' or 'U'")
        _formula = re.findall(r'([A-Z][a-z]*)(\d*)', layer)
        _elements = []
        for _element in _formula:
            _single_element = list(_element)[0]
            _elements.append(_single_element)
        if len(_elements) > 1:
            raise ValueError("Please enter element as layer in string. Example: 'Gd' or 'U'")

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

# def _forminfo(layer, thickness_mm, density_gcm3=None):
#     info = {}
#     if type(layer) is list:
#         if density_gcm3 is not None:
#             if type(density_gcm3) is not list:
#                 raise ValueError("Please enter density in (g/cm^3) for each layer you entered as list."
#                                  " Example: [8.6, 7.9]")
#         if type(thickness_mm) is not list:
#             raise ValueError("Please enter thickness in (mm) for each layer you entered as list in same length."
#                              " Example: [0.15, 0.3]")
#         if len(layer) != len(thickness_mm):
#             raise ValueError("Please enter thickness in (mm) for each layer you entered.")
#         for _layer_number in range(len(layer)):
#             print(_layer_number)
#             if density_gcm3 is not None:
#                 if density_gcm3[_layer_number] is not None:
#                     info[layer[_layer_number]] = {'layer': layer[_layer_number],
#                                                          'thickness': {'value': thickness_mm[_layer_number],
#                                                                        'units': 'mm'},
#                                                          'density': {'value': density_gcm3[_layer_number],
#                                                                      'units': 'g/cm3'}
#                                                          }
#
#                 else:
#                     info[layer[_layer_number]] = {'layer': layer[_layer_number],
#                                                          'thickness': {'value': thickness_mm[_layer_number],
#                                                                        'units': 'mm'},
#                                                          'density': {'value': np.NaN,
#                                                                      'units': 'g/cm3'}
#                                                          }
#     else:
#         if type(density_gcm3) is list:
#             raise ValueError("Density entered can not be a list for single layer '{}' you entered."
#                              " Density input example: 8.6 ".format(layer))
#         if type(thickness_mm) is list:
#             raise ValueError("Thickness entered can not be a list for single layer '{}' you entered."
#                              " Thickness input example: 0.5 ".format(layer))
#         if density_gcm3 is not None:
#             info[layer] = {'layer': layer,
#                                   'thickness': {'value': thickness_mm,
#                                                 'units': 'mm'},
#                                   'density': {'value': density_gcm3,
#                                               'units': 'g/cm3'}
#                                   }
#
#
#         else:
#             info[layer] = {'layer': layer,
#                                   'thickness': {'value': thickness_mm,
#                                                 'units': 'mm'},
#                                   'density': {'value': np.NaN,
#                                               'units': 'g/cm3'}
#                                   }
#
#     if info == {}:
#         raise ValueError("No info has been passed to form layer info stack.")
#
#     return info
