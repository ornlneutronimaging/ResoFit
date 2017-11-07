import itertools
import os
import pprint
import re
from math import isclose

import lmfit
import numpy as np
import pandas as pd
import peakutils as pku
from ImagingReso.resonance import Resonance
import ImagingReso._utilities as reso_util
from cerberus import Validator


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


def set_plt(plt, x_max, fig_title, grid=False):
    plt.set_xlim([0, x_max])
    plt.set_ylim(ymax=1.01)
    plt.set_title(fig_title)
    plt.set_xlabel('Energy (eV)')
    plt.set_ylabel('Neutron attenuation')
    plt.legend(loc='best')
    # ax1.legend(bbox_to_anchor=(1., 1), loc=2, borderaxespad=0.)
    # ax1.legend(bbox_to_anchor=(0, 0.93, 1., .102), loc=3, borderaxespad=0.)
    if grid is True:
        # ax1.set_xticks(np.arange(0, 100, 10))
        # ax1.set_yticks(np.arange(0, 1., 0.1))
        plt.grid()


# def peak_plt(plt, peak_df, peak_map_indexed, peak):
#     plt.plot(self.peak_df_scaled['x'],
#                      self.peak_df_scaled['y'],
#                      'kx', label='_nolegend_')


def rm_baseline(y, deg=7, max_it=None, tol=None):
    baseline = pku.baseline(y=y, deg=deg, max_it=max_it, tol=tol)
    return y - baseline


class Items(object):
    def __init__(self, o_reso):
        self.o_reso = o_reso
        self.shaped_list = None

    def shaped(self, items_list):
        _shaped_list = []
        for _raw_path_to_plot in items_list:
            if type(_raw_path_to_plot) is not list:
                if '*' in _raw_path_to_plot:
                    _shaped_list = _shaped_list + _fill_iso_to_items(_raw_path_to_plot, self.o_reso.stack)
                else:
                    _shaped_list.append(_shape_items(_raw_path_to_plot))
            else:
                if len(_raw_path_to_plot) == 1:
                    _raw_path_to_plot = _shape_items(_raw_path_to_plot[0])
                _shaped_list.append(_raw_path_to_plot)
        # Clean duplicates in list
        _shaped_list = _rm_duplicated_items(_shaped_list)
        self.shaped_list = _shaped_list
        return _shaped_list

    def values(self, y_axis_type='attenuation'):
        # plot specified from 'items_to_plot'
        if y_axis_type != 'sigma':
            _stack = self.o_reso.stack_signal
        else:
            _stack = self.o_reso.stack_sigma
            y_axis_type = 'sigma_b'
        y_axis_tag = y_axis_type
        _y_axis_dict = {}
        for _each_path in self.shaped_list:
            _label = _each_path[-1]
            if len(_each_path) == 3:
                _y_axis_dict[_label] = _stack[_each_path[0]][_each_path[1]][_each_path[2]][y_axis_tag]
            elif len(_each_path) == 2:
                _y_axis_dict[_label] = _stack[_each_path[0]][_each_path[1]][y_axis_tag]
            else:
                raise ValueError("Format error of '{}', should be in the form of "
                                 "['layer', 'element'] or ['layer', 'element', 'isotope']")
        return _y_axis_dict


def _shape_items(name):
    # input is not structured as required by ImagingReso
    if type(name) is not str:
        raise ValueError("'{}' entered is not a string.".format(name))
    if len(name) == 0:
        raise ValueError("'{}' entered has no length.".format(name))
    _path_of_input = []

    if any(str.isdigit(i) for i in name) is True:
        # isotopes
        _parsed = re.findall(r'([A-Z][a-z]*)(\d*)', name)
        _element_str = _parsed[0][0]
        _number_str = re.findall('\d+', name)[0]
        _isotope_str = _number_str + '-' + _element_str
        _path_of_input.append(_element_str)
        _path_of_input.append(_element_str)
        _path_of_input.append(_isotope_str)
    else:
        # elements
        if len(name) > 2:
            raise ValueError("'{}' entered is not a single element symbol.".format(name))
        if len(name) == 1:
            if name.isupper() is False:
                name = name.upper()
            _path_of_input.append(name)
            _path_of_input.append(name)
        if len(name) == 2:
            if name[0].isupper() and name[1].islower() is True:
                _path_of_input.append(name)
                _path_of_input.append(name)
            else:
                raise ValueError("'{}' entered is not a valid element symbol.".format(name))
    return _path_of_input


def _fill_iso_to_items(name, stack=None):
    if '*' not in name:
        raise ValueError("'*' is needed to retrieve all isotopes of '{}' ".format(name))
    else:
        ele_name = name.replace('*', '')
        if stack is None:
            o_reso = Resonance()
            o_reso.add_layer(formula=ele_name,
                             thickness=1)
            stack = o_reso.stack
        iso_list = stack[ele_name][ele_name]['isotopes']['list']
        _path_to_iso = []
        for _each_iso in iso_list:
            _path_to_iso.append(_shape_items(_each_iso))
    return _path_to_iso


def _rm_duplicated_items(raw):
    raw.sort()
    cleaned_list = list(raw for raw, _ in itertools.groupby(raw))
    return cleaned_list


# def almostequatl


class Layer(object):
    def __init__(self):
        self.info = {}

    def add_layer(self, layer, thickness_mm, density_gcm3=None):

        # Input Validation
        _input = {'layer': layer,
                  'thickness': thickness_mm,
                  'density': density_gcm3,
                  }

        schema = {'layer': {'type': 'string',
                            'required': True,
                            },
                  'thickness': {'type': 'number',
                                'required': True,
                                },
                  'density': {'type': 'number',
                              'required': True,
                              'nullable': True,
                              },
                  }

        v = Validator(schema)
        if v.validate(_input) is False:
            raise ValueError(v.errors)

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
                                              'units': 'mm',
                                              },
                                'density': {'value': density_gcm3,
                                            'units': 'g/cm3',
                                            },
                                'molar_mass': {'value': None,
                                               'units': None,
                                               },
                                'molar_conc': {'value': None,
                                               'units': None,
                                               },
                                }
        else:
            self.info[layer] = {'layer': layer,
                                'thickness': {'value': thickness_mm,
                                              'units': 'mm',
                                              },
                                'density': {'value': np.NaN,
                                            'units': 'g/cm3',
                                            },
                                'molar_mass': {'value': None,
                                               'units': None,
                                               },
                                'molar_conc': {'value': None,
                                               'units': None,
                                               },
                                }

    def pprint(self):
        pprint.pprint(self.info)


def find_peak(y, x=None, x_name='x', y_name='y', thres=0.015, min_dist=1, impr_reso=False):
    if x is None:
        x = np.array(range(0, len(y)))
    _index = pku.indexes(y=y, thres=thres, min_dist=min_dist)
    _peak_y = list(y[_index])
    if impr_reso is False:
        _peak_x = list(x[_index])
    else:
        _peak_x = list(pku.interpolate(x, y, ind=_index))
    peak_df = pd.DataFrame()
    peak_df[y_name] = _peak_y
    peak_df[x_name] = _peak_x
    peak_df.sort_values([x_name], inplace=True)
    peak_df.reset_index(inplace=True, drop=True)
    return peak_df


def index_peak(peak_df, peak_map, x_name='x', rel_tol=3.5e-3):
    # if type(peak_map) == dict is False:
    num_peak_detected = len(peak_df[x_name])
    num_peak_indexed = 0
    peak_indexed = []
    _names = peak_map.keys()
    peak_map_indexed = {}
    for _peak_name in _names:
        _df = pd.DataFrame()
        _df_ideal = pd.DataFrame()
        peak_map_indexed[_peak_name] = {}
        _peak_x = peak_map[_peak_name]['peak']['x']
        _peak_y = peak_map[_peak_name]['peak']['y']
        _x_indexed_list = []
        _x_num_indexed_list = []
        _x_s_indexed_list = []
        _y_indexed_list = []
        _x_ideal_list = []
        _y_ideal_list = []
        for _i in range(len(peak_df[x_name])):
            for _j in range(len(_peak_x)):
                if isclose(_peak_x[_j], peak_df[x_name][_i], rel_tol=rel_tol):
                    _y_indexed_list.append(peak_df['y'][_i])
                    _x_indexed_list.append(peak_df[x_name][_i])
                    _x_ideal_list.append(_peak_x[_j])
                    _y_ideal_list.append(_peak_y[_j])
                    if 'x_num' in peak_df.columns:
                        _x_num_indexed_list.append(peak_df['x_num'][_i])
                    if 'x_s' in peak_df.columns:
                        _x_s_indexed_list.append(peak_df['x_s'][_i])
        num_peak_indexed += len(_x_indexed_list)
        _df[x_name] = _x_indexed_list
        _df['y'] = _y_indexed_list
        if 'x_num' in peak_df.columns:
            _df['x_num'] = _x_num_indexed_list
        if 'x_s' in peak_df.columns:
            _df['x_s'] = _x_s_indexed_list

        _df_ideal['x'] = _x_ideal_list
        _df_ideal['y'] = _y_ideal_list
        peak_map_indexed[_peak_name]['exp'] = _df
        peak_map_indexed[_peak_name]['ideal'] = _df_ideal
    # pprint.pprint(peak_map_indexed)
    return peak_map_indexed


class Peak(object):
    def __init__(self):
        self.thres = None
        self.min_dist = None
        self.impr_reso = None

        self.peak_df = None
        self.peak_df_scaled = None
        self.peak_map_full = None
        self.peak_map_indexed = None

    def find(self, y, x=None, x_name='x', y_name='y', thres=0.015, min_dist=1, impr_reso=False):
        self.thres = thres
        self.min_dist = min_dist
        self.impr_reso = impr_reso
        self.peak_df = find_peak(y=y, x=x, x_name=x_name, y_name=y_name,
                                 thres=thres, min_dist=min_dist, impr_reso=impr_reso)

    def add_x_col(self, y, x=None, x_name='x', y_name='y', thres=0.015, min_dist=1, impr_reso=False):
        _peak_df = find_peak(y=y, x=x, x_name=x_name, y_name=y_name,
                             thres=thres, min_dist=min_dist, impr_reso=impr_reso)
        _x_name = x_name
        if _x_name in self.peak_df.columns:
            _x_name += '0'
        _len_before = len(self.peak_df.columns)
        self.peak_df[_x_name] = _peak_df[x_name]
        _len_after = len(self.peak_df.columns)
        assert _len_after > _len_before

    def add_ev_and_scale(self, calibrated_source_to_detector_m, calibrated_offset_us,
                         energy_min, energy_max,
                         _from='x_s', _to='x'):
        self.peak_df[_to] = reso_util.s_to_ev(array=self.peak_df[_from],
                                              source_to_detector_m=calibrated_source_to_detector_m,
                                              offset_us=calibrated_offset_us)
        _peak_df_scaled = self.peak_df
        _peak_df_scaled.drop(_peak_df_scaled[_peak_df_scaled.x < energy_min].index, inplace=True)
        _peak_df_scaled.drop(_peak_df_scaled[_peak_df_scaled.x > energy_max].index, inplace=True)
        _peak_df_scaled.reset_index(drop=True, inplace=True)
        self.peak_df_scaled = _peak_df_scaled

    def index(self, peak_map, rel_tol=5e-3):
        if self.peak_df is None:
            raise ValueError("Please identify use 'Peak.find()' before indexing peak.")
        assert self.peak_df_scaled is not None
        self.peak_map_indexed = index_peak(peak_df=self.peak_df_scaled, peak_map=peak_map, rel_tol=rel_tol)

    def analyze(self):
        model = lmfit.models.GaussianModel()
        pass

# def a_new_decorator(a_func):
#     @wraps(a_func)
#     def wrapTheFunction():
#         print("I am doing some boring work before executing a_func()")
#         a_func()
#         print("I am doing some boring work after executing a_func()")
#
#     return wrapTheFunction
#
#
# @a_new_decorator
# def a_function_requiring_decoration():
#     """Hey yo! Decorate me!"""
#     print("I am the function which needs some decoration to "
#           "remove my foul smell")
#
#
# class Plot(object):
#     def __init__(self, logfile='out.log'):
#         self.logfile = logfile
#
#     def __call__(self, func):
#         log_string = func.__name__ + " was called"
#         print(log_string)
#         # Open the logfile and append
#         with open(self.logfile, 'a') as opened_file:
#             # Now we log to the specified logfile
#             opened_file.write(log_string + '\n')
#         # Now, send a notification
#         self.notify()
#
#     def notify(self):
#         # logit only logs, no more
#         pass
#
#
# class Export(object):
#     def __init__(self, logfile='out.log'):
#         self.logfile = logfile
#
#     def __call__(self, func):
#         log_string = func.__name__ + " was called"
#         print(log_string)
#         # Open the logfile and append
#         with open(self.logfile, 'a') as opened_file:
#             # Now we log to the specified logfile
#             opened_file.write(log_string + '\n')
#         # Now, send a notification
#         self.notify()
#
#     def notify(self):
#         # logit only logs, no more
#         pass
#
#
# class Logit(object):
#     def __init__(self, logfile='out.log'):
#         self.logfile = logfile
#
#     def __call__(self, func):
#         log_string = func.__name__ + " was called"
#         print(log_string)
#         # Open the logfile and append
#         with open(self.logfile, 'a') as opened_file:
#             # Now we log to the specified logfile
#             opened_file.write(log_string + '\n')
#         # Now, send a notification
#         self.notify()
#
#     def notify(self):
#         # logit only logs, no more
#         pass
