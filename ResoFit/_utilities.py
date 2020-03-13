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
import matplotlib.pyplot as plt

x_type_list = ['energy', 'lambda', 'time', 'number']
y_type_list = ['transmission', 'attenuation']
t_unit_list = ['s', 'ms', 'us', 'ns']
peak_type_list = ['indexed', 'all']
# peak_type_list = ['indexed', 'all', 'none']
index_level_list = ['iso', 'ele']
peak_model_list = ['Gaussian', 'Lorentzian']


def check_if_in_list(name, name_list):
    if name not in name_list:
        raise ValueError("'{}' is not valid, only support: '{}'".format(name, name_list))


def convert_energy_to(x_type, x,
                      offset_us=None, source_to_detector_m=None, t_unit='us',
                      num_offset=0,
                      time_resolution_us=None,
                      t_start_us=None):
    check_if_in_list(x_type, x_type_list)
    check_if_in_list(t_unit, t_unit_list)
    if x_type == 'lambda':
        x = reso_util.ev_to_angstroms(x)
    if x_type == 'time':
        if offset_us is None:
            raise ValueError("'offset_us=' is required when x_type='time'")
        if source_to_detector_m is None:
            raise ValueError("'source_to_detector_m=' is required when x_type='time'")
        x = reso_util.ev_to_s(offset_us=offset_us,
                              source_to_detector_m=source_to_detector_m,
                              array=x)
        x = convert_s(x=x, t_unit=t_unit)
    if x_type == 'number':
        if time_resolution_us is not None:
            x = reso_util.ev_to_image_number(offset_us=offset_us,
                                             source_to_detector_m=source_to_detector_m,
                                             array=x,
                                             time_resolution_us=time_resolution_us,
                                             t_start_us=t_start_us)
        else:
            x = np.array(range(len(x))) + num_offset
    return x


def get_peak_tag(x_type):
    tag = 'x'
    if x_type == 'lambda':
        tag = 'x_A'
    if x_type == 'time':
        tag = 'x_s'
    if x_type == 'number':
        tag = 'x_num_o'
    return tag


def get_df_col_name(x_type):
    tag = '_E'
    if x_type == 'lambda':
        tag = '_A'
    if x_type == 'time':
        tag = '_t'
    if x_type == 'number':
        tag = '_#'
    return tag


def convert_attenuation_to(y_type, y):
    check_if_in_list(y_type, y_type_list)
    if y_type == 'transmission':
        y = 1 - y
    return np.array(y)


def convert_s(x, t_unit):
    if t_unit == 'ns':
        _x = x * 1e9
    elif t_unit == 'us':
        _x = x * 1e6
    elif t_unit == 'ms':
        _x = x * 1e3
    else:
        _x = x
    return _x


def convert_exp_peak_df(peak_df: pd.DataFrame, x_type, t_unit):
    check_if_in_list(x_type, x_type_list)
    check_if_in_list(t_unit, t_unit_list)
    if x_type == 'energy':
        assert 'x' in peak_df.columns
        _x = peak_df['x']
    elif x_type == 'lambda':
        assert 'x_A' in peak_df.columns
        _x = peak_df['x_A']
    elif x_type == 'time':
        assert 'x_s' in peak_df.columns
        _x = peak_df['x_s']
        _x = convert_s(x=_x, t_unit=t_unit)
    else:
        assert 'x_num_o' in peak_df.columns
        _x = peak_df['x_num_o']
    return _x.values  # np.array


def check_and_make_dir(current_path, name):
    _dir_path = os.path.join(current_path, name)

    if not os.path.exists(_dir_path):
        os.makedirs(_dir_path)
        print("Folder: '{}' has been created ".format(_dir_path))

    return _dir_path


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


def set_plt(ax, fig_title, grid, x_type, y_type, t_unit, logx, logy):
    check_if_in_list(x_type, x_type_list)
    check_if_in_list(y_type, y_type_list)
    ax.set_title(fig_title)
    if x_type == 'energy':
        ax.set_xlabel('Energy (eV)')
    elif x_type == 'lambda':
        ax.set_xlabel('Wavelength (\u212B)')
    elif x_type == 'number':
        ax.set_xlabel('Image number (#)')
    else:
        check_if_in_list(t_unit, t_unit_list)
        if t_unit == 'us':
            ax.set_xlabel('Time of flight (\u03BCs)')
        else:
            ax.set_xlabel('Time of flight ({})'.format(t_unit))

    if y_type == 'attenuation':
        ax.set_ylabel('Neutron attenuation')
    else:
        ax.set_ylabel('Neutron transmission')
    ax.legend(loc='best')
    # ax1.legend(bbox_to_anchor=(1., 1), loc=2, borderaxespad=0.)
    # ax1.legend(bbox_to_anchor=(0, 0.93, 1., .102), loc=3, borderaxespad=0.)
    if grid:
        # ax1.set_xticks(np.arange(0, 100, 10))
        # ax1.set_yticks(np.arange(0, 1., 0.1))
        ax.grid()
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')
    return ax


def rm_envelope(y, deg, max_it=None, tol=None):
    envelope = pku.envelope(y=y, deg=deg, max_it=max_it, tol=tol)
    # return y + y.max() - envelope
    return y / envelope


class Items(object):
    """
    A easier way to specify layers/elements/isotopes for in plot()/export()

    """

    def __init__(self, o_reso, database='ENDF_VIII'):
        self.o_reso = o_reso
        self.shaped_list = None
        self.database = database

    def shaped(self, items_list):
        _shaped_list = []
        for _raw_path_to_plot in items_list:
            if type(_raw_path_to_plot) is not list:
                if '*' in _raw_path_to_plot:
                    _shaped_list = _shaped_list + _fill_iso_to_items(name=_raw_path_to_plot,
                                                                     stack=self.o_reso.stack,
                                                                     database=self.database)
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
        if self.shaped_list is None:
            raise ValueError("'.shaped_list' is empty, please run '.shaped()' first.")
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


def _fill_iso_to_items(name, stack=None, database='ENDF_VII'):
    if '*' not in name:
        raise ValueError("'*' is needed to retrieve all isotopes of '{}' ".format(name))
    else:
        ele_name = name.replace('*', '')
        if stack is None:
            o_reso = Resonance(database=database)
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

    def add_Layer(self, layers):
        for _each_layer in list(layers.info.keys()):
            self.add_layer(layer=layers.info[_each_layer]['layer'],
                           thickness_mm=layers.info[_each_layer]['thickness'],
                           density_gcm3=layers.info[_each_layer]['density'])

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


def find_peak(y, x=None, x_name='x_num', y_name='y', thres=0.015, min_dist=1, imprv_reso=False):
    if x is None:
        x = np.array(range(len(y)))
    # x_num_gap = 0
    # Note: weirdly, indexes have to be reset here to get correct peak locations
    x = np.array(x)
    y = np.array(y)
    _index = pku.indexes(y=y, thres=thres, min_dist=min_dist)
    if len(_index) != 0:
        _peak_y = list(y[_index])
        if imprv_reso is False:
            _peak_x = list(x[_index])
        else:
            _peak_x = list(pku.interpolate(x, y, ind=_index))
    else:
        # No peaks detected
        _peak_y = []
        _peak_x = []

    peak_df = pd.DataFrame()
    peak_df[y_name] = _peak_y
    # peak_df[x_name] = [x_num + x_num_gap for x_num in _peak_x]
    peak_df[x_name] = _peak_x
    peak_df.sort_values([x_name], inplace=True)
    peak_df.reset_index(inplace=True, drop=True)
    return peak_df


def index_peak(peak_df, peak_map, rel_tol, x_name='x'):
    num_peak_indexed = 0
    _names = peak_map.keys()
    peak_map_indexed = {}
    for _peak_name in _names:
        _df = pd.DataFrame()
        _df_ideal = pd.DataFrame()
        peak_map_indexed[_peak_name] = {}
        _peak_x = peak_map[_peak_name]['ideal']['x']
        _peak_y = peak_map[_peak_name]['ideal']['y']
        _x_indexed_list = []
        _x_num_indexed_list = []
        _x_num_o_indexed_list = []
        _x_s_indexed_list = []
        _x_A_indexed_list = []
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
                    if 'x_num_o' in peak_df.columns:
                        _x_num_o_indexed_list.append(peak_df['x_num_o'][_i])
                    if 'x_s' in peak_df.columns:
                        _x_s_indexed_list.append(peak_df['x_s'][_i])
                    if 'x_A' in peak_df.columns:
                        _x_A_indexed_list.append(peak_df['x_A'][_i])
                        # _peak_name_list.append(_peak_name)
        num_peak_indexed += len(_x_indexed_list)
        _df[x_name] = _x_indexed_list
        _df['y'] = _y_indexed_list
        _df_ideal['x'] = _x_ideal_list
        _df_ideal['y'] = _y_ideal_list
        if 'x_num' in peak_df.columns:
            _df['x_num'] = _x_num_indexed_list
            _df_ideal['x_num'] = _x_num_indexed_list
        if 'x_num_o' in peak_df.columns:
            _df['x_num_o'] = _x_num_o_indexed_list
            _df_ideal['x_num_o'] = _x_num_o_indexed_list
        if 'x_s' in peak_df.columns:
            _df['x_s'] = _x_s_indexed_list
            _df_ideal['x_s'] = _x_s_indexed_list
        if 'x_A' in peak_df.columns:
            _df['x_A'] = _x_A_indexed_list
            _df_ideal['x_A'] = reso_util.ev_to_angstroms(np.array(_x_ideal_list))

        peak_map_indexed[_peak_name]['exp'] = _df
        peak_map_indexed[_peak_name]['ideal'] = _df_ideal
    return peak_map_indexed


class ResoPeak(object):
    def __init__(self, y, x, y_type, x_type):
        """
        Initialization

        """
        self.peak_dict = {}

        self.y = y
        self.x = x
        self.y_type = y_type
        self.x_type = x_type

        self.shape_report = None
        self.prefix_list = None

        self.x_num_gap = 0

    def find_peak(self, thres, min_dist, imprv_reso: bool):
        _peak_dict = self._find_peak(y=self.y, x=self.x, thres=thres, min_dist=min_dist, imprv_reso=imprv_reso)
        _peak_dict['x_type'] = self.x_type
        _peak_dict['y_type'] = self.y_type
        self.peak_dict = _peak_dict
        return _peak_dict

    def _find_peak(self, y: np.array, thres, min_dist, imprv_reso: bool, x=None):
        """"""
        if x is None:
            x = np.array(range(len(y)))
        else:
            x = np.array(x)
            if x.shape != y.shape:
                raise ValueError("The length ({}) of 'x' is not equal the length ({}) of 'y'".format(len(x), len(y)))
        peak_index = pku.indexes(y=y, thres=thres, min_dist=min_dist)
        if len(peak_index) != 0:
            _peak_y = y[peak_index]
            if imprv_reso:
                _peak_x = pku.interpolate(x, y, ind=peak_index)
            else:
                _peak_x = x[peak_index]
        else:
            # No peaks detected
            _peak_x = []
            _peak_y = []

        peak_dict = {
            'x': _peak_x,
            'y': _peak_y,
        }

        return peak_dict

    # def scale_with_ev(self, energy_min, energy_max):
    #     assert self.peak_dict != {}
    #     assert self.

    def convert_peak_dict_type(self, x_type_target, y_type_target):
        _peak_dict = self.peak_dict



class Peak(object):
    def __init__(self):
        """
        Initialization

        """
        self.peak_df = None
        self.peak_df_scaled = None
        self.peak_map_full = None
        self.peak_map_indexed = None
        self.y = None
        self.x = None

        self.shape_report = None
        self.prefix_list = None

        self.x_num_gap = 0

    def find(self, x, y, y_name='y', thres=0.015, min_dist=1, impr_reso=False):
        """
        find peaks in 1d data and return peaks detected using pd.DataFrame
        :param y: 1d data
        :type y: array
        :param x:
        :type x: array (optional)
        :param y_name:
        :type y_name:
        :param thres:
        :type thres:
        :param min_dist:
        :type min_dist:
        :param impr_reso:
        :type impr_reso:
        :return:
        :rtype:
        """
        if x is None:
            x = np.array(range(len(y)))
        self.x = x
        self.y = y
        peak_df = find_peak(y=y, x=None, x_name='x_num', y_name=y_name,
                            thres=thres, min_dist=min_dist, imprv_reso=impr_reso)
        _peak_df = find_peak(y=y, x=x, x_name='x_s', y_name=y_name,
                             thres=thres, min_dist=min_dist, imprv_reso=impr_reso)
        if type(y) == pd.core.series.Series:
            if y.index[0] != 1:
                self.x_num_gap = y.index[0]
                peak_df['x_num_o'] = peak_df['x_num'] + self.x_num_gap
        peak_df['x_s'] = _peak_df['x_s']
        self.peak_df = peak_df
        return self.peak_df

    # def add_x_s(self, y, x=None, x_name='x', y_name='y', thres=0.15, min_dist=1, impr_reso=False):
    #     if x_name == 'x_s':
    #         self.x_s = x
    #
    #     _peak_df = find_peak(y=y, x=x, x_name=x_name, y_name=y_name,
    #                          thres=thres, min_dist=min_dist, impr_reso=impr_reso)
    #     _x_name = x_name
    #     if _x_name in self.peak_df.columns:
    #         _x_name += '0'
    #     _len_before = len(self.peak_df.columns)
    #     self.peak_df[_x_name] = _peak_df[x_name]
    #     _len_after = len(self.peak_df.columns)
    #     assert _len_after > _len_before

    def _extend_x_cols(self, offset_us, source_to_detector_m):
        assert 'x_s' in self.peak_df.keys()
        _peak_df = self.peak_df.copy()
        _peak_df['x'] = reso_util.s_to_ev(array=_peak_df['x_s'],
                                          offset_us=offset_us,
                                          source_to_detector_m=source_to_detector_m)
        _peak_df['x_A'] = reso_util.ev_to_angstroms(_peak_df['x'])
        self.peak_df = _peak_df

    def _scale_peak_df(self, energy_min, energy_max):
        _peak_df_scaled = self.peak_df.copy()
        _peak_df_scaled.drop(_peak_df_scaled[_peak_df_scaled.x < energy_min].index, inplace=True)
        _peak_df_scaled.drop(_peak_df_scaled[_peak_df_scaled.x > energy_max].index, inplace=True)
        _peak_df_scaled.reset_index(drop=True, inplace=True)
        self.peak_df_scaled = _peak_df_scaled

    def index_peak(self, peak_map, rel_tol=5e-3):
        if self.peak_df is None:
            raise ValueError("Please identify use 'Peak.find()' before indexing peak.")
        assert self.peak_df_scaled is not None
        self.peak_map_indexed = index_peak(peak_df=self.peak_df_scaled, peak_map=peak_map, rel_tol=rel_tol)

    def analyze(self, report=False, fit_model='Lorentzian'):
        check_if_in_list(fit_model, peak_model_list)
        _y = self.y
        _x = self.x
        _peak_map_indexed = self.peak_map_indexed
        model = lmfit.models.GaussianModel(prefix='bkg_')
        pars = model.guess(_y, x=_x)
        self.prefix_list = []
        for _ele in _peak_map_indexed.keys():
            if '-' not in _ele:
                for _ind in range(len(_peak_map_indexed[_ele]['exp'])):
                    _prefix = _ele + '_' + str(_ind) + '_'
                    if fit_model == 'Gaussian':
                        _model = lmfit.models.GaussianModel(prefix=_prefix)
                    else:  # fit_model == 'Lorentzian':
                        _model = lmfit.models.LorentzianModel(prefix=_prefix)
                    _center = _peak_map_indexed[_ele]['exp']['x_num'][_ind]
                    pars.update(_model.make_params())
                    pars[_prefix + 'amplitude'].value = 3.0
                    pars[_prefix + 'center'].set(_center, min=_center - 10, max=_center + 10)
                    pars[_prefix + 'sigma'].set(2.0, min=0.5, max=50)
                    model += _model
                    self.prefix_list.append(_prefix)
        _out = model.fit(_y, pars, x=_x)
        self.shape_report = _out
        self.__fwhm()
        self.__fill_img_num_to_peak_map_indexed()
        print("+------------ Peak analysis ------------+\n{} peak fitting:".format(fit_model))
        print("{}\n".format(self.fwhm_df))

        if report is True:
            print(_out.fit_report())

    def plot_fit(self):
        if self.shape_report is not None:
            self.shape_report.plot()
            plt.show()
        else:
            print("Peaks have not been fitted. Please run 'Peak.analyze()' before plotting.")

    def __fwhm(self):
        _fwhm_df = pd.DataFrame()
        # generate ele list for _fwhm_df
        _ele_list = [_ele_name.split('_')[0] for _ele_name in self.prefix_list]
        _prefix_list = self.prefix_list
        _values = self.shape_report.__dict__['params'].valuesdict()
        pars_center_name = [_i + 'center' for _i in _prefix_list]
        pars_fwhm_name = [_i + 'fwhm' for _i in _prefix_list]
        pars_center_value = [_values[_name] for _name in pars_center_name]
        pars_fwhm_value = [_values[_name] for _name in pars_fwhm_name]
        _fwhm_df['ele_name'] = _ele_list
        _fwhm_df['center_val'] = pars_center_value
        _fwhm_df['fwhm_val'] = pars_fwhm_value
        _fwhm_df.sort_values(['center_val'], inplace=True)
        _fwhm_df.reset_index(inplace=True, drop=True)
        self.fwhm_df = _fwhm_df

    def __fill_img_num_to_peak_map_indexed(self):
        assert 'x_s' in self.peak_df.keys()
        # if self.x_s is None:
        #     raise ValueError("Column of x in time (s) has not been added.")
        _peak_map_indexed = self.peak_map_indexed
        _fwhm_df = self.fwhm_df
        for _ele in _peak_map_indexed.keys():
            _peak_map_indexed[_ele]['peak_span'] = {}
            _img_num_list = []
            _peak_span_df = pd.DataFrame()
            for _ind in range(len(_fwhm_df)):
                if _fwhm_df['ele_name'][_ind] == _ele:
                    half_fwhm = _fwhm_df['fwhm_val'][_ind] / 2
                    _min = _fwhm_df['center_val'][_ind] - half_fwhm + self.x_num_gap
                    _max = _fwhm_df['center_val'][_ind] + half_fwhm + self.x_num_gap
                    _min = int(np.floor(_min))
                    _max = int(np.ceil(_max)) + 1
                    _img_num_list += [a for a in range(_min, _max)]
            _peak_span_df['img_num'] = _img_num_list
            _peak_map_indexed[_ele]['peak_span'] = _peak_span_df
        self.peak_map_indexed = _peak_map_indexed

    def fill_peak_span(self, offset_us, source_to_detector_m):
        assert 'x_s' in self.peak_df.keys()
        for _keys in self.peak_map_indexed.keys():
            _live_location = self.peak_map_indexed[_keys]['peak_span']
            _img_num_list = _live_location['img_num']
            # _live_location['time_s'] = list(self.x_s.reindex(_img_num_list))
            _live_location['time_s'] = list(self.peak_df['x_s'].reindex(_img_num_list))
            # _live_location['energy_ev'] = reso_util.s_to_ev(array=list(self.x_s.reindex(_img_num_list)),
            _live_location['energy_ev'] = reso_util.s_to_ev(array=_live_location['time_s'],
                                                            offset_us=offset_us,
                                                            source_to_detector_m=source_to_detector_m)
            _live_location['y'] = list(self.y.reindex(_img_num_list))
            # _live_location['y'] = list(self.y.loc[_img_num_list])

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
