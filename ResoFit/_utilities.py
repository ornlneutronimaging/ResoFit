import re
import numpy as np
import pandas as pd
from ImagingReso import _utilities
from ImagingReso.resonance import Resonance
import os


class Simulation(object):
    _energy_min = np.NaN
    _energy_max = np.NaN
    _energy_step = np.NaN
    # Input sample name or names as str, case sensitive
    layer_1 = ''
    thickness_1 = np.NaN  # mm
    layers = []

    def __init__(self, layer_1, thickness_1, density_1=np.NaN,
                 _energy_min=1e-5, _energy_max=1000, _energy_step=0.01):
        o_reso = Resonance(energy_min=_energy_min, energy_max=_energy_max, energy_step=_energy_step)
        o_reso.add_layer(formula=layer_1, thickness=thickness_1, density=density_1)
        self.layers.append(layer_1)
        self.o_reso = o_reso

    def add_layer(self, layer_added, thickness, density=np.NaN):
        o_reso = self.o_reso
        o_reso.add_layer(formula=layer_added, thickness=thickness, density=density)
        self.layers.append(layer_added)
        self.o_reso = o_reso

    def set_isotopic_ratio(self, layer, element, new_isotopic_ratio_list=[]):
        # Check if layer exist
        layers = self.layers
        if layer not in layers:
            raise ValueError('Layer {} does not exist.'.format(layer))
        _formula = re.findall(r'([A-Z][a-z]*)(\d*)', layer)
        # Check if element exist
        _elements = []
        for _element in _formula:
            _single_element = list(_element)[0]
            _elements.append(_single_element)
        if element not in _elements:
            raise ValueError('Element {} specified does not exist in {} layer.'.format(element, layer))
        o_reso = self.o_reso
        o_reso.set_isotopic_ratio(compound=layer, element=element, list_ratio=new_isotopic_ratio_list)
        self.o_reso = o_reso

    def x(self, layer_name, element_name, angstrom=False):
        o_reso = self.o_reso
        _x = o_reso.stack_sigma[layer_name][element_name]['energy_eV']
        if angstrom is True:
            _x = _utilities.ev_to_angstroms(_x)
        return _x

    def y(self, layer_name, element_name, transmission=False):
        o_reso = self.o_reso
        if transmission is True:
            _y = o_reso.stack_signal[layer_name][element_name]['attenuation']
        else:
            _y = o_reso.stack_signal[layer_name][element_name]['transmission']
        return _y


class Experiment(object):
    # folder = ''
    # spectra = ''
    # data = ''
    delay_us = np.NaN
    source_to_detector_m = np.NaN
    repeat = np.int
    # why need to define these outside __init__

    def __init__(self, spectra, data, folder='data', repeat=1, delay_us=0, source_to_detector_m=16.12):
        _file_path = os.path.abspath(os.path.dirname(__file__))
        _folder_path = os.path.join(_file_path, folder)
        if os.path.isdir(_folder_path) is False:
            raise ValueError('Folder specified does not exist')
        # Spectra file
        self.spectra_path = os.path.join(_folder_path, spectra)
        if os.path.exists(self.spectra_path) is False:
            raise ValueError("Can not find spectra file '{}' in '/{}' folder".format(spectra, folder))
        path_to_spectra, spectra_format = os.path.splitext(self.spectra_path)
        if spectra_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        # Data file
        self.data_path = os.path.join(_folder_path, data)
        if os.path.exists(self.data_path) is False:
            raise ValueError("Can not find data file '{}' in '/{}' folder".format(data, folder))
        path_to_data, date_format = os.path.splitext(self.data_path)
        if date_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        # Force repeat be an int >=1
        if isinstance(repeat, int) is False:
            raise ValueError("Repeat value must be an integer!")
        if repeat < 1:
            raise ValueError("Repeat value must be an integer >= 1 !")

        self.source_to_detector_m = source_to_detector_m
        self.delay_us = delay_us
        self.repeat = repeat

    def x(self, angstrom=False):
        spectra = pd.read_csv(self.spectra_path, sep='\t', header=None)
        x_in_s = np.array(spectra[0])
        y_in_counts = np.array(spectra[1])
        delay_us = self.delay_us
        source_to_detector_m = self.source_to_detector_m
        if angstrom is True:
            x = _utilities.s_to_angstroms(x_in_s,
                                          delay_us=delay_us,
                                          source_to_detector_m=source_to_detector_m)
        else:
            x = _utilities.s_to_ev(x_in_s,
                                   delay_us=delay_us,
                                   source_to_detector_m=source_to_detector_m)
        return x

    def ob_y(self):
        spectra = pd.read_csv(self.spectra_path, sep='\t', header=None)
        y_in_counts = np.array(spectra[1])
        return y_in_counts

    def y(self, transmission=False):
        data = pd.read_csv(self.data_path, sep='\t', header=None)
        if np.array(data[0])[:3] == [1, 2, 3, 4]:
            y = np.array(data[1]) / self.repeat
        else:
            y = np.array(data[0]) / self.repeat
        if transmission is False:
            y = 1 - y
        return y

