import re
import numpy as np
import pandas as pd
from ImagingReso import _utilities
from ImagingReso.resonance import Resonance
import os
from lmfit import Parameters


class Simulation(object):
    energy_min = np.NaN
    energy_max = np.NaN
    energy_step = np.NaN
    # Input sample name or names as str, case sensitive
    layer_1 = ''
    thickness_1 = np.NaN  # mm
    layers = []

    def __init__(self, layer_1, thickness_1, density_1=np.NaN,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01):

        self.o_reso = Resonance(energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
        self.o_reso.add_layer(formula=layer_1, thickness=thickness_1, density=density_1)
        self.layers.append(layer_1)

    def add_layer(self, layer_to_add, layer_thickness, layer_density=np.NaN):
        self.o_reso.add_layer(formula=layer_to_add,
                              thickness=layer_thickness,
                              density=layer_density)
        self.layers.append(layer_to_add)

    def set_isotopic_ratio(self, layer, element, new_isotopic_ratio_list=[]):
        # Check if layer exist
        if layer not in self.layers:
            raise ValueError('Layer {} does not exist.'.format(layer))
        _formula = re.findall(r'([A-Z][a-z]*)(\d*)', layer)
        # Check if element exist
        _elements = []
        for _element in _formula:
            _single_element = list(_element)[0]
            _elements.append(_single_element)
        if element not in _elements:
            raise ValueError('Element {} specified does not exist in {} layer.'.format(element, layer))
        self.o_reso.set_isotopic_ratio(compound=layer, element=element, list_ratio=new_isotopic_ratio_list)

    def x(self, angstrom=False):
        _x = self.o_reso.total_signal['energy_eV']
        if angstrom is True:
            _x = _utilities.ev_to_angstroms(_x)
        return _x

    def y(self, transmission=False):
        if transmission is True:
            _y = self.o_reso.total_signal['transmission']
        else:
            _y = self.o_reso.total_signal['attenuation']
        return _y

    def xy_simu(self, angstrom=False, transmission=False):
        _x = self.o_reso.total_signal['energy_eV']
        if angstrom is True:
            _x = _utilities.ev_to_angstroms(_x)
        if transmission is True:
            _y = self.o_reso.total_signal['transmission']
        else:
            _y = self.o_reso.total_signal['attenuation']
        return _x, _y

    def to_csv(self, filename='simulation.csv', angstrom=False, transmission=False):
        """
        Output x and y values to .csv file
        :param filename:
        :param angstrom:
        :param transmission:
        :return: .csv file
        """
        _x = self.o_reso.total_signal['energy_eV']
        _x_tag = 'x (eV)'
        if angstrom is True:
            _x = _utilities.ev_to_angstroms(_x)
            _x_tag = 'x (\u212B)'
        if transmission is True:
            _y = self.o_reso.total_signal['transmission']
            _y_tag = 'y (transmission)'
        else:
            _y = self.o_reso.total_signal['attenuation']
            _y_tag = 'y (attenuation)'

        df = pd.DataFrame(_x, index=None)
        df.rename(columns={0: _x_tag}, inplace=True)
        df[_y_tag] = _y
        df.to_csv(filename)

    # def x_layer(self, layer, angstrom=False):
    #     _x = self.o_reso.total_signal[layer]['energy_eV']
    #     if angstrom is True:
    #         _x = _utilities.ev_to_angstroms(_x)
    #     # pprint.pprint(o_reso.stack_sigma)
    #     # pprint.pprint(o_reso)
    #     return _x
    #
    # def y_layer(self, layer, transmission=False):
    #     if transmission is True:
    #         _y = self.o_reso.total_signal[layer]['transmission']
    #     else:
    #         _y = self.o_reso.stack_signal[layer]['attenuation']
    #     return _y
