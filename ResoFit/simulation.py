import re
import numpy as np
import pandas as pd
import ImagingReso._utilities as reso_util
from ImagingReso.resonance import Resonance
import ResoFit._utilities as fit_util


class Simulation(object):

    # Input sample name or names as str, case sensitive

    def __init__(self, energy_min=1e-5, energy_max=1000, energy_step=0.01):
        """
        initialize the a Simulation() using the Resonance() in ImagingReso

        :param energy_min:
        :param energy_max:
        :param energy_step:
        """
        self.nergy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step

        self.o_reso = Resonance(energy_min=energy_min, energy_max=energy_max, energy_step=energy_step)
        # self.o_reso.add_layer(formula=layer_1, thickness=thickness_1, density=density_1)
        # self.layer_1 = layer_1
        # self.layers.append(layer_1)
        # self.simu_x = self.o_reso.total_signal['energy_eV']
        # self.simu_y = self.o_reso.total_signal['attenuation']
        self.simu_x = None
        self.simu_y = None
        self.layer_list = []

    def add_layer(self, layer, layer_thickness_mm, layer_density_gcm3=np.NaN):
        """
        Add layers and update x y values to pass
        :param layer:
        :param layer_thickness_mm:
        :param layer_density_gcm3: can be omitted same as Resonance() in ImagingReso
        :return: x in eV
                 y in attenuation
        """
        self.o_reso.add_layer(formula=layer,
                              thickness=layer_thickness_mm,
                              density=layer_density_gcm3)
        self.layer_list.append(layer)
        self.simu_x = self.o_reso.total_signal['energy_eV']
        self.simu_y = self.o_reso.total_signal['attenuation']

    def set_isotopic_ratio(self, layer, element, new_isotopic_ratio_list):
        """
        Set isotopic ratios for picked element and update x y values to pass
        :param layer:
        :param element:
        :param new_isotopic_ratio_list:
        :return: x in eV
                 y in attenuation
        """
        if type(new_isotopic_ratio_list) is not list:
            raise ValueError("{} is not a list".format(new_isotopic_ratio_list))
        # Check if layer exist
        if layer not in self.layer_list:
            raise ValueError('Layer {} does not exist.'.format(layer))
        # Check if element exist
        _formula = re.findall(r'([A-Z][a-z]*)(\d*)', layer)
        _elements = []
        for _element in _formula:
            _single_element = list(_element)[0]
            _elements.append(_single_element)
        if element not in _elements:
            raise ValueError('Element {} specified does not exist in {} layer.'.format(element, layer))
        self.o_reso.set_isotopic_ratio(compound=layer, element=element, list_ratio=new_isotopic_ratio_list)
        self.simu_x = self.o_reso.total_signal['energy_eV']
        self.simu_y = self.o_reso.total_signal['attenuation']

    def x_angstrom(self):
        """
        Convert x to angstrom
        :return: x in angstrom
        """
        _x = reso_util.ev_to_angstroms(self.o_reso.total_signal['energy_eV'])
        return _x

    def y_transmission(self):
        """
        Convert y to transmission
        :return: x in transmission
        """
        _y = self.o_reso.total_signal['transmission']
        return _y

    def xy_simu(self, angstrom=False, transmission=False):
        """
        Get x and y arrays
        :param angstrom: bool to determine the output x
        :param transmission: bool to determine the output y
        :return: x and y arrays
        """
        _x = self.o_reso.total_signal['energy_eV']
        if angstrom is True:
            _x = reso_util.ev_to_angstroms(_x)
        if transmission is True:
            _y = self.o_reso.total_signal['transmission']
        else:
            _y = self.o_reso.total_signal['attenuation']
        return _x, _y

    def __export(self, filename=None, x_axis='energy', y_axis='attenuation',
                 all_layers=False, all_elements=False, all_isotopes=False, items_to_export=None,
                 offset_us=0., source_to_detector_m=16.,
                 t_start_us=1, time_resolution_us=0.16, time_unit='us'):

        if items_to_export is not None:
            # Shape items
            items = fit_util.Items(o_reso=self.o_reso)
            items_to_export = items.shaped(items_list=items_to_export)

        self.o_reso.export(filename=filename,
                           x_axis=x_axis,
                           y_axis=y_axis,
                           all_layers=all_layers,
                           all_elements=all_elements,
                           all_isotopes=all_isotopes,
                           items_to_export=items_to_export,
                           offset_us=offset_us,
                           source_to_detector_m=source_to_detector_m,
                           t_start_us=t_start_us,
                           time_resolution_us=time_resolution_us,
                           time_unit=time_unit)
