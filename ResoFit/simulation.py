import re
import numpy as np
import pandas as pd
import ImagingReso._utilities as reso_util
from ImagingReso.resonance import Resonance
import ResoFit._utilities as fit_util
import pprint
import peakutils as pku
from ResoFit._pulse_shape import NeutronPulse


class Simulation(object):
    # Input sample name or names as str, case sensitive

    def __init__(self, energy_min=1e-5, energy_max=1000, energy_step=0.01, database='ENDF_VIII'):
        """
        initialize the a Simulation() using the Resonance() in ImagingReso

        :param energy_min:
        :type energy_min:
        :param energy_max:
        :type energy_max:
        :param energy_step:
        :type energy_step:
        :param database:
        :type database:
        """
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.database = database

        self.o_reso = Resonance(energy_min=energy_min, energy_max=energy_max, energy_step=energy_step,
                                database=database)

        self.x_simu = None
        self.y_simu = None
        self.layer_list = []

        self.x_tof_us = None
        self.y_att = None

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
        self.x_simu = self.o_reso.total_signal['energy_eV']
        self.y_simu = self.o_reso.total_signal['attenuation']

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
        self.x_simu = self.o_reso.total_signal['energy_eV']
        self.y_simu = self.o_reso.total_signal['attenuation']

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

    def xy_simu(self, x_type='energy', y_type='attenuation'):
        """
        Get x and y arrays

        :param y_type:
        :type y_type:
        :param x_type:
        :type x_type:

        :return: x and y arrays
        :rtype: array
        """
        _x = self.o_reso.total_signal['energy_eV']
        if x_type == 'lambda':
            _x = reso_util.ev_to_angstroms(_x)
        _y = self.o_reso.total_signal[y_type]

        return _x, _y

    def _convolve_beam_shape(self, convolve_proton):
        path1 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_1.dat'
        path2 = '/Users/y9z/Dropbox (ORNL)/Postdoc_Research/neutron_beam_shape/SNS/neutron_pulse/source_section_2.dat'

        neutron_pulse = NeutronPulse(path1, model_index=1)
        neutron_pulse.load_shape_each(path2)
        neutron_pulse.fit_shape(e_min=1, e_max=500,
                                drop=False, norm=True,
                                check_each=False,
                                save_fig=False,
                                overwrite_csv=False)
        neutron_pulse.fit_params(check_each=False, loglog_fit=True, overwrite_csv=False)

        e_list = self.x_simu
        # t_new = np.linspace(0.1, 30, 300)
        neutron_pulse.make_shape(e_ev=e_list, t_interp=None, for_sum=True, norm=False,
                                 convolve_proton=convolve_proton, overwrite_csv=False)
        self.neutron_pulse = neutron_pulse

        tof_beam_shape_df = neutron_pulse.shape_tof_df_interp.set_index('tof_us')
        tof_trans_df = tof_beam_shape_df * self.o_reso.total_signal['transmission']

        tof_beam_shape_df['sum'] = tof_beam_shape_df.sum(axis=1)
        tof_trans_df['sum'] = tof_trans_df.sum(axis=1)
        print(tof_beam_shape_df)
        print(tof_trans_df)

        self.x_tof_us = np.array(tof_beam_shape_df.index)
        self.y_att = 1 - np.array(tof_trans_df['sum'] / tof_beam_shape_df['sum'])

        # tof_beam_shape_df_proton = neutron_pulse.shape_tof_df_interp_proton.set_index('tof_us')
        # tof_trans_df_proton = tof_beam_shape_df_proton * self.o_reso.total_signal['transmission']
        #
        # tof_beam_shape_df_proton['sum'] = tof_beam_shape_df.sum(axis=1)
        # tof_trans_df_proton['sum'] = tof_trans_df.sum(axis=1)
        # print(tof_beam_shape_df)
        # print(tof_trans_df)
        #
        # self.y_att_proton = 1 - np.array(tof_trans_df_proton['sum'] / tof_beam_shape_df_proton['sum'])

    def peak_map(self, thres, min_dist, impr_reso=True, isotope=False):
        """
        Get peak map (eV and sigma) for each element and/or nuclide

        :param thres:
        :type thres: float
        :param min_dist:
        :type min_dist:
        :param impr_reso:
        :type impr_reso:
        :param isotope:
        :type isotope:
        :return:
        :rtype:
        """
        if len(self.layer_list) == 0:
            raise ValueError("No layer has been added.")
        _stack_sigma = self.o_reso.stack_sigma
        _layer_list = self.layer_list
        _x_energy = _stack_sigma[_layer_list[0]][_layer_list[0]]['energy_eV']
        # peak_dict = {'energy_eV': _x_energy}
        peak_dict = {}
        for _ele in _layer_list:
            _ele_sigma = _stack_sigma[_ele][_ele]['sigma_b']
            # peak_dict[_ele] = {'sigma_b': _ele_sigma, }
            peak_dict[_ele] = {}
            _peak_df = fit_util.find_peak(x=_x_energy, y=_ele_sigma,
                                          thres=thres, min_dist=min_dist, impr_reso=impr_reso)
            peak_dict[_ele]['peak'] = _peak_df
            if isotope is True:
                for _iso in self.o_reso.stack[_ele][_ele]['isotopes']['list']:
                    _iso_sigma = _stack_sigma[_ele][_ele][_iso]['sigma_b']
                    peak_dict[_ele][_iso] = {'sigma_b': _iso_sigma, }
                    _peak_df = fit_util.find_peak(x=_x_energy, y=_iso_sigma,
                                                  thres=0.5, min_dist=50, impr_reso=impr_reso)
                    peak_dict[_ele][_iso]['peak'] = _peak_df
        # print(peak_dict)
        # pprint.pprint(peak_dict)
        return peak_dict

    def plot_simu(self, y_type='attenuation', x_type='energy', mixed=True, all_layers=False, all_elements=False,
                  all_isotopes=False, items_to_plot=None, time_unit='us', offset_us=0., time_resolution_us=0.16,
                  source_to_detector_m=16., lambda_max_angstroms=1, t_start_us=1):
        if len(self.layer_list) == 0:
            raise ValueError("No layer has been added.")
        if items_to_plot is not None:
            # shape format of items
            items = fit_util.Items(o_reso=self.o_reso, database=self.database)
            items_to_plot = items.shaped(items_list=items_to_plot)

        self.o_reso.plot(y_axis=y_type, x_axis=x_type, mixed=mixed,
                         all_layers=all_layers, all_elements=all_elements,
                         all_isotopes=all_isotopes, items_to_plot=items_to_plot,
                         lambda_max_angstroms=lambda_max_angstroms,
                         source_to_detector_m=source_to_detector_m,
                         offset_us=offset_us,
                         time_resolution_us=time_resolution_us,
                         time_unit=time_unit,
                         t_start_us=t_start_us)

    def _export_simu(self, filename=None, x_axis='energy', y_axis='attenuation',
                     all_layers=False, all_elements=False, all_isotopes=False, items_to_export=None,
                     offset_us=0., source_to_detector_m=16.,
                     t_start_us=1, time_resolution_us=0.16, time_unit='us'):
        if items_to_export is not None:
            # Shape items
            items = fit_util.Items(o_reso=self.o_reso, database=self.database)
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
