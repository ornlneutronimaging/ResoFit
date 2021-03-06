import os
import re

import numpy as np
from ImagingReso.resonance import Resonance

import ResoFit._utilities as fit_util
from ResoFit._pulse_shape import NeutronPulse
from ResoFit._utilities import Layer


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

        self.o_reso = Resonance(energy_min=energy_min,
                                energy_max=energy_max,
                                energy_step=energy_step,
                                database=database)
        self.neutron_pulse = None

        self.layer_list = []
        self.layer = fit_util.Layer()

        self.x_tof_us = None
        self.y_att = None

    def add_layer(self, layer: str, thickness_mm: float, density_gcm3=np.NaN):
        """

        :param layer:
        :type layer:
        :param thickness_mm:
        :type thickness_mm:
        :param density_gcm3:
        :type density_gcm3:
        :return:
        :rtype:
        """
        self.o_reso.add_layer(formula=layer,
                              thickness=thickness_mm,
                              density=density_gcm3)
        self.layer_list.append(layer)
        self.layer.add_layer(layer=layer, thickness_mm=thickness_mm, density_gcm3=density_gcm3)

    def add_Layer(self, layer: Layer):
        """
        Add layer using Layer class

        :param layer:
        """
        for _each_layer in list(layer.info.keys()):
            self.add_layer(layer=_each_layer,
                           thickness_mm=layer.info[_each_layer]['thickness']['value'],
                           density_gcm3=layer.info[_each_layer]['density']['value'])

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
        # self.x_simu = np.array(self.o_reso.total_signal['energy_eV']).round(5)
        # self.y_simu = np.array(self.o_reso.total_signal['attenuation'])

    def get_x(self, x_type, offset_us=None, source_to_detector_m=None, t_unit='us',
              t_start_us=None, time_resolution_us=None, num_offset=0):
        """
        Get x by specified type

        :param x_type:
        :type x_type:
        :param offset_us:
        :type offset_us:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :return: x in specified type
        :rtype: np.array
        """
        fit_util.check_if_in_list(x_type, fit_util.x_type_list)
        _x = np.array(self.o_reso.total_signal['energy_eV']).round(5)
        x = fit_util.convert_energy_to(x=_x,
                                       x_type=x_type,
                                       offset_us=offset_us,
                                       source_to_detector_m=source_to_detector_m,
                                       t_unit=t_unit,
                                       t_start_us=t_start_us,
                                       time_resolution_us=time_resolution_us,
                                       num_offset=num_offset)
        return x

    def get_y(self, y_type):
        """
        Get x by specified type

        :param y_type:
        :type y_type:
        :return: y in specified type
        :rtype: np.array
        """
        fit_util.check_if_in_list(y_type, fit_util.y_type_list)
        y = self.o_reso.total_signal[y_type]
        return y

    def _convolve_beam_shapes(self, source_to_detector_m, conv_proton, proton_params={}, model_index=1):
        _file_path = os.path.abspath(os.path.dirname(__file__))
        _rel_path_to_neutron1 = 'data/_data_for_tutorial/neutron_pulse/source_section_1.dat'
        _rel_path_to_neutron2 = 'data/_data_for_tutorial/neutron_pulse/source_section_2.dat'
        path1 = os.path.join(_file_path, _rel_path_to_neutron1)
        path2 = os.path.join(_file_path, _rel_path_to_neutron2)

        self.neutron_pulse = NeutronPulse(path1, model_index=model_index)
        self.neutron_pulse.load_shape_each(path2)
        self.neutron_pulse.fit_shape(e_min=1, e_max=500,
                                     drop=False, norm=True,
                                     check_each=False,
                                     save_fig=False,
                                     overwrite_csv=False)
        self.neutron_pulse.fit_params(check_each=False, loglog_fit=True, overwrite_csv=False)

        self.neutron_pulse.make_shape(e_ev=self.get_x(x_type='energy'), t_interp=None, for_sum=True, norm=False,
                                      source_to_detector_m=source_to_detector_m,
                                      conv_proton=conv_proton, proton_params=proton_params,
                                      overwrite_csv=False)

        tof_beam_shape_df = self.neutron_pulse.shape_tof_df_interp.set_index('tof_us')
        tof_trans_df = tof_beam_shape_df * self.o_reso.total_signal['transmission']

        tof_beam_shape_df['sum'] = tof_beam_shape_df.sum(axis=1)
        tof_trans_df['sum'] = tof_trans_df.sum(axis=1)
        # print(tof_beam_shape_df)
        # print(tof_trans_df)

        self.x_tof_us = np.array(tof_beam_shape_df.index)
        self.y_att = 1 - np.array(tof_trans_df['sum'] / tof_beam_shape_df['sum'])

    def peak_map(self, x_type, y_type, thres=0.15, min_dist=20, impr_reso=True,
                 offset_us=None, source_to_detector_m=None, t_unit='us',
                 t_start_us=None, time_resolution_us=None, num_offset=0,
                 ):
        """
        Get peak map for each element and/or isotope

        :param thres:
        :type thres:
        :param min_dist:
        :type min_dist:
        :param impr_reso:
        :type impr_reso:
        :return:
        :rtype:
        """
        if len(self.layer_list) == 0:
            raise ValueError("No layer has been added.")

        _stack_signal = self.o_reso.stack_signal
        _layer_list = self.layer_list
        # _x_energy = _stack_signal[_layer_list[0]][_layer_list[0]]['energy_eV']
        _x = self.get_x(x_type=x_type,
                        offset_us=offset_us,
                        source_to_detector_m=source_to_detector_m,
                        t_unit=t_unit,
                        t_start_us=t_start_us,
                        time_resolution_us=time_resolution_us,
                        num_offset=num_offset
                        )
        # _x = sorted(_x)
        peak_map = {}
        for _ele in _layer_list:
            # Isotope
            for _iso in self.o_reso.stack[_ele][_ele]['isotopes']['list']:
                peak_map[_iso] = {}
                _iso_y = _stack_signal[_ele][_ele][_iso]['attenuation']
                _peak_df = fit_util.find_peak(x=_x, y=_iso_y, x_name='x',
                                              thres=thres, min_dist=min_dist,
                                              imprv_reso=impr_reso)
                if y_type == 'transmission':
                    _peak_df['y'] = 1 - _peak_df['y']
                peak_map[_iso]['ideal'] = _peak_df
            # Element
            peak_map[_ele] = {}
            _ele_y = _stack_signal[_ele][_ele]['attenuation']
            _peak_df = fit_util.find_peak(x=_x, y=_ele_y, x_name='x',
                                          thres=thres, min_dist=min_dist,
                                          imprv_reso=impr_reso)
            if y_type == 'transmission':
                _peak_df['y'] = 1 - _peak_df['y']
            peak_map[_ele]['ideal'] = _peak_df
        peak_map_dict = {'peak_map': peak_map,
                         'x_type': x_type,
                         'y_type': y_type}
        return peak_map_dict

    def plot(self, y_type='attenuation', x_type='energy',
             logx=False, logy=False,
             mixed=True, all_layers=False, all_elements=False,
             all_isotopes=False, items_to_plot=None, time_unit='us', offset_us=0.,
             source_to_detector_m=16.,
             t_start_us=1,
             time_resolution_us=0.16,
             ax_mpl=None,
             fmt='-', ms=2, lw=1.5, alpha=1):
        if len(self.layer_list) == 0:
            raise ValueError("No layer has been added.")
        if items_to_plot is not None:
            # shape format of items
            items = fit_util.Items(o_reso=self.o_reso, database=self.database)
            items_to_plot = items.shaped(items_list=items_to_plot)

        ax = self.o_reso.plot(y_axis=y_type, x_axis=x_type, mixed=mixed,
                              all_layers=all_layers, all_elements=all_elements,
                              all_isotopes=all_isotopes, items_to_plot=items_to_plot,
                              source_to_detector_m=source_to_detector_m,
                              offset_us=offset_us,
                              time_resolution_us=time_resolution_us,
                              time_unit=time_unit,
                              t_start_us=t_start_us,
                              ax_mpl=ax_mpl,
                              logx=logx,
                              logy=logy,
                              # plotly=plotly
                              fmt=fmt,
                              ms=ms,
                              lw=lw,
                              alpha=alpha)
        return ax

    def export(self, output_type='clip', filename=None, x_type='energy', y_type='attenuation',
               all_layers=False, all_elements=False, all_isotopes=False, items_to_export=None,
               offset_us=0., source_to_detector_m=16.,
               t_start_us=1, time_resolution_us=0.16, time_unit='us'):
        if items_to_export is not None:
            # Shape items
            items = fit_util.Items(o_reso=self.o_reso, database=self.database)
            items_to_export = items.shaped(items_list=items_to_export)

        _df = self.o_reso.export(output_type=output_type,
                                 filename=filename,
                                 x_axis=x_type,
                                 y_axis=y_type,
                                 all_layers=all_layers,
                                 all_elements=all_elements,
                                 all_isotopes=all_isotopes,
                                 items_to_export=items_to_export,
                                 offset_us=offset_us,
                                 source_to_detector_m=source_to_detector_m,
                                 t_start_us=t_start_us,
                                 time_resolution_us=time_resolution_us,
                                 time_unit=time_unit)
        return _df
