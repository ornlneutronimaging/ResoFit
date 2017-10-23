import matplotlib.pyplot as plt
import peakutils as pku
from lmfit import Parameters
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
import numpy as np
from lmfit import minimize
import re
from ResoFit._gap_functions import y_gap_for_fitting
from ResoFit._gap_functions import y_gap_for_iso_fitting
import periodictable as pt
from ResoFit._utilities import Layer
from ResoFit._utilities import shape_item_to_plot
import pandas as pd
import pprint


class FitResonance(Experiment):
    fit_result = None
    fitted_density_gcm3 = None
    fitted_thickness_mm = None
    fitted_residual = None
    fitted_gap = None
    fitted_fjac = None
    fitted_layer = None
    fitted_simulation = None
    layer_list = None
    raw_layer = None
    fitted_iso_result = None
    fitted_iso_residual = None
    params_for_fit = None
    isotope_stack = {}

    def __init__(self, spectra_file, data_file,
                 calibrated_offset_us, calibrated_source_to_detector_m,
                 folder, repeat=1, baseline=False,
                 norm_to_file=None, slice_start=None, slice_end=None,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01):
        super().__init__(spectra_file=spectra_file, data_file=data_file, folder=folder, repeat=repeat)
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.calibrated_offset_us = calibrated_offset_us
        self.calibrated_source_to_detector_m = calibrated_source_to_detector_m
        self.raw_layer = None
        self.slice(slice_start=slice_start, slice_end=slice_end)
        self.baseline = baseline
        if norm_to_file is not None:
            self.norm_to(norm_to_file)
        self.exp_x_interp, self.exp_y_interp = self.xy_scaled(energy_min=self.energy_min,
                                                              energy_max=self.energy_max,
                                                              energy_step=self.energy_step,
                                                              angstrom=False, transmission=False,
                                                              offset_us=self.calibrated_offset_us,
                                                              source_to_detector_m=self.calibrated_source_to_detector_m,
                                                              baseline=self.baseline)

    def fit(self, raw_layer, vary='density', isotope=False, isotope_vary_tag_list=[], each_step=False):
        if vary not in ['density', 'thickness', 'none']:
            raise ValueError("'vary=' can only be one of ['density', 'thickness', 'none']")
        # Default vary is: 'density'
        thickness_vary_tag = False
        density_vary_tag = True
        if vary == 'thickness':
            thickness_vary_tag = True
            density_vary_tag = False
        if vary == 'none':
            density_vary_tag = False
        self.raw_layer = raw_layer

        '''Load params'''

        self.layer_list = list(raw_layer.info.keys())
        self.params_for_fit = Parameters()
        for _each_layer in self.layer_list:
            if self.raw_layer.info[_each_layer]['density']['value'] is np.NaN:
                self.raw_layer.info[_each_layer]['density']['value'] = pt.elements.isotope(_each_layer).density
            self.params_for_fit.add('thickness_mm_' + _each_layer,
                                    value=self.raw_layer.info[_each_layer]['thickness']['value'],
                                    vary=thickness_vary_tag,
                                    min=0)
            self.params_for_fit.add('density_gcm3_' + _each_layer,
                                    value=self.raw_layer.info[_each_layer]['density']['value'],
                                    vary=density_vary_tag,
                                    min=0)
        # Print before
        print("Params before '{}' fitting:".format(vary))
        self.params_for_fit.pretty_print()
        # Fitting
        self.fit_result = minimize(y_gap_for_fitting, self.params_for_fit, method='leastsq',
                                   args=(self.exp_x_interp, self.exp_y_interp, self.layer_list,
                                         self.energy_min, self.energy_max, self.energy_step, each_step))
        # Print after
        print("Params after '{}' fitting:".format(vary))
        self.fit_result.__dict__['params'].pretty_print()
        # Print chi^2
        self.fitted_residual = self.fit_result.__dict__['residual']
        print("Fitting chi^2 : {}\n".format(sum(self.fitted_residual ** 2)))

        '''Export fitted params as Layer()'''

        # Save the fitted 'density' or 'thickness' in Layer()
        self.fitted_layer = Layer()
        for _each_layer in self.layer_list:
            self.fitted_layer.add_layer(layer=_each_layer,
                                        thickness_mm=self.fit_result.__dict__['params'].valuesdict()[
                                            'thickness_mm_' + _each_layer],
                                        density_gcm3=self.fit_result.__dict__['params'].valuesdict()[
                                            'density_gcm3_' + _each_layer])
        # self.fitted_fjac = self.fit_result.__dict__['fjac']
        # print(self.fit_result.__dict__['fjac'][0])

        '''Create fitted simulation'''

        self.fitted_simulation = Simulation(energy_min=self.energy_min,
                                            energy_max=self.energy_max,
                                            energy_step=self.energy_step)
        for each_layer in self.layer_list:
            self.fitted_simulation.add_layer(layer=each_layer,
                                             layer_thickness_mm=self.fitted_layer.info[each_layer]['thickness'][
                                                 'value'],
                                             layer_density_gcm3=self.fitted_layer.info[each_layer]['density']['value'])
        return self.fit_result

    def fit_iso(self, layer, each_step=False):
        """

        :param layer:
        :type layer:
        :param each_step:
        :type each_step:
        :return:
        :rtype:
        """
        params_for_iso_fit = Parameters()
        self.isotope_stack[layer] = {'list': self.fitted_simulation.o_reso.stack[layer][layer]['isotopes']['list'],
                                     'ratios': self.fitted_simulation.o_reso.stack[layer][layer]['isotopes'][
                                         'isotopic_ratio']}
        _formatted_isotope_list = []
        _params_name_list = []
        # _restriction_list = []
        for _isotope_index in range(len(self.isotope_stack[layer]['list'])):
            _formatted_isotope_name = self.isotope_stack[layer]['list'][_isotope_index].replace('-', '_')
            _formatted_isotope_list.append(_formatted_isotope_name)
            _params_name_list.append('isotope_ratio_' + _formatted_isotope_name)

        for _name_index in range(len(_params_name_list)):
            params_for_iso_fit.add(_params_name_list[_name_index],
                                   value=self.isotope_stack[layer]['ratios'][_name_index],
                                   min=0,
                                   max=1)

        # Restriction is not working as expected ###
        # Create expr to restrict the total to be 1.0
        # for _each_param_name in _params_name_list:
        #     _params_name_list_temp = _params_name_list[:]
        #     _params_name_list_temp.remove(_each_param_name)
        #     _params_name_list_temp.insert(0, '1.0')
        #     _restriction_list.append(' - '.join(_params_name_list_temp))
        #     print(_params_name_list_temp)
        # print(_restriction_list)
        # params_for_iso_fit.pretty_print()
        # for _i in range(len(_restriction_list)):
        #     params_for_iso_fit[_params_name_list[_i]].set(expr=_restriction_list[_i], vary=True)

        # Print params before
        print("Params before 'isotope' fitting:")
        params_for_iso_fit.pretty_print()
        # Fitting
        self.fitted_iso_result = minimize(y_gap_for_iso_fitting, params_for_iso_fit, method='leastsq',
                                          args=(self.exp_x_interp, self.exp_y_interp, layer, _formatted_isotope_list,
                                                self.fitted_simulation, each_step))
        # Print params after
        print("Params after 'isotope' fitting:")
        self.fitted_iso_result.__dict__['params'].pretty_print()
        # Print chi^2
        self.fitted_iso_residual = self.fitted_iso_result.__dict__['residual']
        print("Fit iso chi^2 : {}\n".format(sum(self.fitted_iso_residual ** 2)))

        return

    def molar_conc(self):
        molar_conc_units = 'mol/cm3'
        print("Molar-conc. ({})\tBefore_fit\tAfter_fit".format(molar_conc_units))
        for _each_layer in self.layer_list:
            molar_mass_value = self.fitted_simulation.o_reso.stack[_each_layer][_each_layer]['molar_mass']['value']
            molar_mass_units = self.fitted_simulation.o_reso.stack[_each_layer][_each_layer]['molar_mass']['units']
            # Adding molar_mass to fitted_layer info
            self.fitted_layer.info[_each_layer]['molar_mass']['value'] = molar_mass_value
            self.fitted_layer.info[_each_layer]['molar_mass']['units'] = molar_mass_units
            # Adding molar_mass to raw_layer info
            self.raw_layer.info[_each_layer]['molar_mass']['value'] = molar_mass_value
            self.raw_layer.info[_each_layer]['molar_mass']['units'] = molar_mass_units
            # Adding molar_concentration to fitted_layer info
            molar_conc_value = self.fitted_layer.info[_each_layer]['density']['value'] / molar_mass_value
            self.fitted_layer.info[_each_layer]['molar_conc']['value'] = molar_conc_value
            self.fitted_layer.info[_each_layer]['molar_conc']['units'] = molar_conc_units
            # Calculate starting molar_concentration and fitted_layer info
            start_molar_conc_value = self.raw_layer.info[_each_layer]['density']['value'] / molar_mass_value
            self.raw_layer.info[_each_layer]['molar_conc']['value'] = start_molar_conc_value
            self.raw_layer.info[_each_layer]['molar_conc']['units'] = molar_conc_units
            # molar_conc_output[_each_layer] = {'Before_fit': start_molar_conc_value,
            #                                   'After_fit': molar_conc_value}
            print("{}\t{}\t{}".format(_each_layer, start_molar_conc_value, molar_conc_value))

        return self.fitted_layer.info

    def plot(self, error=True, table=True, grid=True, before=False, interp=False,
             all_elements=False, all_isotopes=False, items_to_plot=None):
        """"""
        # Form signals from fitted_layer
        if self.fitted_simulation is None:
            self.fitted_simulation = Simulation(energy_min=self.energy_min,
                                                energy_max=self.energy_max,
                                                energy_step=self.energy_step)
            for each_layer in self.layer_list:
                self.fitted_simulation.add_layer(layer=each_layer,
                                                 layer_thickness_mm=self.fitted_layer.info[each_layer]['thickness'][
                                                     'value'],
                                                 layer_density_gcm3=self.fitted_layer.info[each_layer]['density'][
                                                     'value'])
        simu_x, simu_y = self.fitted_simulation.xy_simu(angstrom=False, transmission=False)

        # Get plot labels
        simu_label = 'Fit'
        simu_before_label = 'Fit_init'
        exp_label = 'Exp'
        exp_interp_label = 'Exp_interp'
        sample_name = ' & '.join(self.layer_list)
        fig_title = 'Fitting result of sample ' + '(' + sample_name + ')'

        if table is True:
            # plot table + graph
            ax1 = plt.subplot2grid(shape=(10, 10), loc=(0, 1), rowspan=8, colspan=8)
        else:
            # plot graph only
            ax1 = plt.subplot(111)

        # Plot after fitting
        ax1.plot(simu_x, simu_y, 'b-', label=simu_label, linewidth=1)

        """Plot options"""

        # 1.
        if before is True:
            # Plot before fitting
            # Form signals from raw_layer
            simulation = Simulation(energy_min=self.energy_min,
                                    energy_max=self.energy_max,
                                    energy_step=self.energy_step)
            for each_layer in self.layer_list:
                simulation.add_layer(layer=each_layer,
                                     layer_thickness_mm=self.raw_layer.info[each_layer]['thickness']['value'],
                                     layer_density_gcm3=self.raw_layer.info[each_layer]['density']['value'])
            simu_x, simu_y_before = simulation.xy_simu(angstrom=False, transmission=False)
            ax1.plot(simu_x, simu_y_before,
                     'c-.', label=simu_before_label, linewidth=1)
        # 2.
        if interp is True:
            # Plot exp. data (interpolated)
            x_interp, y_interp = self.xy_scaled(energy_max=self.energy_max, energy_min=self.energy_min,
                                                energy_step=self.energy_step,
                                                angstrom=False, transmission=False, baseline=self.baseline,
                                                offset_us=self.calibrated_offset_us,
                                                source_to_detector_m=self.source_to_detector_m)
            ax1.plot(x_interp, y_interp, 'r-.', label=exp_interp_label, linewidth=1)
        else:
            # Plot exp. data (raw)
            ax1.plot(self.x_raw(angstrom=False, offset_us=self.calibrated_offset_us,
                                source_to_detector_m=self.source_to_detector_m),
                     self.y_raw(transmission=False, baseline=self.baseline),
                     'rx', label=exp_label, markersize=2)
        # 3.
        if error is True:
            # Plot fitting differences
            ax1.plot(simu_x, self.fitted_residual - 0.2, 'g-', label='Diff.', linewidth=1, alpha=1)
        # 4.
        if all_elements is True:
            # show signal from each elements
            _stack_signal = self.fitted_simulation.o_reso.stack_signal
            _stack = self.fitted_simulation.o_reso.stack
            y_axis_tag = 'attenuation'
            for _layer in _stack.keys():
                for _element in _stack[_layer]['elements']:
                    _y_axis = _stack_signal[_layer][_element][y_axis_tag]
                    ax1.plot(simu_x, _y_axis, label="{}".format(_element), linewidth=1, alpha=0.85)
        # 4.
        if all_isotopes is True:
            # show signal from each isotopes
            _stack_signal = self.fitted_simulation.o_reso.stack_signal
            _stack = self.fitted_simulation.o_reso.stack
            y_axis_tag = 'attenuation'
            for _layer in _stack.keys():
                for _element in _stack[_layer]['elements']:
                    for _isotope in _stack[_layer][_element]['isotopes']['list']:
                        _y_axis = _stack_signal[_layer][_element][_isotope][y_axis_tag]
                        ax1.plot(simu_x, _y_axis, label="{}".format(_isotope), linewidth=1, alpha=1)
        # 5.
        if items_to_plot is not None:
            # plot specified from 'items_to_plot'
            _stack_signal = self.fitted_simulation.o_reso.stack_signal
            y_axis_tag = 'attenuation'

            for _path_to_plot in items_to_plot:
                if type(_path_to_plot) is not list:
                    _path_to_plot = shape_item_to_plot(_path_to_plot)
                _path_to_plot = list(_path_to_plot)
                _live_path = _stack_signal
                _label = _path_to_plot[-1]#"/".join(_path_to_plot)
                while _path_to_plot:
                    _item = _path_to_plot.pop(0)
                    _live_path = _live_path[_item]
                _y_axis = _live_path[y_axis_tag]
                ax1.plot(simu_x, _y_axis, '--', label=_label, linewidth=1, alpha=1)

        ax1.set_xlim([0, self.energy_max])
        ax1.set_ylim(ymax=1.01)
        ax1.set_title(fig_title)
        ax1.set_xlabel('Energy (eV)')
        ax1.set_ylabel('Attenuation')
        ax1.legend(loc='best')
        if grid is True:
            # ax1.set_xticks(np.arange(0, 100, 10))
            # ax1.set_yticks(np.arange(0, 1., 0.1))
            ax1.grid()

        # Plot table
        if table is True:
            columns = self.fit_result.__dict__['var_names']
            rows = ['Before', 'After']
            _row_before = []
            _row_after = []
            for _each in columns:
                _row_after.append(self.fit_result.__dict__['params'].valuesdict()[_each])
                _row_before.append(self.params_for_fit.valuesdict()[_each])
            table = ax1.table(rowLabels=rows, colLabels=columns, cellText=[_row_before, _row_after], loc='upper right',
                              bbox=[0, -0.33, 1.0, 0.18])
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            plt.tight_layout()

        plt.show()
        # plt.savefig('test.tiff')

    # def export(self, filename=None):
    #     _x_axis = self.total_signal['energy_eV']
    #     x_axis_label = None
    #     df = pd.DataFrame()
    #
    #     """X-axis"""
    #     # determine values and labels for x-axis with options from
    #     # 'energy(eV)' & 'lambda(A)' & 'time(us)' & 'image number(#)'
    #     x_axis_label = 'Energy (eV)'
    #     df[x_axis_label] = _x_axis
    #
    #     """Y-axis"""
    #     df['Total_'+y_axis_tag] = _y_axis
    #             # export based on specified level : layer|element|isotope
    #             if all_layers:
    #                 for _compound in _stack.keys():
    #                     _y_axis = _stack_signal[_compound][y_axis_tag]
    #                     df[_compound] = _y_axis
    #
    #             if all_elements:
    #                 for _compound in _stack.keys():
    #                     for _element in _stack[_compound]['elements']:
    #                         _y_axis = _stack_signal[_compound][_element][y_axis_tag]
    #                         df[_compound + '/' + _element] = _y_axis
    #
    #             if all_isotopes:
    #                 for _compound in _stack.keys():
    #                     for _element in _stack[_compound]['elements']:
    #                         for _isotope in _stack[_compound][_element]['isotopes']['list']:
    #                             _y_axis = _stack_signal[_compound][_element][_isotope][y_axis_tag]
    #                             df[_compound + '/' + _element + '/' + _isotope] = _y_axis
    #         else:
    #             # export specified transmission or attenuation
    #             for _path_to_export in items_to_export:
    #                 _path_to_export = list(_path_to_export)
    #                 _live_path = _stack_signal
    #                 _label = "/".join(_path_to_export)
    #                 while _path_to_export:
    #                     _item = _path_to_export.pop(0)
    #                     _live_path = _live_path[_item]
    #                 _y_axis = _live_path[y_axis_tag]
    #                 df[_label] = _y_axis
    #     else:
    #         # export sigma
    #         _stack_sigma = self.stack_sigma
    #         y_axis_tag = 'sigma_b'
    #         if items_to_export is None:
    #             for _compound in _stack.keys():
    #                 for _element in _stack[_compound]['elements']:
    #                     _y_axis = _stack_sigma[_compound][_element][y_axis_tag]
    #                     df[_compound + '/' + _element + '/atoms_per_cm3'] = _stack[_compound]['atoms_per_cm3'][_element]
    #                     df[_compound + '/' + _element] = _y_axis
    #                     if all_isotopes:
    #                         for _isotope in _stack[_compound][_element]['isotopes']['list']:
    #                             _y_axis = _stack_sigma[_compound][_element][_isotope][y_axis_tag]
    #                             df[_compound + '/' + _element + '/' + _isotope] = _y_axis
    #         else:
    #             # export specified sigma
    #             for _path_to_export in items_to_export:
    #                 if len(_path_to_export) == 1:
    #                     raise ValueError(
    #                         "Getting total sigma of '{}' at layer level is not supported. "
    #                         "If it is a single element layer, please follow ['layer', 'element'] format.".format(
    #                             _path_to_export[0]))
    #                 _path_to_export = list(_path_to_export)
    #                 _live_path = _stack_sigma
    #                 _label = "/".join(_path_to_export)
    #                 while _path_to_export:
    #                     _item = _path_to_export.pop(0)
    #                     _live_path = _live_path[_item]
    #                 _y_axis = _live_path[y_axis_tag]
    #                 df[_label] = _y_axis
    #
    #     if filename is None:
    #         df.to_clipboard(excel=True)
    #     else:
    #         df.to_csv(filename)
    #
