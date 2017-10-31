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
import ResoFit._utilities as fit_util
import pandas as pd
import pprint


class FitResonance(Experiment):
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

        self.fit_result = None
        self.fitted_density_gcm3 = None
        self.fitted_thickness_mm = None
        self.fitted_residual = None
        self.fitted_gap = None
        self.fitted_fjac = None
        self.fitted_layer = None
        self.fitted_simulation = None
        self.layer_list = None
        self.raw_layer = None
        self.fitted_iso_result = None
        self.fitted_iso_residual = None
        self.params_for_fit = None
        self.params_for_iso_fit = None
        self.isotope_stack = {}
        self.sample_vary = None
        self.df = None

    def fit(self, raw_layer, vary='density', each_step=False):
        if vary not in ['density', 'thickness', 'none']:
            raise ValueError("'vary=' can only be one of ['density', 'thickness', 'none']")
        # Default vary is: 'density'
        self.sample_vary = vary
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
        print("-------Fitting ({})-------\nParams before:".format(vary))
        self.params_for_fit.pretty_print()
        # Fitting
        self.fit_result = minimize(y_gap_for_fitting, self.params_for_fit, method='leastsq',
                                   args=(self.exp_x_interp, self.exp_y_interp, self.layer_list,
                                         self.energy_min, self.energy_max, self.energy_step, each_step))
        # Print after
        print("Params after:")
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
        self.params_for_iso_fit = Parameters()
        self.isotope_stack[layer] = {'list': self.fitted_simulation.o_reso.stack[layer][layer]['isotopes']['list'],
                                     'ratios': self.fitted_simulation.o_reso.stack[layer][layer]['isotopes'][
                                         'isotopic_ratio']}
        _formatted_isotope_list = []
        _params_name_list = []
        # Form list of param name
        for _isotope_index in range(len(self.isotope_stack[layer]['list'])):
            _split = self.isotope_stack[layer]['list'][_isotope_index].split('-')
            _flip = _split[::-1]
            _formatted_isotope_name = ''.join(_flip)
            # _formatted_isotope_name = self.isotope_stack[layer]['list'][_isotope_index].replace('-', '_')
            _formatted_isotope_list.append(_formatted_isotope_name)
            _params_name_list = _formatted_isotope_list
        # Form Parameters() for fitting
        for _name_index in range(len(_params_name_list)):
            self.params_for_iso_fit.add(_params_name_list[_name_index],
                                        value=self.isotope_stack[layer]['ratios'][_name_index],
                                        min=0,
                                        max=1)
        # Constrain sum of isotope ratios to be 1

        # _params_name_list_temp = _params_name_list[:]
        # _constraint = '+'.join(_params_name_list_temp)
        # self.params_for_iso_fit.add('sum', expr=_constraint)

        _constraint_param = _params_name_list[-1]
        _params_name_list_temp = _params_name_list[:]
        _params_name_list_temp.remove(_constraint_param)

        _constraint = '-'.join(_params_name_list_temp)
        _constraint = '1-' + _constraint
        self.params_for_iso_fit[_constraint_param].set(expr=_constraint)

        # Print params before
        print("-------Fitting (isotopic at.%)-------\nParams before:")
        self.params_for_iso_fit.pretty_print()
        # Fitting
        self.fitted_iso_result = minimize(y_gap_for_iso_fitting, self.params_for_iso_fit, method='leastsq',
                                          args=(self.exp_x_interp, self.exp_y_interp, layer, _formatted_isotope_list,
                                                self.fitted_simulation, each_step))
        # Print params after
        print("Params after:")
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
        print('\n')

        return self.fitted_layer.info

    def plot(self, error=True, table=True, grid=True, before=False, interp=False,
             all_elements=False, all_isotopes=False, items_to_plot=None, save_fig=False):
        """

        :param error:
        :type error:
        :param table:
        :type table:
        :param grid:
        :type grid:
        :param before:
        :type before:
        :param interp:
        :type interp:
        :param all_elements:
        :type all_elements:
        :param all_isotopes:
        :type all_isotopes:
        :param items_to_plot:
        :type items_to_plot:
        :param save_fig:
        :type save_fig:
        :return:
        :rtype:
        """
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
        if self.sample_vary is None:
            raise ValueError("Vary type ['density'|'thickness'] is not set.")
        fig_title = 'Fitting result of sample (' + sample_name + ')'

        # Create pd.DataFrame
        self.df = pd.DataFrame()

        # Clear any left plt
        plt.close()
        if table is True:
            # plot table + graph
            ax1 = plt.subplot2grid(shape=(10, 10), loc=(0, 1), rowspan=8, colspan=8)
        else:
            # plot graph only
            ax1 = plt.subplot(111)

        # Plot after fitting
        ax1.plot(simu_x, simu_y, 'b-', label=simu_label, linewidth=1)
        # Save to df
        _live_df_x_label = simu_label + '_eV'
        _live_df_y_label = simu_label + '_attenuation'
        self.df[_live_df_x_label] = simu_x
        self.df[_live_df_y_label] = simu_y
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
            # Save to df
            _live_df_x_label = simu_before_label + '_eV'
            _live_df_y_label = simu_before_label + '_attenuation'
            self.df[_live_df_x_label] = simu_x
            self.df[_live_df_y_label] = simu_y_before
        # 2.
        if interp is True:
            # Plot exp. data (interpolated)
            x_interp, y_interp = self.xy_scaled(energy_max=self.energy_max, energy_min=self.energy_min,
                                                energy_step=self.energy_step,
                                                angstrom=False, transmission=False, baseline=self.baseline,
                                                offset_us=self.calibrated_offset_us,
                                                source_to_detector_m=self.source_to_detector_m)
            ax1.plot(x_interp, y_interp, 'r-.', label=exp_interp_label, linewidth=1)
            # Save to df
            _live_df_x_label = exp_interp_label + '_eV'
            _live_df_y_label = exp_interp_label + '_attenuation'
            self.df[_live_df_x_label] = x_interp
            self.df[_live_df_y_label] = y_interp
        else:
            # Plot exp. data (raw)
            exp_x = self.x_raw(angstrom=False, offset_us=self.calibrated_offset_us,
                               source_to_detector_m=self.source_to_detector_m)
            exp_y = self.y_raw(transmission=False, baseline=self.baseline)
            ax1.plot(exp_x, exp_y, 'rx', label=exp_label, markersize=2)
            # Save to df
            _df = pd.DataFrame()
            _live_df_x_label = exp_label + '_eV'
            _live_df_y_label = exp_label + '_attenuation'
            _df[_live_df_x_label] = exp_x
            _df[_live_df_y_label] = exp_y
            # Concatenate since the length of raw and simu are not the same
            self.df = pd.concat([self.df, _df], axis=1)

        # 3.
        if error is True:
            # Plot fitting differences
            error_label = 'Diff.'
            _move_below_by = 0.2
            moved_fitted_residual = self.fitted_residual - _move_below_by
            ax1.plot(simu_x, moved_fitted_residual, 'g-', label=error_label, linewidth=1, alpha=1)
            # Save to df
            _live_df_x_label = error_label + '_eV'
            _live_df_y_label = error_label + '_attenuation'
            self.df[_live_df_x_label] = simu_x
            self.df[_live_df_y_label] = moved_fitted_residual
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
                    # Save to df
                    _live_df_x_label = _element + '_eV'
                    _live_df_y_label = _element + '_attenuation'
                    self.df[_live_df_x_label] = simu_x
                    self.df[_live_df_y_label] = _y_axis
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
                        # Save to df
                        _live_df_x_label = _isotope + '_eV'
                        _live_df_y_label = _isotope + '_attenuation'
                        self.df[_live_df_x_label] = simu_x
                        self.df[_live_df_y_label] = _y_axis
        # 5.
        if items_to_plot is not None:
            # plot specified from 'items_to_plot'
            y_axis_tag = 'attenuation'
            items = fit_util.Items(o_reso=self.fitted_simulation.o_reso)
            shaped_items = items.shaped(items_list=items_to_plot)
            _signal_dict = items.values(y_axis_type=y_axis_tag)
            for _each_label in list(_signal_dict.keys()):
                ax1.plot(simu_x, _signal_dict[_each_label], '--', label=_each_label, linewidth=1, alpha=1)
                # Save to df
                _live_df_x_label = _each_label + '_eV'
                _live_df_y_label = _each_label + '_attenuation'
                self.df[_live_df_x_label] = simu_x
                self.df[_live_df_y_label] = _signal_dict[_each_label]

        # Set plot limit and captions
        fit_util.set_plt(ax1, x_max=self.energy_max, fig_title=fig_title, grid=grid)

        # Plot table
        if table is True:
            if self.fitted_iso_result is None:
                columns = list(self.fit_result.__dict__['params'].valuesdict().keys())
            else:
                columns = self.fit_result.__dict__['var_names']

            columns_to_show_dict = {}
            for _each in columns:
                _split = _each.split('_')
                if _split[0] == 'thickness':
                    _name_to_show = r'$d_{\rm{' + _split[-1] + '}}$' + ' (mm)'
                else:
                    _name_to_show = r'$\rho_{\rm{' + _split[-1] + '}}$' + ' (g/cm$^3$)'
                columns_to_show_dict[_each] = _name_to_show
            columns_to_show = list(columns_to_show_dict.values())
            rows = ['Before', 'After']
            _row_before = []
            _row_after = []
            for _each in columns:
                _row_after.append(round(self.fit_result.__dict__['params'].valuesdict()[_each], 3))
                _row_before.append(round(self.params_for_fit.valuesdict()[_each], 3))

            if self.fitted_iso_result is not None:
                _iso_columns = list(self.fitted_iso_result.__dict__['params'].valuesdict().keys())
                columns = columns + _iso_columns
                _iso_columns_to_show_dict = {}
                for _each_iso in _iso_columns:
                    _num_str = re.findall('\d+', _each_iso)[0]
                    _name_str = _each_iso[0]
                    _sup_name = r"$^{" + _num_str + "}$" + _name_str
                    _iso_columns_to_show_dict[_each_iso] = _sup_name
                _iso_columns_to_show = list(_iso_columns_to_show_dict.values())
                columns_to_show = columns_to_show + _iso_columns_to_show
                for _each in _iso_columns:
                    _row_after.append(round(self.fitted_iso_result.__dict__['params'].valuesdict()[_each], 3))
                    _row_before.append(round(self.params_for_iso_fit.valuesdict()[_each], 3))
            table = ax1.table(rowLabels=rows, colLabels=columns_to_show, cellText=[_row_before, _row_after],
                              loc='upper right',
                              bbox=[0, -0.33, 1.0, 0.18])
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            plt.tight_layout()

        if save_fig:
            _sample_name = '_'.join(self.layer_list)
            _filename = 'fitting_' + _sample_name + '.tiff'
            plt.savefig(_filename, dpi=600, transparent=True)
        else:
            plt.show()

    def export(self, filename=None):
        if self.df is None:
            raise ValueError("pd.DataFrame is empty, please run required step.")
        elif filename is None:
            self.df.to_clipboard(excel=True)
        else:
            self.df.to_csv(filename)
