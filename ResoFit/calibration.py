import ImagingReso._utilities as reso_util
import matplotlib.pyplot as plt
from lmfit import Parameters
from lmfit import minimize

import ResoFit._utilities as fit_util
from ResoFit._gap_functions import y_gap_for_calibration
# from ResoFit._gap_functions import y_gap_for_adv_calibration
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation
from math import isclose
import pandas as pd
from ResoFit._utilities import Peak


class Calibration(Simulation):
    def __init__(self, spectra_file, data_file, layer,
                 energy_min=1e-5, energy_max=1000, energy_step=0.01,
                 repeat=1, folder='data', baseline=False,
                 database='ENDF_VII'):
        """
        Initialization with passed file location and sample info

        :param spectra_file:
        :type spectra_file:
        :param data_file:
        :type data_file:
        :param layer: Layer()
        :type layer:
        :param energy_min:
        :type energy_min:
        :param energy_max:
        :type energy_max:
        :param energy_step:
        :type energy_step:
        :param repeat:
        :type repeat:
        :param folder:
        :type folder:
        :param baseline: True -> to remove baseline/background by detrend
        :type baseline: boolean
        """
        super().__init__(energy_min=energy_min,
                         energy_max=energy_max,
                         energy_step=energy_step,
                         database=database)
        for _each_layer in list(layer.info.keys()):
            self.add_layer(layer=_each_layer,
                           layer_thickness_mm=layer.info[_each_layer]['thickness']['value'],
                           layer_density_gcm3=layer.info[_each_layer]['density']['value'])
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.experiment = Experiment(spectra_file=spectra_file, data_file=data_file, repeat=repeat, folder=folder)
        self.repeat = repeat
        self.data_file = data_file
        self.spectra_file = spectra_file
        self.init_source_to_detector_m = None
        self.init_offset_us = None
        self.calibrated_offset_us = None
        self.calibrated_source_to_detector_m = None
        self.calibrate_result = None
        self.exp_x_raw_calibrated = None
        self.exp_y_raw_calibrated = None
        self.exp_x_interp_calibrated = None
        self.exp_y_interp_calibrated = None
        self.baseline = baseline
        # self.calibrated_residual = None
        self.params_to_calibrate = None
        self.raw_layer = layer
        self.database = database

        # self.peak_df_scaled = None
        # self.peak_map_full = None
        # self.peak_map_indexed = None

    def norm_to(self, file):
        if file is not None:
            self.experiment.norm_to(file=file)

    def slice(self, slice_start=None, slice_end=None):
        self.experiment.slice(slice_start=slice_start, slice_end=slice_end)

    def calibrate(self, source_to_detector_m, offset_us, vary='all', each_step=False):
        """
        calibrate the instrumental parameters: source-to-detector-distance & detector delay
        :param each_step: boolean. True -> show values and chi^2 of each step
        :param source_to_detector_m: estimated distance in m
        :param offset_us: estimated time offset in us
        :param vary: vary one of or both of 'source_to_detector' and 'offset' to calibrate (default: 'all')

        :return: lmfit MinimizerResult
        """
        self.init_source_to_detector_m = source_to_detector_m
        self.init_offset_us = offset_us
        if vary not in ['source_to_detector', 'offset', 'all', 'none']:
            raise ValueError("'vary=' can only be one of ['source_to_detector', 'offset', 'all' 'none']")
        simu_x = self.x_simu
        simu_y = self.y_simu

        source_to_detector_vary_tag = True
        offset_vary_tag = True
        if vary == 'source_to_detector':
            offset_vary_tag = False
        if vary == 'offset':
            source_to_detector_vary_tag = False
        if vary == 'none':
            source_to_detector_vary_tag = False
            offset_vary_tag = False
        self.params_to_calibrate = Parameters()
        self.params_to_calibrate.add('source_to_detector_m', value=source_to_detector_m,
                                     vary=source_to_detector_vary_tag)
        self.params_to_calibrate.add('offset_us', value=offset_us, vary=offset_vary_tag)
        # Print before
        print("+----------------- Calibration -----------------+\nParams before:")
        self.params_to_calibrate.pretty_print()
        # Use lmfit to obtain 'source_to_detector_m' & 'offset_us' to minimize 'y_gap_for_calibration'
        self.calibrate_result = minimize(y_gap_for_calibration,
                                         self.params_to_calibrate,
                                         method='leastsq',
                                         args=(simu_x, simu_y,
                                               self.energy_min, self.energy_max, self.energy_step,
                                               self.experiment, self.baseline, each_step))
        # Print after
        print("\nParams after:")
        self.calibrate_result.__dict__['params'].pretty_print()
        # Print chi^2
        # self.calibrated_residual = self.calibrate_result.__dict__['residual']
        print("Calibration chi^2 : {}\n".format(self.calibrate_result.__dict__['chisqr']))
        self.calibrated_offset_us = self.calibrate_result.__dict__['params'].valuesdict()['offset_us']
        self.calibrated_source_to_detector_m = \
            self.calibrate_result.__dict__['params'].valuesdict()['source_to_detector_m']

        # Save the calibrated experimental x & y in Calibration class
        self.exp_x_raw_calibrated = self.experiment.x_raw(x_type='energy',
                                                          offset_us=self.calibrated_offset_us,
                                                          source_to_detector_m=self.calibrated_source_to_detector_m)
        self.exp_y_raw_calibrated = self.experiment.y_raw(y_type='attenuation', baseline=self.baseline)

        self.exp_x_interp_calibrated, self.exp_y_interp_calibrated = self.experiment.xy_scaled(
            energy_min=self.energy_min,
            energy_max=self.energy_max,
            energy_step=self.energy_step,
            offset_us=self.calibrated_offset_us,
            source_to_detector_m=self.calibrated_source_to_detector_m,
            baseline=self.baseline)

        return self.calibrate_result

    def __find_peak(self, thres=0.15, min_dist=2):
        # load detected peak with x in image number
        if self.calibrate_result is None:
            raise ValueError("Instrument params have not been calibrated.")
        self.experiment.find_peak(thres=thres, min_dist=min_dist)
        self.experiment.scale_peak_with_ev(energy_min=self.energy_min,
                                           energy_max=self.energy_max,
                                           calibrated_offset_us=self.calibrated_offset_us,
                                           calibrated_source_to_detector_m=self.calibrated_source_to_detector_m)
        assert self.experiment.o_peak.peak_df_scaled is not None
        return self.experiment.o_peak.peak_df_scaled

    def index_peak(self, thres, min_dist, map_thres=0.01, map_min_dist=20, rel_tol=5e-3, impr_reso=True, isotope=False):
        if self.experiment.o_peak is None:
            self.__find_peak(thres=thres, min_dist=min_dist)
        # find peak map using Simulation.peak_map()
        _peak_map = self.peak_map(thres=map_thres, min_dist=map_min_dist, impr_reso=impr_reso, isotope=isotope)
        # pass peak map to Peak()
        self.experiment.o_peak.peak_map_full = _peak_map
        # index using Peak()
        self.experiment.o_peak.index(_peak_map, rel_tol=rel_tol)

        return self.experiment.o_peak.peak_map_indexed

    def analyze_peak(self, report=False, show_fit=False):
        if self.experiment.o_peak is None:
            raise AttributeError("Please run 'Calibration.index_peak()' before peak analysis.")
        self.experiment.o_peak.analyze(report=report)
        self.experiment.o_peak.fill_peak_span(offset_us=self.calibrated_offset_us,
                                              source_to_detector_m=self.calibrated_source_to_detector_m)
        if show_fit:
            self.experiment.o_peak.plot_fit()

        return self.experiment.o_peak.peak_map_indexed

    # def calibrate_peak_pos(self, thres=0.15, min_dist=2, vary='all', each_step=False):
    #     """
    #     calibrate the instrumental parameters: source-to-detector-distance & detector delay
    #     based on peak positions obtained from the instrument parameters after Calibration.calibrate().
    #
    #     :param thres:
    #     :type thres:
    #     :param min_dist:
    #     :type min_dist:
    #     :param vary: vary one of or both of 'source_to_detector' and 'offset' to calibrate (default: 'all')
    #     :type vary:
    #     :param each_step: True -> show values and chi^2 of each step
    #     :type each_step: boolean.
    #     :return: calibration result
    #     :rtype: lmfit MinimizerResult
    #     """
    #     if self.peak_map_indexed is None:
    #         raise ValueError('Calibrate must be run before running advanced calibration.')
    #     # self.init_source_to_detector_m = source_to_detector_m
    #     # self.init_offset_us = offset_us
    #     if vary not in ['source_to_detector', 'offset', 'all', 'none']:
    #         raise ValueError("'vary=' can only be one of ['source_to_detector', 'offset', 'all' 'none']")
    #     ideal_x = []
    #     for _ele in self.peak_map_indexed.keys():
    #         ideal_x = ideal_x + list(self.peak_map_indexed[_ele]['ideal']['x'])
    #     sorted(ideal_x)
    #     print(ideal_x)
    #
    #     source_to_detector_vary_tag = True
    #     offset_vary_tag = True
    #     if vary == 'source_to_detector':
    #         offset_vary_tag = False
    #     if vary == 'offset':
    #         source_to_detector_vary_tag = False
    #     if vary == 'none':
    #         source_to_detector_vary_tag = False
    #         offset_vary_tag = False
    #     self.params_to_calibrate = Parameters()
    #     self.params_to_calibrate.add('source_to_detector_m',
    #                                  value=self.calibrated_source_to_detector_m,
    #                                  vary=source_to_detector_vary_tag)
    #     self.params_to_calibrate.add('offset_us',
    #                                  value=self.calibrated_offset_us,
    #                                  vary=offset_vary_tag)
    #     # Print before
    #     print("-------Calibration(advanced)-------\nParams before:")
    #     self.params_to_calibrate.pretty_print()
    #     # Use lmfit to obtain 'source_to_detector_m' & 'offset_us' to minimize 'y_gap_for_calibration'
    #     self.calibrate_result = minimize(y_gap_for_adv_calibration,
    #                                      self.params_to_calibrate,
    #                                      method='leastsq',
    #                                      args=(ideal_x, thres, min_dist,
    #                                            self.experiment, each_step))
    #     # Print after
    #     print("Params after:")
    #     self.calibrate_result.__dict__['params'].pretty_print()
    #     # Print chi^2
    #     self.calibrated_residual = self.calibrate_result.__dict__['residual']
    #     print("Calibration chi^2 : {}\n".format(sum(self.calibrated_residual ** 2)))
    #     self.calibrated_offset_us = self.calibrate_result.__dict__['params'].valuesdict()['offset_us']
    #     self.calibrated_source_to_detector_m = \
    #         self.calibrate_result.__dict__['params'].valuesdict()['source_to_detector_m']
    #
    #     # Save the calibrated experimental x & y in Calibration class
    #     self.exp_x_raw_calibrated = self.experiment.x_raw(angstrom=False,
    #                                                       offset_us=self.calibrated_offset_us,
    #                                                       source_to_detector_m=self.calibrated_source_to_detector_m)
    #     self.exp_y_raw_calibrated = self.experiment.y_raw(transmission=False, baseline=self.baseline)
    #
    #     self.exp_x_interp_calibrated, self.exp_y_interp_calibrated = self.experiment.xy_scaled(
    #         energy_min=self.energy_min,
    #         energy_max=self.energy_max,
    #         energy_step=self.energy_step,
    #         offset_us=self.calibrated_offset_us,
    #         source_to_detector_m=self.calibrated_source_to_detector_m,
    #         baseline=self.baseline)
    #
    #     return self.calibrate_result

    def plot(self, table=True, grid=True, before=False, interp=False, total=False,
             all_elements=False, all_isotopes=False, items_to_plot=None,
             peak_mark=True, peak_id='indexed',
             save_fig=False):
        """

        :param peak_mark:
        :type peak_mark:
        :param peak_id:
        :type peak_id:
        :param total:
        :type total:
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
        if all_elements is True:
            if len(self.layer_list) == 1:
                raise ValueError("'all_elements=True' has not effect on the plot if only one element was involved. ")
        if peak_id not in ['indexed', 'all']:
            raise ValueError("'peak=' must be one of ['indexed', 'all'].")
        simu_label = 'Ideal'
        exp_label = 'Exp'
        exp_before_label = 'Exp_init'
        exp_interp_label = 'Exp_interp'
        sample_name = ' & '.join(self.layer_list)
        fig_title = 'Calibration result of sample (' + sample_name + ')'

        # clear any left plt
        plt.close()

        # plot table + graph
        if table is True:
            ax1 = plt.subplot2grid(shape=(10, 10), loc=(0, 1), rowspan=8, colspan=8)
        # plot graph only
        else:
            ax1 = plt.subplot(111)

        # Plot simulated total signal
        if total is True:
            ax1.plot(self.x_simu, self.y_simu, 'b-', label=simu_label, linewidth=1)

        """Plot options"""

        # 1.
        if before is True:
            # Plot the raw data before fitting
            ax1.plot(self.experiment.x_raw(offset_us=self.init_offset_us,
                                           source_to_detector_m=self.init_source_to_detector_m),
                     self.experiment.y_raw(baseline=self.baseline),
                     linestyle='-', linewidth=1,
                     marker='o', markersize=2,
                     color='c', label=exp_before_label)
        # 2.
        if interp is True:
            # plot the interpolated raw data
            ax1.plot(self.exp_x_interp_calibrated, self.exp_y_interp_calibrated,
                     'y:', label=exp_interp_label, linewidth=1)
        else:
            # plot the calibrated raw data
            ax1.plot(self.exp_x_raw_calibrated, self.exp_y_raw_calibrated,
                     linestyle='-', linewidth=1,
                     marker='o', markersize=2,
                     color='r', label=exp_label)

        # 3.
        if all_elements is True:
            # show signal from each elements
            _stack_signal = self.o_reso.stack_signal
            _stack = self.o_reso.stack
            y_axis_tag = 'attenuation'
            for _layer in _stack.keys():
                for _element in _stack[_layer]['elements']:
                    _y_axis = _stack_signal[_layer][_element][y_axis_tag]
                    ax1.plot(self.x_simu, _y_axis, label="{}".format(_element), linewidth=1, alpha=0.85)
        # 4.
        if all_isotopes is True:
            # show signal from each isotopes
            _stack_signal = self.o_reso.stack_signal
            _stack = self.o_reso.stack
            y_axis_tag = 'attenuation'
            for _layer in _stack.keys():
                for _element in _stack[_layer]['elements']:
                    for _isotope in _stack[_layer][_element]['isotopes']['list']:
                        _y_axis = _stack_signal[_layer][_element][_isotope][y_axis_tag]
                        ax1.plot(self.x_simu, _y_axis, label="{}".format(_isotope), linewidth=1, alpha=1)
        # 5.
        if items_to_plot is not None:
            # plot specified from 'items_to_plot'
            y_axis_tag = 'attenuation'
            items = fit_util.Items(o_reso=self.o_reso, database=self.database)
            items.shaped(items_list=items_to_plot)
            _signal_dict = items.values(y_axis_type=y_axis_tag)
            for _each_label in list(_signal_dict.keys()):
                ax1.plot(self.x_simu, _signal_dict[_each_label], '--', label=_each_label, linewidth=1, alpha=1)

        # plot peaks detected and indexed
        if self.experiment.o_peak and self.experiment.o_peak.peak_map_indexed is not None:
            _peak_df_scaled = self.experiment.o_peak.peak_df_scaled
            _peak_map_indexed = self.experiment.o_peak.peak_map_indexed
            _peak_map_full = self.experiment.o_peak.peak_map_full
            if peak_mark is True:
                ax1.plot(_peak_df_scaled['x'],
                         _peak_df_scaled['y'],
                         'kx', label='_nolegend_')
            ax1.set_ylim(bottom=-0.1)
            for _ele_name in _peak_map_indexed.keys():
                if peak_id == 'all':
                    ax1.plot(_peak_map_full[_ele_name]['peak']['x'],
                             [-0.05] * len(_peak_map_full[_ele_name]['peak']['x']),
                             '|', ms=10,
                             label=_ele_name)
                elif peak_id == 'indexed':
                    ax1.plot(_peak_map_indexed[_ele_name]['exp']['x'],
                             [-0.05] * len(_peak_map_indexed[_ele_name]['exp']['x']),
                             '|', ms=8,
                             label=_ele_name)
                if 'peak_span' in _peak_map_indexed[_ele_name].keys():
                    _data_point_x = _peak_map_indexed[_ele_name]['peak_span']['energy_ev']
                    _data_point_y = _peak_map_indexed[_ele_name]['peak_span']['y']
                    ax1.scatter(_data_point_x,
                                _data_point_y,
                                label='_nolegend_')

        # Set plot limit and captions
        fit_util.set_plt(ax1, x_max=self.energy_max, fig_title=fig_title, grid=grid)

        # Plot table
        if table is True:
            # ax2 = plt.subplot2grid(shape=(10, 7), loc=(0, 1), rowspan=4, colspan=5)
            # ax2.axis('off')
            columns = list(self.calibrate_result.__dict__['params'].valuesdict().keys())
            columns_to_show = [r'$L$ (m)', r'$\Delta t$ ($\rm{\mu}$s)']
            rows = ['Before', 'After']
            _row_before = []
            _row_after = []
            for _each in columns:
                _row_after.append(self.calibrate_result.__dict__['params'].valuesdict()[_each])
                _row_before.append(self.params_to_calibrate.valuesdict()[_each])
            table = ax1.table(rowLabels=rows, colLabels=columns_to_show,  # colWidths=
                              cellText=[_row_before, _row_after],  # rows of data values
                              bbox=[0, -0.33, 1.0, 0.18]  # [left,bottom,width,height]
                              )
            # table.scale(0.5, 1)
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            plt.tight_layout()

        if save_fig:
            _sample_name = '_'.join(self.layer_list)
            _filename = 'calibration_' + _sample_name + '.png'
            plt.savefig(_filename, dpi=600, transparent=True)
            plt.close()
        # else:
        #     plt.show()

    def export_simu(self, filename=None, x_axis='energy', y_axis='attenuation',
                    all_layers=False, all_elements=False, all_isotopes=False, items_to_export=None,
                    t_start_us=1, time_resolution_us=0.16, time_unit='us'):
        if items_to_export is not None:
            # Shape items
            items = fit_util.Items(o_reso=self.o_reso, database=self.database)
            items_to_export = items.shaped(items_list=items_to_export)

        self._export_simu(filename=filename,
                          x_axis=x_axis,
                          y_axis=y_axis,
                          all_layers=all_layers,
                          all_elements=all_elements,
                          all_isotopes=all_isotopes,
                          items_to_export=items_to_export,
                          offset_us=self.calibrated_offset_us,
                          source_to_detector_m=self.calibrated_source_to_detector_m,
                          t_start_us=t_start_us,
                          time_resolution_us=time_resolution_us,
                          time_unit=time_unit)
