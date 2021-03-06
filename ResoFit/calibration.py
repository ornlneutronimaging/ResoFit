import ImagingReso._utilities as reso_util
import matplotlib.pyplot as plt
from lmfit import Parameters
from lmfit import minimize
import pandas as pd
from itertools import cycle
import pprint

import ResoFit._utilities as fit_util
from ResoFit._gap_functions import y_gap_for_calibration
# from ResoFit._gap_functions import y_gap_for_adv_calibration
from ResoFit.experiment import Experiment
from ResoFit.simulation import Simulation

_exp_time_offset_us = 5.2


class Calibration(object):
    def __init__(self,
                 # Initialize ResoFit.experiment
                 spectra_file: str,
                 data_file: str,
                 folder: str,
                 exp_source_to_detector_m, exp_offset_us,
                 baseline: bool,
                 baseline_deg: int,
                 # Initialize ResoFit.simulation
                 layer: fit_util.Layer,
                 energy_min, energy_max, energy_step,
                 database: str,
                 x_type: str,
                 y_type: str):
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
        :param folder:
        :type folder:
        :param baseline: True -> to remove baseline/background by detrend
        :type baseline: boolean
        """
        self.x_type = x_type
        self.y_type = y_type
        self.energy_min = energy_min
        self.energy_max = energy_max
        self.energy_step = energy_step
        self.simulation = Simulation(energy_min=energy_min,
                                     energy_max=energy_max,
                                     energy_step=energy_step,
                                     database=database)
        self.simulation.add_Layer(layer=layer)
        self.experiment = Experiment(spectra_file=spectra_file,
                                     data_file=data_file,
                                     folder=folder,
                                     source_to_detector_m=exp_source_to_detector_m,
                                     offset_us=exp_offset_us,
                                     baseline=baseline,
                                     baseline_deg=baseline_deg)
        self.experiment.t_start_us = self.experiment.t_start_us + _exp_time_offset_us
        self.init_source_to_detector_m = exp_source_to_detector_m
        self.init_offset_us = exp_offset_us
        self.calibrated_offset_us = None
        self.calibrated_source_to_detector_m = None
        self.calibrate_result = None
        self.params_to_calibrate = None

    def calibrate(self, source_to_detector_m=None, offset_us=None, vary='all',
                  each_step=False):
        """
        calibrate the instrumental parameters: source-to-detector-distance & detector delay
        :param each_step: boolean. True -> show values and chi^2 of each step
        :param source_to_detector_m: estimated distance in m
        :param offset_us: estimated time offset in us
        :param vary: vary one of or both of 'source_to_detector' and 'offset' to calibrate (default: 'all')

        :return: lmfit MinimizerResult
        """
        # Overwrite init values if input detected
        if source_to_detector_m is None:
            source_to_detector_m = self.init_source_to_detector_m
        if offset_us is None:
            offset_us = self.init_offset_us

        vary_type_list = ['source_to_detector', 'offset', 'all', 'none']
        if vary not in vary_type_list:
            raise ValueError("'vary=' can only be one of '{}'".format(vary_type_list))
        simu_x = self.simulation.get_x(x_type='energy', offset_us=offset_us, source_to_detector_m=source_to_detector_m)
        simu_y = self.simulation.get_y(y_type='attenuation')
        _run = True
        if vary == 'all':
            source_to_detector_vary_tag = True
            offset_vary_tag = True
        elif vary == 'source_to_detector':
            source_to_detector_vary_tag = True
            offset_vary_tag = False
        elif vary == 'offset':
            source_to_detector_vary_tag = False
            offset_vary_tag = True
        else:  # vary == 'none':
            source_to_detector_vary_tag = False
            offset_vary_tag = False
            _run = False

        self.params_to_calibrate = Parameters()
        self.params_to_calibrate.add('source_to_detector_m',
                                     value=source_to_detector_m,
                                     vary=source_to_detector_vary_tag)
        self.params_to_calibrate.add('offset_us',
                                     value=offset_us,
                                     vary=offset_vary_tag)
        # Print before
        print("+----------------- Calibration -----------------+\nParams before:")
        self.params_to_calibrate.pretty_print()
        # Use lmfit to obtain 'source_to_detector_m' & 'offset_us' to minimize 'y_gap_for_calibration'
        if _run:
            self.calibrate_result = minimize(y_gap_for_calibration,
                                             self.params_to_calibrate,
                                             method='leastsq',
                                             args=(simu_x, simu_y,
                                                   self.energy_min, self.energy_max, self.energy_step,
                                                   self.experiment, 'energy', 'attenuation', each_step))
            # Print after
            print("\nParams after:")
            self.calibrate_result.__dict__['params'].pretty_print()
            # Print chi^2
            # self.calibrated_residual = self.calibrate_result.__dict__['residual']
            print("Calibration chi^2 : {}\n".format(self.calibrate_result.__dict__['chisqr']))
            self.calibrated_offset_us = self.calibrate_result.__dict__['params'].valuesdict()['offset_us']
            self.calibrated_source_to_detector_m = \
                self.calibrate_result.__dict__['params'].valuesdict()['source_to_detector_m']
            return self.calibrate_result
        else:
            self.calibrated_offset_us = offset_us
            self.calibrated_source_to_detector_m = source_to_detector_m
            print("\ncalibrate() was not run as requested, input values used:\n"
                  "calibrated_offset_us = {}\ncalibrated_source_to_detector_m = {}".format(offset_us,
                                                                                           source_to_detector_m))
            # self.experiment.xy_scaled(energy_min=self.energy_min,
            #                           energy_max=self.energy_max,
            #                           energy_step=self.energy_step,
            #                           x_type='energy',
            #                           y_type='attenuation',
            #                           offset_us=offset_us,
            #                           source_to_detector_m=source_to_detector_m,
            #                           )

    def __find_peak(self, thres, min_dist):
        # load detected peak with x in image number
        # if self.calibrate_result is None:
        if self.calibrated_source_to_detector_m is None or self.calibrated_offset_us is None:
            raise ValueError("Instrument params have not been calibrated.")
        self.experiment.find_peak(x_type=self.x_type, y_type=self.y_type,
                                  thres=thres, min_dist=min_dist)
        # self.experiment.o_peak._scale_peak_df(energy_min=self.energy_min, energy_max=self.energy_max,
        #                                       )
        return self.experiment.o_peak.peak_dict

    def index_peak(self, thres_exp, min_dist_exp, thres_map, min_dist_map, rel_tol, impr_reso=True):
        if self.experiment.o_peak is None:
            self.__find_peak(thres=thres_exp, min_dist=min_dist_exp)
        # find peak map using Simulation.peak_map()
        _peak_map_dict = self.simulation.peak_map(thres=thres_map, min_dist=min_dist_map, impr_reso=impr_reso,
                                                  x_type=self.x_type, y_type=self.y_type,
                                                  offset_us=self.calibrated_offset_us,
                                                  source_to_detector_m=self.calibrated_source_to_detector_m,
                                                  t_unit=self.experiment.t_unit,
                                                  t_start_us=self.experiment.t_start_us,
                                                  time_resolution_us=self.experiment.time_resolution_us,
                                                  num_offset=self.experiment.img_start)
        # pass peak map to Peak()
        assert _peak_map_dict['x_type'] == self.experiment.o_peak.peak_dict['x_type']
        assert _peak_map_dict['y_type'] == self.experiment.o_peak.peak_dict['y_type']
        self.experiment.o_peak.peak_map_full = _peak_map_dict['peak_map']
        # index using Peak()
        self.experiment.o_peak.index_peak(_peak_map_dict, rel_tol=rel_tol)
        # return self.experiment.o_peak.peak_map_indexed

    def analyze_peak(self, fit_model, report=False, show_fit=False):
        if self.experiment.o_peak is None:
            raise AttributeError("Please run 'Calibration.index_peak()' before peak analysis.")
        self.experiment.o_peak.analyze(report=report, fit_model=fit_model)
        if show_fit:
            self.experiment.o_peak.plot_fit()

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

    def plot(self, x_type=None, y_type=None, t_unit='us',
             index_level='iso', peak_id='indexed', peak_exp='indexed',
             peak_height=True,
             before=False, interp=False, mixed=False,
             logx=True, logy=False, table=True, grid=True, save_fig=False):
        """"""
        fit_util.check_if_in_list(peak_id, fit_util.peak_type_list)
        fit_util.check_if_in_list(peak_exp, fit_util.peak_type_list)
        fit_util.check_if_in_list(index_level, fit_util.index_level_list)

        if x_type is None:
            x_type = self.x_type
        if y_type is None:
            y_type = self.y_type

        old_colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                      '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                      '#bcbd22', '#17becf']
        marker_styles = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
        color_cycle = cycle(new_colors)
        # color_cycle_2 = cycle(new_colors)
        # color_cycle_3 = cycle(new_colors)
        # color_cycle_4 = cycle(new_colors)
        style_cycle = cycle(marker_styles)

        simu_label = 'Ideal'
        exp_label = 'Exp'
        exp_before_label = 'Exp_init'
        exp_interp_label = 'Exp_interp'
        sample_name = ' & '.join(self.simulation.layer_list)
        fig_title = "Calibration result of sample ('{}')".format(sample_name)
        fig = plt.Figure()

        # plot table + graph
        if table:
            ax1 = plt.subplot2grid(shape=(10, 10), loc=(0, 1), rowspan=8, colspan=8)
        # plot graph only
        else:
            ax1 = plt.subplot(111)

        # Plot simulated total signal
        if mixed:
            _x = self.simulation.get_x(x_type=x_type,
                                       t_unit=t_unit,
                                       offset_us=self.calibrated_offset_us,
                                       source_to_detector_m=self.calibrated_source_to_detector_m,
                                       t_start_us=self.experiment.t_start_us,
                                       time_resolution_us=self.experiment.time_resolution_us,
                                       num_offset=self.experiment.slice_start
                                       )
            _y = self.simulation.get_y(y_type=y_type)
            ax1.plot(_x, _y, 'b-', label=simu_label, linewidth=1)
        """Plot options"""
        # 1.
        if before:
            # Plot the raw data before fitting
            _x_init = self.experiment.get_x(x_type=x_type,
                                            t_unit=t_unit,
                                            offset_us=self.init_offset_us,
                                            source_to_detector_m=self.init_source_to_detector_m,
                                            )
            _y_init = self.experiment.get_y(y_type=y_type)
            ax1.plot(_x_init,
                     _y_init,
                     linestyle='-', linewidth=1,
                     marker='o', markersize=2,
                     color='c', label=exp_before_label)

        # 2.
        if interp:
            _exp_x_interp_calibrated, _exp_y_interp_calibrated = self.experiment.xy_scaled(
                x_type=x_type,
                y_type=y_type,
                energy_min=self.energy_min,
                energy_max=self.energy_max,
                energy_step=self.energy_step,
                t_unit=t_unit,
                offset_us=self.calibrated_offset_us,
                source_to_detector_m=self.calibrated_source_to_detector_m,
            )
            # plot the interpolated raw data
            ax1.plot(_exp_x_interp_calibrated,
                     _exp_y_interp_calibrated,
                     'r:', label=exp_interp_label, linewidth=1)
        else:
            # plot the calibrated raw data
            _x_cali = self.experiment.get_x(x_type=x_type,
                                            t_unit=t_unit,
                                            offset_us=self.calibrated_offset_us,
                                            source_to_detector_m=self.calibrated_source_to_detector_m)
            _y_cali = self.experiment.get_y(y_type=y_type)
            ax1.plot(_x_cali,
                     _y_cali,
                     linestyle='-', linewidth=1,
                     marker='o', markersize=2,
                     color='r', label=exp_label)
        if peak_exp == 'all':
            # _peak_x_exp = fit_util.convert_exp_peak_df(x_type=x_type, peak_df=_peak_df_scaled, t_unit=t_unit)
            _peak_df_scaled = self.experiment.o_peak.peak_dict['df']
            _peak_x_exp = _peak_df_scaled['x']
            # if x_type == 'time':
            #     _peak_x_exp = fit_util.convert_s(x=_peak_x_exp, t_unit=t_unit)
            # _peak_y_exp = fit_util.convert_attenuation_to(y_type=y_type, y=_peak_df_scaled['y'])
            _peak_y_exp = _peak_df_scaled['y']
            ax1.scatter(_peak_x_exp,
                        _peak_y_exp,
                        c='k',
                        marker='x',
                        # s=30,
                        # marker='o',
                        # facecolors='none',
                        # edgecolors='k',
                        label='_nolegend_')
        # plot peaks detected and indexed
        if self.experiment.o_peak is not None:
            if self.experiment.o_peak.peak_map_indexed_dict is not None:
                if y_type == 'transmission':
                    _start_point = 1
                    ax1.set_ylim(top=1.1, bottom=-0.01)
                    _pos = 1.05
                else:
                    _start_point = 0
                    ax1.set_ylim(top=1.01, bottom=-0.1)
                    _pos = -0.05
                _peak_map_indexed = self.experiment.o_peak.peak_map_indexed_dict['peak_map_indexed']
                _peak_map_full = self.experiment.o_peak.peak_map_full

                if index_level == 'iso':
                    _peak_name_list = [_name for _name in _peak_map_indexed.keys() if '-' in _name]
                else:
                    _peak_name_list = [_name for _name in _peak_map_indexed.keys() if '-' not in _name]

                if peak_id == 'all':
                    _current_peak_map = _peak_map_full
                    # _tag = 'ideal'
                else:  # peak_id == 'indexed'
                    _current_peak_map = _peak_map_indexed
                _tag = 'ideal'

                for _peak_name in _peak_name_list:
                    if len(_current_peak_map[_peak_name][_tag]) > 0:

                        _peak_x = _current_peak_map[_peak_name][_tag]['x']
                        _peak_y = _current_peak_map[_peak_name][_tag]['y']
                        if peak_exp == 'indexed':
                            _legend_name = '_nolegend_'
                        else:
                            _legend_name = _peak_name
                        _current_color = next(color_cycle)
                        _current_style = next(style_cycle)
                        ax1.plot(_peak_x,
                                 [_pos] * len(_peak_x),
                                 '|', ms=10,
                                 color=_current_color,
                                 label=_legend_name)
                        if peak_height:
                            ax1.plot(_peak_x,
                                     _peak_y,
                                     '_',
                                     # marker=next(style_cycle_1),
                                     # ms=4,
                                     color=_current_color,
                                     label='_nolegend_')
                            ax1.vlines(_peak_x,
                                       _start_point,
                                       _peak_y,
                                       color=_current_color,
                                       alpha=1,
                                       label='_nolegend_')

                        if peak_exp == 'indexed':
                            _peak_x_exp = _peak_map_indexed[_peak_name]['exp']['x']
                            _peak_y_exp = _peak_map_indexed[_peak_name]['exp']['y']
                            ax1.scatter(_peak_x_exp,
                                        _peak_y_exp,
                                        marker=_current_style,
                                        # ms=4,
                                        color=_current_color,
                                        label=_peak_name)

                        if 'peak_span' in _peak_map_indexed[_peak_name].keys():
                            if len(_peak_map_indexed[_peak_name]['exp']) > 0:
                                _data_point_x = _peak_map_indexed[_peak_name]['peak_span']['x']
                                _data_point_y = _peak_map_indexed[_peak_name]['peak_span']['y']
                                ax1.scatter(_data_point_x,
                                            _data_point_y,
                                            label='_nolegend_')

        # Set plot limit and captions
        ax1 = fit_util.set_plt(ax1, fig_title=fig_title, grid=grid,
                               x_type=x_type, y_type=y_type, t_unit=t_unit,
                               logx=logx, logy=logy)

        # Plot table
        if table:
            # ax2 = plt.subplot2grid(shape=(10, 7), loc=(0, 1), rowspan=4, colspan=5)
            # ax2.axis('off')
            # columns = list(self.calibrate_result.__dict__['params'].valuesdict().keys())
            columns_to_show = [r'$L$ (m)', r'$\Delta t$ ($\rm{\mu}$s)']
            rows = ['Before', 'After']
            _row_before = [self.init_source_to_detector_m, self.init_offset_us]
            _row_after = [self.calibrated_source_to_detector_m, self.calibrated_offset_us]
            # for _each in columns:
            #     _row_after.append(self.calibrate_result.__dict__['params'].valuesdict()[_each])
            #     _row_before.append(self.params_to_calibrate.valuesdict()[_each])
            table = ax1.table(rowLabels=rows, colLabels=columns_to_show,  # colWidths=
                              cellText=[_row_before, _row_after],  # rows of data values
                              bbox=[0, -0.33, 1.0, 0.18]  # [left,bottom,width,height]
                              )
            # table.scale(0.5, 1)
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            plt.tight_layout()

        if save_fig:
            _sample_name = '_'.join(self.simulation.layer_list)
            _filename = 'calibration_' + _sample_name + '.png'
            plt.savefig(_filename, dpi=600, transparent=True)
            plt.close()
        return ax1

    def export(self, x_type='energy', y_type='attenuation', t_unit='us',
               index_level='iso', peak_id='indexed',
               before=False, interp=False, mixed=True):

        simu_label = 'ideal'
        exp_label = 'exp_raw'
        exp_before_label = 'exp_init'
        exp_interp_label = 'exp_interp'
        _df = pd.DataFrame()

        _col_suffix = fit_util.get_df_col_name(x_type=x_type)
        # Simulated total signal
        if mixed:
            _x = self.simulation.get_x(x_type=x_type,
                                       t_unit=t_unit,
                                       offset_us=self.calibrated_offset_us,
                                       source_to_detector_m=self.calibrated_source_to_detector_m,
                                       t_start_us=self.experiment.t_start_us,
                                       time_resolution_us=self.experiment.time_resolution_us)
            _y = self.simulation.get_y(y_type=y_type)
            _df['x_' + simu_label] = _x
            _df['y_' + simu_label] = _y

        """Plot options"""
        # Raw data before fitting
        if before:
            _x_init = self.experiment.get_x(x_type=x_type,
                                            t_unit=t_unit,
                                            offset_us=self.init_offset_us,
                                            source_to_detector_m=self.init_source_to_detector_m)
            _y_init = self.experiment.get_y(y_type=y_type,
                                            baseline=self.baseline)
            _df['x_' + exp_before_label] = _x_init
            _df['y_' + exp_before_label] = _y_init

        # 2.
        if interp:
            _exp_x_interp_calibrated, _exp_y_interp_calibrated = self.experiment.xy_scaled(
                x_type=x_type,
                y_type=y_type,
                energy_min=self.energy_min,
                energy_max=self.energy_max,
                energy_step=self.energy_step,
                t_unit=t_unit,
                offset_us=self.calibrated_offset_us,
                source_to_detector_m=self.calibrated_source_to_detector_m,
                baseline=self.baseline)
            # Interpolated raw data
            _df['x_' + exp_interp_label + _col_suffix] = _exp_x_interp_calibrated
            _df['y_' + exp_interp_label] = _exp_y_interp_calibrated
        else:
            # plot the calibrated raw data
            _x_cali = self.experiment.get_x(x_type=x_type,
                                            t_unit=t_unit,
                                            offset_us=self.calibrated_offset_us,
                                            source_to_detector_m=self.calibrated_source_to_detector_m)
            _y_cali = self.experiment.get_y(y_type=y_type,
                                            baseline=self.baseline)
            _df['x_' + exp_label + _col_suffix] = pd.Series(_x_cali)
            _df['y_' + exp_label] = pd.Series(_y_cali)

        # plot peaks detected and indexed
        if self.experiment.o_peak and self.experiment.o_peak.peak_map_indexed is not None:
            _peak_df_scaled = self.experiment.o_peak.peak_df_scaled
            _peak_map_indexed = self.experiment.o_peak.peak_map_indexed
            _peak_map_full = self.experiment.o_peak.peak_map_full
            _x_peak_exp_all = fit_util.convert_exp_peak_df(x_type=x_type, peak_df=_peak_df_scaled, t_unit=t_unit)
            _y_peak_exp_all = fit_util.convert_attenuation_to(y_type=y_type, y=_peak_df_scaled['y'])
            # _df = pd.concat([_df, _peak_df_scaled], axis=1)

            _df['x_peak_exp_all'] = pd.Series(_x_peak_exp_all)
            _df['y_peak_exp_all'] = pd.Series(_y_peak_exp_all)

            x_tag = fit_util.get_peak_tag(x_type=x_type)
            for _peak_name in _peak_map_indexed.keys():
                if len(_peak_map_full[_peak_name]['ideal']) > 0:
                    _x_peak_ideal_all = _peak_map_full[_peak_name]['ideal'][x_tag]
                    _y_peak_ideal_all = _peak_map_full[_peak_name]['ideal']['y']
                    _df['x_peak_ideal_all(' + _peak_name + ')'] = _x_peak_ideal_all
                    _df['y_peak_ideal_all(' + _peak_name + ')'] = _y_peak_ideal_all
                if len(_peak_map_indexed[_peak_name]['ideal']) > 0:
                    _x_peak_ideal_indexed = _peak_map_indexed[_peak_name]['ideal'][x_tag]
                    _y_peak_ideal_indexed = _peak_map_indexed[_peak_name]['ideal']['y']
                    _x_peak_exp_indexed = _peak_map_indexed[_peak_name]['exp'][x_tag]
                    _y_peak_exp_indexed = _peak_map_indexed[_peak_name]['exp']['y']
                    _df['x_peak_exp(' + _peak_name + ')'] = _x_peak_exp_indexed
                    _df['y_peak_exp(' + _peak_name + ')'] = _y_peak_exp_indexed
                    _df['x_peak_ideal(' + _peak_name + ')'] = _x_peak_ideal_indexed
                    _df['y_peak_ideal(' + _peak_name + ')'] = _y_peak_ideal_indexed

        _df.to_clipboard(index=False)

        return _df
    # def export_simu(self, filename=None, x_axis='energy', y_axis='attenuation',
    #                 all_layers=False, all_elements=False, all_isotopes=False, items_to_export=None,
    #                 t_start_us=1, time_resolution_us=0.16, time_unit='us'):
    #     if items_to_export is not None:
    #         # Shape items
    #         items = fit_util.Items(o_reso=self.simulation.o_reso, database=self.database)
    #         items_to_export = items.shaped(items_list=items_to_export)
    #
    #     self.simulation._export(filename=filename,
    #                             x_axis=x_axis,
    #                             y_axis=y_axis,
    #                             all_layers=all_layers,
    #                             all_elements=all_elements,
    #                             all_isotopes=all_isotopes,
    #                             items_to_export=items_to_export,
    #                             offset_us=self.calibrated_offset_us,
    #                             source_to_detector_m=self.calibrated_source_to_detector_m,
    #                             t_start_us=t_start_us,
    #                             time_resolution_us=time_resolution_us,
    #                             time_unit=time_unit)
