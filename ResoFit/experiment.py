import os

import ImagingReso._utilities as reso_util
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import ResoFit._utilities as fit_util
from ResoFit._utilities import load_txt_csv


class Experiment(object):
    def __init__(self, spectra_file, data_file, folder, baseline=False):
        """
        Load experiment data from 'YOUR_FILE_NAME.csv' or 'YOUR_FILE_NAME.txt' files
        :param spectra_file: data file stores the time-of-flight
        :param data_file: data file of neutron transmission
        :param folder: folder name in str under /ResoFit directory
        """
        _file_path = os.path.abspath(os.path.dirname(__file__))
        self.folder_path = os.path.join(_file_path, folder)

        # Error for 'folder' existence
        if os.path.isdir(self.folder_path) is False:
            raise ValueError("Folder '{}' specified does not exist".format(folder))

        self.spectra_path = os.path.join(self.folder_path, spectra_file)
        self.data_path = os.path.join(self.folder_path, data_file)
        self.spectra_file = spectra_file
        self.data_file = data_file
        # Error for 'repeat' int & >=1
        # if isinstance(norm_factor, int) is False:
        #     raise ValueError("Repeat value must be an integer!")
        # elif norm_factor < 1:
        #     raise ValueError("Repeat value must be an integer >= 1 !")

        self.spectra = load_txt_csv(self.spectra_path)
        self.data = load_txt_csv(self.data_path)
        # self.norm_factor = norm_factor
        # self.data[0] = self.data[0] / norm_factor
        self.img_start = 0
        # assert type(self.norm_factor) is int

        # detector position (m) for the actual measurement
        self.source_to_detector_m = 16.
        # offset in time (us) for the actual measurement
        self.offset_us = 0.

        self.baseline = baseline
        self.slice_start = None
        self.slice_end = None
        self.o_peak = None

        # Error loading data and spectra
        if type(self.spectra[0][0]) is str:
            if self.spectra[0][0].islower() or self.spectra[0][0].isupper() is True:
                raise ValueError("Remove the axis descriptions in '{}' before loading ".format(spectra_file))
            else:
                raise ValueError("The file '{}' columns must be separated with 'tab' or ',' ".format(spectra_file))

        if type(self.data[0][0]) is str:
            if self.data[0][0].islower() or self.data[0][0].isupper() is True:
                raise ValueError("Remove the axis descriptions in '{}' before loading ".format(data_file))
            else:
                raise ValueError("The file '{}' columns must be separated with 'tab' or ',' ".format(data_file))
        if list(self.data[0][:4]) == [1, 2, 3, 4]:
            raise ValueError(
                "Duplicated index column was found in '{}', please remove duplicated column".format(data_file))
        # store raw data (df)
        self.data_raw = self.data
        self.spectra_raw = self.spectra
        self.t_start_us = fit_util.convert_s(self.spectra[0][0], t_unit='us')
        # self.t_start_us = 5.2
        self.time_resolution_us = fit_util.convert_s(self.spectra[0][2] - self.spectra[0][1], t_unit='us')
        # self.time_resolution_us = 0.160
        # convert transmission into attenuation
        # self.data[0] = 1 - self.data[0]

        # raw image number saved
        self.img_num = self.data.index.values

    def get_x(self, x_type='energy', offset_us=None, source_to_detector_m=None, t_unit='us'):
        """
        Get the 'x' in eV or angstrom with experimental parameters

        :param t_unit:
        :type t_unit:
        :param x_type:
        :type x_type:
        :param offset_us:
        :type offset_us:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :return:
        :rtype:
        """
        fit_util.check_if_in_list(x_type, fit_util.x_type_list)
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m

        _x_exp_raw = np.array(self.spectra[0])  # For x_type == 'time' (x in seconds)
        x_e = np.array(reso_util.s_to_ev(array=_x_exp_raw,
                                         offset_us=self.offset_us,
                                         source_to_detector_m=self.source_to_detector_m))
        x_exp_raw = fit_util.convert_energy_to(x=x_e, x_type=x_type, offset_us=self.offset_us,
                                               source_to_detector_m=self.source_to_detector_m, t_unit=t_unit,
                                               num_offset=self.img_start)
        return x_exp_raw

    def get_y(self, y_type='attenuation', baseline=None, deg=7):
        """
        Get the 'y' in eV or angstrom with experimental parameters
        :param deg:
        :type deg:
        :param y_type: bool to switch between transmission and attenuation
        :param baseline: boolean to remove baseline/background by detrend
        :return: array
        """
        fit_util.check_if_in_list(y_type, fit_util.y_type_list)
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline
        assert type(baseline) == bool
        # if norm_factor is None:
        #     _norm_factor = 1
        # else:
        #     _norm_factor = norm_factor

        y_exp_raw = np.array(self.data[0])

        if y_type == 'attenuation':
            y_exp_raw = 1 - y_exp_raw
            if _baseline is True:
                y_exp_raw = fit_util.rm_baseline(y_exp_raw, deg=deg)
        if y_type == 'transmission':
            if _baseline is True:
                y_exp_raw = fit_util.rm_envelope(y_exp_raw, deg=deg)

        return y_exp_raw

    def xy_scaled(self, energy_min, energy_max, energy_step,
                  x_type='energy', y_type='attenuation', t_unit='us',
                  offset_us=None, source_to_detector_m=None, baseline=None, deg=7):
        """
        Get interpolated x & y within the scaled range same as simulation

        :param energy_min:
        :type energy_min:
        :param energy_max:
        :type energy_max:
        :param energy_step:
        :type energy_step:
        :param x_type:
        :type x_type:
        :param y_type:
        :type y_type:
        :param t_unit:
        :type t_unit:
        :param offset_us:
        :type offset_us:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :param baseline:
        :type baseline:
        :param deg:
        :type deg:
        :return:
        :rtype:
        """
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline

        x_exp_raw = self.get_x(x_type='energy',
                               offset_us=self.offset_us,
                               source_to_detector_m=self.source_to_detector_m)

        _x_max_energy = x_exp_raw[0]
        _x_min_energy = x_exp_raw[-1]

        if energy_min < _x_min_energy:
            raise ValueError(
                "'Energy min' ({} eV) used for interpolation is beyond 'data min' ({} eV) ".format(energy_min,
                                                                                                   _x_min_energy))
        if energy_max > _x_max_energy:
            raise ValueError(
                "'Energy max' ({} eV) used for interpolation is beyond 'data max' ({} eV) ".format(energy_max,
                                                                                                   _x_max_energy))

        y_exp_raw = self.get_y(y_type=y_type, baseline=_baseline, deg=deg)

        nbr_point = int((energy_max - energy_min) / energy_step + 1)
        _x_interp = np.linspace(energy_min, energy_max, nbr_point)
        # y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='slinear')
        y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='cubic')
        y_interp = y_interp_function(_x_interp)
        x_interp = fit_util.convert_energy_to(x_type=x_type, x=_x_interp, t_unit=t_unit,
                                              offset_us=self.offset_us, source_to_detector_m=self.source_to_detector_m)
        return x_interp, y_interp

    def slice(self, start=None, end=None):
        """
        Slice the signal by image number

        :param start: start image
        :param end: end image
        :param reset_index: True -> reset pd.Dataframe indexes after slicing
        :return: pd.Dataframe. sliced self.spectra and self.data
        """
        if start and end is not None:
            if start > end:
                raise ValueError(
                    "The image number of 'start' ({}) can NOT be greater than 'end' ({}).".format(start,
                                                                                                  end))
            if end == start:
                raise ValueError(
                    "The image number of 'start' ({}) and 'end' ({}) can not be the same.".format(start,
                                                                                                  end))
        if end is not None:
            self.data.drop(self.data.index[end:], inplace=True)
            self.spectra.drop(self.spectra.index[end:], inplace=True)
            # No 'index reset needed' after drop
            self.slice_end = end
            # raw image number saved
            self.img_num = self.data.index.values
        if start is not None:
            self.img_start = start
            self.data.drop(self.data.index[:start], inplace=True)
            self.spectra.drop(self.spectra.index[:start], inplace=True)
            self.slice_start = start
            # raw image number saved
            self.img_num = self.data.index.values
            # Disabled reset_index #
            # if reset_index is True:
            #     self.spectra.reset_index(drop=True, inplace=True)
            #     self.data.reset_index(drop=True, inplace=True)
            #     self.img_start = 0

    def norm_to(self, file, norm_factor=1, reset_index=False):
        """
        Use specified file for normalization and save normalized data signal in self.data

        :param file: string. filename with suffix. ex: 'your_data.csv' inside the folder specified in __init__
        :param norm_factor:
        :type norm_factor:
        :param reset_index: True -> reset pd.Dataframe indexes after slicing
        :return: pd.DataFrame in place. normalized data signal in self.data
        """
        if file is not None:
            # Load file
            _full_path = os.path.join(self.folder_path, file)
            df = load_txt_csv(_full_path)
            # Resize length
            if len(self.data) != len(df):
                if self.slice_start is None and self.slice_end is None:
                    raise ValueError("The length of the 'norm_to_file' is not equal to the length of the data file.")
                else:
                    if self.slice_end is not None:
                        df.drop(df.index[self.slice_end:], inplace=True)
                    if self.slice_start is not None:
                        df.drop(df.index[:self.slice_start], inplace=True)
                        if reset_index is True:
                            df.reset_index(drop=True, inplace=True)
            self.data[0] = self.data[0] / df[0]
        # Apply norm_factor
        self.data[0] = self.data[0] / norm_factor

    def find_peak(self, thres=0.15, min_dist=2, baseline=None, deg=7):
        """
        find and return x and y of detected peak in pd.DataFrame
        x is image number from data file. type (int)
        y is attenuation
        Note: impr_reso for finding peak is disabled here to make sure the output image_num is integer
        :param baseline:
        :type baseline:
        :param thres:
        :type thres:
        :param min_dist:
        :type min_dist:
        :param deg:
        :type deg:

        :return:
        :rtype:
        """
        _y = self.data[0][:]
        _x = self.spectra[0][:]  # slicing is needed here to leave self.spectra[0] untouched
        _y = 1 - _y  # force to peaks
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline
        if _baseline:
            _y = fit_util.rm_baseline(_y, deg=deg)  # force to remove baseline
        self.o_peak = fit_util.Peak()
        self.o_peak.find(x=_x, y=_y, y_name='y', thres=thres, min_dist=min_dist, impr_reso=False)
        if len(self.o_peak.peak_df) < 1:
            raise ValueError("No peak has been detected.")
        return self.o_peak.peak_df

    def _scale_peak_with_ev(self, energy_min, energy_max, offset_us=None, source_to_detector_m=None):
        """

        :param energy_min:
        :type energy_min:
        :param energy_max:
        :type energy_max:
        :param calibrated_source_to_detector_m:
        :type calibrated_source_to_detector_m:
        :param calibrated_offset_us:
        :type calibrated_offset_us:
        :return:
        :rtype:
        """
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        self.o_peak._extend_x_cols(offset_us=self.offset_us,
                                   source_to_detector_m=self.source_to_detector_m)
        self.o_peak._scale_peak_df(energy_min=energy_min, energy_max=energy_max)
        assert self.o_peak.peak_df_scaled is not None
        return self.o_peak.peak_df_scaled

    def plot(self,
             y_type='transmission', baseline=None, deg=7,
             x_type='time', t_unit='us', offset_us=None, source_to_detector_m=None,
             logx=False, logy=False, ax_mpl=None, fmt='.', ms=2, lw=1.5, alpha=1, grid=False):
        """
        Display the loaded signal from data and spectra files.
        """
        fit_util.check_if_in_list(x_type, fit_util.x_type_list)
        fit_util.check_if_in_list(y_type, fit_util.y_type_list)
        fit_util.check_if_in_list(t_unit, fit_util.t_unit_list)
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline
        fig_title = 'Experimental data'

        if ax_mpl is None:
            fig, ax_mpl = plt.subplots()
        """X-axis"""
        x_exp_raw = self.get_x(x_type=x_type,
                               t_unit=t_unit,
                               offset_us=self.offset_us,
                               source_to_detector_m=self.source_to_detector_m)
        """Y-axis"""
        y_exp_raw = self.get_y(y_type=y_type, baseline=_baseline, deg=deg)
        # if y_type == 'transmission':
        #     y_axis_label = 'Neutron Transmission'
        #     ax_mpl.set_ylim(top=1.01 * max(y_exp_raw), bottom=-0.01)
        # else:
        #     y_axis_label = 'Neutron Attenuation'
        #     ax_mpl.set_ylim(top=1.01, bottom=0.99 * min(y_exp_raw))

        # # Export
        # if filename is None:
        #     df.to_clipboard(excel=True)
        # else:
        #     df.to_csv(filename)

        # Plot
        ax_mpl.plot(x_exp_raw, y_exp_raw, fmt, label=self.data_file.split('.')[0] + '_data',
                    ms=ms, lw=lw, alpha=alpha)
        if self.o_peak is not None:
            if self.o_peak.peak_df_scaled is not None:
                _current_peak_df = self.o_peak.peak_df_scaled
            elif self.o_peak.peak_df is not None:
                _current_peak_df = self.o_peak.peak_df
            else:
                _current_peak_df = None
            if _current_peak_df is not None:
                _x_tag = fit_util.get_peak_tag(x_type=x_type)
                ax_mpl.scatter(_current_peak_df[_x_tag],
                               _current_peak_df['y'],
                               c='k',
                               marker='x',
                               # s=30,
                               # marker='o',
                               # facecolors='none',
                               # edgecolors='k',
                               label='_nolegend_')
        ax_mpl = fit_util.set_plt(ax=ax_mpl, fig_title=fig_title, grid=grid,
                                  x_type=x_type, y_type=y_type, t_unit=t_unit, logx=logx, logy=logy)
        return ax_mpl
