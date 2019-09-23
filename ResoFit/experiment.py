import os

import ImagingReso._utilities as reso_util
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import ResoFit._utilities as fit_util
from ResoFit._utilities import load_txt_csv


class Experiment(object):
    def __init__(self, spectra_file: str, data_file: str, folder: str,
                 source_to_detector_m, offset_us):
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

        self.spectra_file = spectra_file
        self.data_file = data_file
        self.source_to_detector_m = source_to_detector_m
        self.offset_us = offset_us
        self.spectra_path = os.path.join(self.folder_path, spectra_file)
        self.data_path = os.path.join(self.folder_path, data_file)

        # Load spectrum and data
        self.spectra = load_txt_csv(self.spectra_path)
        self.data = load_txt_csv(self.data_path)
        self.t_unit = 'us'
        self.img_start = 0

        # Default slice parameter
        self.slice_start = None
        self.slice_end = None

        # Class to store peak info
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
        self.data_raw = self.data[:]
        self.spectra_raw = self.spectra[:]
        self.t_start_us = fit_util.convert_s(self.spectra[0][0], t_unit='us')
        self.time_resolution_us = fit_util.convert_s(self.spectra[0][2] - self.spectra[0][1], t_unit='us')
        # raw image number saved
        self.img_num = self.data.index.values

    def get_x(self, x_type, offset_us=None, source_to_detector_m=None, t_unit='us'):
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
        if t_unit != self.t_unit:
            self.t_unit = t_unit

        _x_exp_raw = np.array(self.spectra[0][:])  # Default x_type == 'time' (x in seconds)
        x_e = np.array(reso_util.s_to_ev(array=_x_exp_raw,
                                         offset_us=self.offset_us,
                                         source_to_detector_m=self.source_to_detector_m))
        x_exp_raw = fit_util.convert_energy_to(x=x_e, x_type=x_type, offset_us=self.offset_us,
                                               source_to_detector_m=self.source_to_detector_m,
                                               t_unit=self.t_unit,
                                               num_offset=self.img_start)
        return x_exp_raw

    def get_y(self, y_type, baseline: bool, baseline_deg=None):
        """
        Get the 'y' in eV or angstrom with experimental parameters
        :param baseline_deg:
        :type baseline_deg:
        :param y_type: bool to switch between transmission and attenuation
        :param baseline: boolean to remove baseline/background by detrend
        :return: array
        """
        fit_util.check_if_in_list(y_type, fit_util.y_type_list)
        y_exp_raw = np.array(self.data[0][:])
        if baseline:
            if type(baseline_deg) is not int:
                raise ValueError("{} 'baseline_deg' only accept 'int' input".format(baseline_deg))
            else:
                y_exp_raw = fit_util.rm_envelope(y_exp_raw, deg=baseline_deg)

        if y_type == 'attenuation':
            y_exp_raw = 1 - y_exp_raw

        return y_exp_raw

    def xy_scaled(self, energy_min, energy_max, energy_step,
                  x_type, y_type, baseline, baseline_deg=None, t_unit='us',
                  offset_us=None, source_to_detector_m=None):
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
        :param baseline_deg:
        :type baseline_deg:
        :return:
        :rtype:
        """
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        if t_unit != self.t_unit:
            self.t_unit = t_unit

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

        y_exp_raw = self.get_y(y_type=y_type, baseline=baseline, baseline_deg=baseline_deg)

        nbr_point = int((energy_max - energy_min) / energy_step + 1)
        _x_interp = np.linspace(energy_min, energy_max, nbr_point)
        # y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='slinear')
        y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='cubic')
        y_interp = y_interp_function(_x_interp)
        x_interp = fit_util.convert_energy_to(x_type=x_type, x=_x_interp, t_unit=self.t_unit,
                                              offset_us=self.offset_us, source_to_detector_m=self.source_to_detector_m)
        return x_interp, y_interp

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

    def find_peak(self, x_type, y_type, thres, min_dist, baseline, baseline_deg=None, imprv_reso=False):
        """
        find and return x and y of detected peak in pd.DataFrame
        x is image number from data file. type (int)
        y is attenuation
        Note: impr_reso for finding peak is disabled here to make sure the output image_num is integer

        :param x_type:
        :type x_type:
        :param y_type:
        :type y_type:
        :param thres:
        :type thres:
        :param min_dist:
        :type min_dist:
        :param baseline:
        :type baseline:
        :param baseline_deg:
        :type baseline_deg:
        :param imprv_reso:
        :type imprv_reso:
        :return:
        :rtype:
        """

        _x = self.get_x(
            x_type=x_type,
            offset_us=self.offset_us,
            source_to_detector_m=self.source_to_detector_m,
            t_unit=self.t_unit
        )
        _y = self.get_y(
            y_type='attenuation',
            baseline=baseline,
            baseline_deg=baseline_deg
        )
        self.o_peak = fit_util.ResoPeak(x=_x, y=_y, x_type=x_type, y_type=y_type)
        self.o_peak.find_peak(thres=thres, min_dist=min_dist, imprv_reso=imprv_reso)
        if len(self.o_peak.peak_dict['y']) < 1:
            raise ValueError("No peak has been detected.")
        if y_type == 'transmission':
            self.o_peak.peak_dict['y'] = 1 - self.o_peak.peak_dict['y']
        return self.o_peak.peak_dict

    def plot(self, x_type, y_type, baseline, baseline_deg=None,
             t_unit='us', offset_us=None, source_to_detector_m=None,
             logx=False, logy=False, ax_mpl=None, fmt='.', ms=2, lw=1.5, alpha=1,
             grid=False, label=None, plot_before=False):
        """
        Display the loaded signal from data and spectra files.
        """
        self.__check_in_list(x_type=x_type, y_type=y_type, t_unit=t_unit)
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        if t_unit != self.t_unit:
            self.t_unit = t_unit

        fig_title = 'Experimental data'
        if label is None:
            _label = self.data_file.split('.')[0] + '_data'
        else:
            _label = label

        if ax_mpl is None:
            fig, ax_mpl = plt.subplots()
        """X-axis"""
        x_exp_raw = self.get_x(x_type=x_type,
                               t_unit=self.t_unit,
                               offset_us=self.offset_us,
                               source_to_detector_m=self.source_to_detector_m)
        """Y-axis"""
        if baseline:
            if plot_before:
                y_exp_raw_before = self.get_y(y_type=y_type, baseline=False, baseline_deg=baseline_deg)
                ax_mpl.plot(x_exp_raw, y_exp_raw_before, '--', label='Before baseline removal', ms=ms, lw=lw,
                            alpha=alpha)
        y_exp_raw = self.get_y(y_type=y_type, baseline=baseline, baseline_deg=baseline_deg)

        # Plot
        assert y_exp_raw.shape == x_exp_raw.shape
        if len(y_exp_raw) - len(x_exp_raw) == 1:
            y_exp_raw = y_exp_raw[:-1]
        ax_mpl.plot(x_exp_raw, y_exp_raw, fmt, label=_label, ms=ms, lw=lw, alpha=alpha)

        if self.o_peak is not None:
            if len(self.o_peak.peak_dict) != 0:
                _x_tag = fit_util.get_peak_tag(x_type=x_type)
                ax_mpl.scatter(self.o_peak.peak_dict['x'],
                               self.o_peak.peak_dict['y'],
                               c='r',
                               marker='x',
                               # s=30,
                               # marker='o',
                               # facecolors='none',
                               # edgecolors='k',
                               label='_nolegend_')
        ax_mpl = fit_util.set_plt(ax=ax_mpl, fig_title=fig_title, grid=grid,
                                  x_type=x_type, y_type=y_type, t_unit=t_unit, logx=logx, logy=logy)
        return ax_mpl

    def export(self, x_type, y_type, baseline, baseline_deg=None,
               t_unit='us', offset_us=None, source_to_detector_m=None):
        self.__check_in_list(x_type=x_type, y_type=y_type, t_unit=t_unit)
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        if t_unit != self.t_unit:
            self.t_unit = t_unit
        _df = pd.DataFrame()
        """X-axis"""
        x_exp_raw = self.get_x(x_type=x_type,
                               t_unit=t_unit,
                               offset_us=self.offset_us,
                               source_to_detector_m=self.source_to_detector_m)
        """Y-axis"""
        y_exp_raw = self.get_y(y_type=y_type, baseline=baseline, baseline_deg=baseline_deg)

        _df['x'] = x_exp_raw
        _df['y'] = y_exp_raw
        _df.to_clipboard(index=False)

        return _df

    def __check_in_list(self, x_type, y_type, t_unit):
        fit_util.check_if_in_list(x_type, fit_util.x_type_list)
        fit_util.check_if_in_list(y_type, fit_util.y_type_list)
        fit_util.check_if_in_list(t_unit, fit_util.t_unit_list)
