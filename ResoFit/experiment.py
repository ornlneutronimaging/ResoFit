import os

import ImagingReso._utilities as reso_util
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import ResoFit._utilities as fit_util
from ResoFit._utilities import load_txt_csv

x_type_list = ['energy', 'lambda', 'time']
y_type_list = ['attenuation', 'transmission']


class Experiment(object):
    def __init__(self, spectra_file, data_file, folder, repeat=1, baseline=False):
        """
        Load experiment data from 'YOUR_FILE_NAME.csv' or 'YOUR_FILE_NAME.txt' files
        :param spectra_file: data file stores the time-of-flight
        :param data_file: data file of neutron transmission
        :param repeat: input is needed only if the exp data is a summed result of multiple runs, default: 1, type: int
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
        if isinstance(repeat, int) is False:
            raise ValueError("Repeat value must be an integer!")
        elif repeat < 1:
            raise ValueError("Repeat value must be an integer >= 1 !")

        self.spectra = load_txt_csv(self.spectra_path)
        self.data = load_txt_csv(self.data_path)
        self.repeat = repeat
        assert type(self.repeat) is int

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

        # convert transmission into attenuation
        # self.data[0] = 1 - self.data[0]

        # raw image number saved
        self.img_num = self.data.index.values

    # def x_raw(self, x_type='energy', **kwargs):
    def x_raw(self, x_type='energy', offset_us=None, source_to_detector_m=None):
        """
        Get the 'x' in eV or angstrom with experimental parameters

        :param x_type:
        :type x_type:
        :param offset_us:
        :type offset_us:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :return:
        :rtype:
        """
        if x_type not in x_type_list:
            raise ValueError("'{}' is not supported. Must be one from {}.".format(x_type, x_type_list))
        # _kwarg_list = ['offset_us', 'source_to_detector_m']
        # for each_kwarg in list(kwargs.keys()):
        #     if each_kwarg not in _kwarg_list:
        #         raise ValueError("'{}' is not a valid **kwargs. Please refer '{}'".format(each_kwarg, _kwarg_list))
        # if 'offset_us' in kwargs.keys():
        #     self.offset_us = kwargs['offset_us']
        # if 'source_to_detector_m' in kwargs.keys():
        #     self.source_to_detector_m = kwargs['source_to_detector_m']
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        # For x_type == 'time' (x in seconds)
        x_exp_raw = np.array(self.spectra[0])
        if x_type == 'energy':
            x_exp_raw = np.array(reso_util.s_to_ev(array=x_exp_raw,
                                                   offset_us=self.offset_us,
                                                   source_to_detector_m=self.source_to_detector_m))
        elif x_type == 'lambda':
            x_exp_raw = np.array(reso_util.s_to_angstroms(array=x_exp_raw,
                                                          offset_us=self.offset_us,
                                                          source_to_detector_m=self.source_to_detector_m))
        # elif x_type == 'number':
        #     x_exp_raw = np.array(range(len(x_exp_raw)))

        return x_exp_raw

    def y_raw(self, y_type='attenuation', baseline=None, deg=7):
        """
        Get the 'y' in eV or angstrom with experimental parameters
        :param y_type: bool to switch between transmission and attenuation
        :param baseline: boolean to remove baseline/background by detrend
        :return: array
        """
        if y_type not in y_type_list:
            raise ValueError("'{}' is not supported. Must be one from {}.".format(y_type, y_type_list))
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline
        assert type(baseline) == bool

        y_exp_raw = np.array(self.data[0] / self.repeat)

        if y_type == 'attenuation':
            y_exp_raw = 1 - y_exp_raw
            if _baseline is True:
                y_exp_raw = fit_util.rm_baseline(y_exp_raw, deg=deg)
        else:  # y_type == 'transmission'
            if _baseline is True:
                y_exp_raw = fit_util.rm_envelope(y_exp_raw, deg=deg)

        return y_exp_raw

    def xy_scaled(self, energy_min, energy_max, energy_step,
                  x_type='energy', y_type='attenuation',
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
        :param baseline:
        :type baseline:
        :param offset_us:
        :type offset_us:
        :param source_to_detector_m:
        :type source_to_detector_m:
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

        x_exp_raw = self.x_raw(x_type=x_type,
                               offset_us=self.offset_us,
                               source_to_detector_m=self.source_to_detector_m)

        if x_type == 'energy':
            _x_max_energy = x_exp_raw[0]
            _x_min_energy = x_exp_raw[-1]
        elif x_type == 'lambda':
            _x_max_energy = reso_util.angstroms_to_ev(x_exp_raw[0])
            _x_min_energy = reso_util.ev_to_angstroms(x_exp_raw[-1])
        elif x_type == 'time':
            _x_max_energy = reso_util.s_to_ev(array=x_exp_raw[0],
                                              offset_us=self.offset_us,
                                              source_to_detector_m=self.source_to_detector_m)
            _x_min_energy = reso_util.s_to_ev(array=x_exp_raw[-1],
                                              offset_us=self.offset_us,
                                              source_to_detector_m=self.source_to_detector_m)
        else:
            raise ValueError("'{}' is not supported for scaling ".format(x_type))

        if energy_min < _x_min_energy:
            raise ValueError(
                "'Energy min' ({} eV) used for interpolation is beyond 'data min' ({} eV) ".format(energy_min,
                                                                                                   _x_min_energy))
        if energy_max > _x_max_energy:
            raise ValueError(
                "'Energy max' ({} eV) used for interpolation is beyond 'data max' ({} eV) ".format(energy_max,
                                                                                                   _x_max_energy))

        y_exp_raw = self.y_raw(y_type=y_type, baseline=_baseline, deg=deg)

        nbr_point = int((energy_max - energy_min) / energy_step + 1)
        x_interp = np.linspace(energy_min, energy_max, nbr_point)
        # y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='slinear')
        y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='cubic')
        y_interp = y_interp_function(x_interp)

        return x_interp, y_interp

    def slice(self, start=None, end=None, reset_index=False):
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
            self.data.drop(self.data.index[:start], inplace=True)
            self.spectra.drop(self.spectra.index[:start], inplace=True)
            self.slice_start = start
            # raw image number saved
            self.img_num = self.data.index.values
            if reset_index is True:
                self.spectra.reset_index(drop=True, inplace=True)
                self.data.reset_index(drop=True, inplace=True)

                # return self.spectra[0], self.data[0]

    def norm_to(self, file, reset_index=False):
        """
        Use specified file for normalization and save normalized data signal in self.data

        :param file: string. filename with suffix. ex: 'your_data.csv' inside the folder specified in __init__
        :param reset_index: True -> reset pd.Dataframe indexes after slicing
        :return: pd.Dataframe in place. normalized data signal in self.data
        """
        _full_path = os.path.join(self.folder_path, file)
        df = load_txt_csv(_full_path)
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
        # convert transmission into attenuation
        self.data[0] = self.data[0] / df[0]

    def find_peak(self, thres=0.15, min_dist=2, deg=7):
        """
        find and return x and y of detected peak in pd.DataFrame
        x is image number from data file. type (int)
        y is attenuation
        Note: impr_reso for finding peak is disable here to make sure the output image_num is integer
        :param thres:
        :type thres:
        :param min_dist:
        :type min_dist:

        :return:
        :rtype:
        """
        _y = self.data[0][:]
        _x = self.spectra[0][:]  # slicing is needed here to leave self.spectra[0] untouched

        _y = 1 - _y  # force to peaks
        _y = fit_util.rm_baseline(_y, deg=deg)  # force to remove baseline

        self.o_peak = fit_util.Peak()
        self.o_peak.find(_y,
                         x_name='x_num', y_name='y',
                         thres=thres, min_dist=min_dist, impr_reso=False)
        self.o_peak.add_x_col(x=_x, y=_y,
                              x_name='x_s', y_name='y',
                              thres=thres, min_dist=min_dist, impr_reso=False)

        if len(self.o_peak.peak_df) < 1:
            raise ValueError("No peak has been detected.")
        return self.o_peak.peak_df

    def scale_peak_with_ev(self, energy_min, energy_max,
                           calibrated_source_to_detector_m, calibrated_offset_us):
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
        self.o_peak.add_ev_and_scale(energy_min=energy_min, energy_max=energy_max,
                                     calibrated_source_to_detector_m=calibrated_source_to_detector_m,
                                     calibrated_offset_us=calibrated_offset_us)
        assert self.o_peak.peak_df_scaled is not None

        return self.o_peak.peak_df_scaled

    def plot(self, energy_xmax=150, lambda_xmax=None,
             y_type='transmission', baseline=None, deg=7,
             x_type='time', time_unit='us', offset_us=None, source_to_detector_m=None,
             logx=False, ax_mpl=None):
        """
        Display the loaded signal from data and spectra files.

        :param energy_xmax:
        :type energy_xmax:
        :param lambda_xmax:
        :type lambda_xmax:
        :param y_type:
        :type y_type:
        :param baseline:
        :type baseline:
        :param deg:
        :type deg:
        :param x_type:
        :type x_type:
        :param time_unit:
        :type time_unit:
        :param offset_us:
        :type offset_us:
        :param source_to_detector_m:
        :type source_to_detector_m:
        :param logx:
        :type logx:
        :param ax_mpl:
        :type ax_mpl:
        :return:
        :rtype:
        """
        _x_type_list = x_type_list[:]
        _x_type_list.append('number')
        _y_type_list = y_type_list
        _time_unit_list = ['s', 'us', 'ns']
        if x_type not in _x_type_list:
            raise ValueError("Please specify the x-axis type using one from '{}'.".format(_x_type_list))
        if y_type not in _y_type_list:
            raise ValueError("Please specify the y-axis type using one from '{}'.".format(_y_type_list))
        if time_unit not in _time_unit_list:
            raise ValueError("Please specify the time unit using one from '{}'.".format(_time_unit_list))
        if offset_us is not None:
            self.offset_us = offset_us
        if source_to_detector_m is not None:
            self.source_to_detector_m = source_to_detector_m
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline

        x_axis_label = None
        x_exp_raw = None
        df = pd.DataFrame()

        if ax_mpl is None:
            fig, ax_mpl = plt.subplots()
        """X-axis"""
        # determine values and labels for x-axis with options from
        # 'energy(eV)' & 'lambda(A)' & 'time(us)' & 'image number(#)'
        if x_type in ['energy', 'lambda']:
            if x_type == 'energy':
                x_axis_label = 'Energy (eV)'
                # angstrom = False
                ax_mpl.set_xlim(xmin=0, xmax=energy_xmax)
            else:
                x_axis_label = u"Wavelength (\u212B)"
                # angstrom = True
                if lambda_xmax is not None:
                    ax_mpl.set_xlim(xmin=0, xmax=lambda_xmax)
            x_exp_raw = self.x_raw(x_type=x_type,
                                   offset_us=self.offset_us,
                                   source_to_detector_m=self.source_to_detector_m)

        if x_type in ['time', 'number']:

            if x_type == 'time':
                if time_unit == 's':
                    x_axis_label = 'Time (s)'
                    x_exp_raw = self.spectra[0][:]
                if time_unit == 'us':
                    x_axis_label = 'Time (us)'
                    x_exp_raw = 1e6 * self.spectra[0][:]
                if time_unit == 'ns':
                    x_axis_label = 'Time (ns)'
                    x_exp_raw = 1e9 * self.spectra[0][:]

            if x_type == 'number':
                x_axis_label = 'Image number (#)'
                x_exp_raw = self.data.index.values

        assert x_axis_label is not None
        df[x_axis_label] = x_exp_raw

        """Y-axis"""
        # Determine to plot transmission or attenuation
        # Determine to put transmission or attenuation words for y-axis
        y_exp_raw = self.y_raw(y_type=y_type, baseline=_baseline, deg=deg)
        if y_type == 'transmission':
            y_axis_label = 'Neutron Transmission'
            ax_mpl.set_ylim(top=1.01 * max(y_exp_raw), bottom=-0.01)
        else:
            y_axis_label = 'Neutron Attenuation'
            ax_mpl.set_ylim(top=1.01, bottom=0.99 * min(y_exp_raw))

        assert y_axis_label is not None
        df[y_axis_label] = y_exp_raw

        # # Export
        # if filename is None:
        #     df.to_clipboard(excel=True)
        # else:
        #     df.to_csv(filename)

        # Plot
        if logx:
            ax_mpl.semilogx(x_exp_raw, y_exp_raw, '-o', label=self.data_file.split('.')[0], markersize=1)
        else:
            ax_mpl.plot(x_exp_raw, y_exp_raw, '-o', label=self.data_file.split('.')[0], markersize=1)
        ax_mpl.set_xlabel(x_axis_label)
        ax_mpl.set_ylabel(y_axis_label)
        ax_mpl.legend(loc='best')
        return ax_mpl
