import os

import ImagingReso._utilities as reso_util
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import ResoFit._utilities as fit_util
from ResoFit._utilities import load_txt_csv


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
        # self.calibrated_offset_us = None
        # self.calibrated_source_to_detector_m = None

        self.baseline = baseline
        self.slice_start = None
        self.slice_end = None
        self.peak_df_raw = None

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

    def x_raw(self, x_type='energy', **kwargs):
        """
        Get the 'x' in eV or angstrom with experimental parameters
        :param x_type: bool to switch between eV and angstrom
        :return: array
        """
        if x_type not in ['energy', 'lambda']:
            raise ValueError("'{}' is not supported. Must be one from ['energy', 'lambda'].")
        if 'offset_us' in kwargs.keys():
            self.offset_us = kwargs['offset_us']
        if 'source_to_detector_m' in kwargs.keys():
            self.source_to_detector_m = kwargs['source_to_detector_m']

        x_exp_raw = np.array(reso_util.s_to_ev(array=self.spectra[0],  # x in seconds
                                               offset_us=self.offset_us,
                                               source_to_detector_m=self.source_to_detector_m))
        if x_type == 'lambda':
            x_exp_raw = np.array(reso_util.ev_to_angstroms(x_exp_raw))
        return x_exp_raw

    def y_raw(self, y_type='attenuation', baseline=None):
        """
        Get the 'y' in eV or angstrom with experimental parameters
        :param y_type: bool to switch between transmission and attenuation
        :param baseline: boolean to remove baseline/background by detrend
        :return: array
        """
        if y_type not in ['attenuation', 'transmission']:
            raise ValueError("'{}' is not supported. Must be one from ['attenuation', 'transmission'].")
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline
        assert type(baseline) == bool

        y_exp_raw = np.array(self.data[0]) / self.repeat

        if y_type == 'attenuation':
            y_exp_raw = 1 - y_exp_raw
            if _baseline is True:  # baseline removal only works for peaks instead of dips currently
                y_exp_raw = fit_util.rm_baseline(y_exp_raw)
        elif y_type == 'transmission':
            if _baseline is True:  # baseline removal only works for peaks instead of dips currently
                raise ValueError("Baseline removal only works for peaks instead of dips!")

        return y_exp_raw

    def xy_scaled(self, energy_min, energy_max, energy_step,
                  x_type='energy', y_type='attenuation', baseline=None, **kwargs):
        """
        Get interpolated x & y within the scaled range same as simulation
        :param baseline: boolean to remove baseline/background by detrend
        :param energy_min:
        :param energy_max:
        :param energy_step:
        :param x_type:
        :param y_type:
        :return: np.array. interpolated x_exp (in eV or angstrom) and y_exp with specified energy range and step
        """
        if 'offset_us' in kwargs.keys():
            self.offset_us = kwargs['offset_us']
        if 'source_to_detector_m' in kwargs.keys():
            self.source_to_detector_m = kwargs['source_to_detector_m']
        if baseline is None:
            _baseline = self.baseline
        else:
            _baseline = baseline
        if x_type not in ['energy', 'lambda']:
            raise ValueError("'{}' is not supported. Must be one from ['energy', 'lambda'].")
        if y_type not in ['attenuation', 'transmission']:
            raise ValueError("'{}' is not supported. Must be one from ['attenuation', 'transmission'].")
        assert type(_baseline) == bool

        x_exp_raw = reso_util.s_to_ev(array=self.spectra[0],  # x in seconds
                                      offset_us=self.offset_us,
                                      source_to_detector_m=self.source_to_detector_m)
        _list = list(x_exp_raw)
        _x_max = _list[0]
        _x_min = _list[-1]
        if energy_min < _x_min:
            raise ValueError(
                "'Energy min' ({} eV) used for interpolation is beyond 'data min' ({} eV) ".format(energy_min, _x_min))
        if energy_max > _x_max:
            raise ValueError(
                "'Energy max' ({} eV) used for interpolation is beyond 'data max' ({} eV) ".format(energy_max, _x_max))

        y_exp_raw = np.array(self.data[0]) / self.repeat
        if y_type == 'attenuation':
            y_exp_raw = 1 - y_exp_raw

        nbr_point = int((energy_max - energy_min) / energy_step + 1)
        x_interp = np.linspace(energy_min, energy_max, nbr_point)
        y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='slinear')
        y_interp = y_interp_function(x_interp)

        if y_type == 'attenuation':
            if _baseline is True:  # baseline removal only works for peaks instead of dips currently
                y_interp = fit_util.rm_baseline(y_interp)
        elif y_type == 'transmission':
            if _baseline is True:
                raise ValueError("Baseline removal only works for peaks instead of dips!")

        if x_type == 'lambda':
            x_interp = reso_util.ev_to_angstroms(x_interp)
        return x_interp, y_interp

    def slice(self, slice_start=None, slice_end=None, reset_index=False):
        """
        Slice the signal by image number
        :param slice_start: start image
        :param slice_end: end image
        :param reset_index: True -> reset pd.Dataframe indexes after slicing
        :return: pd.Dataframe. sliced self.spectra and self.data
        """
        if slice_start and slice_end is not None:
            if slice_start > slice_end:
                raise ValueError(
                    "The image number of 'start' ({}) can NOT be greater than 'end' ({}).".format(slice_start,
                                                                                                  slice_end))
            if slice_end == slice_start:
                raise ValueError(
                    "The image number of 'start' ({}) and 'end' ({}) can not be the same.".format(slice_start,
                                                                                                  slice_end))
        if slice_end is not None:
            self.data.drop(self.data.index[slice_end:], inplace=True)
            self.spectra.drop(self.spectra.index[slice_end:], inplace=True)
            # No 'index reset needed' after drop
            self.slice_end = slice_end
            # raw image number saved
            self.img_num = self.data.index.values
        if slice_start is not None:
            self.data.drop(self.data.index[:slice_start], inplace=True)
            self.spectra.drop(self.spectra.index[:slice_start], inplace=True)
            self.slice_start = slice_start
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

    def find_peak(self, thres=0.15, min_dist=2):
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
        # remove baseline
        _y = fit_util.rm_baseline(self.data[0])
        _x = self.spectra[0]
        _index_gap = 0
        if self.baseline is False:
            _y = fit_util.rm_baseline(self.data[0])
        if self.slice_start is not None:
            _y.reset_index(drop=True, inplace=True)
            _x.reset_index(drop=True, inplace=True)
            _index_gap = self.slice_start
        _peak_num_df = fit_util.find_peak(x=_y.index.values + _index_gap, y=_y,
                                          thres=thres, min_dist=min_dist, impr_reso=False)
        _peak_time_df = fit_util.find_peak(x=_x, y=_y,
                                           thres=thres, min_dist=min_dist, impr_reso=False)
        # assert _peak_num_df['y'] == _peak_time_df['y']
        # _peak_df.drop(_peak_df[_peak_df.x < self.energy_min].index, inplace=True)
        # _peak_df.drop(_peak_df[_peak_df.x > self.energy_max].index, inplace=True)
        # _peak_df.reset_index(drop=True, inplace=True)
        _peak_df = pd.DataFrame()
        _peak_df['y'] = _peak_num_df['y']
        _peak_df['x_num'] = _peak_num_df['x']
        _peak_df['x_s'] = _peak_time_df['x']
        if len(_peak_df['y']) < 1:
            raise ValueError("No peak has been detected.")
        self.peak_df_raw = _peak_df

        return self.peak_df_raw

    def plot_raw(self, energy_xmax=150, lambda_xmax=None,
                 transmission=False, baseline=False,
                 x_axis='energy', time_unit='us', **kwargs):
        """
        Display the loaded signal from data and spectra files.
        :param energy_xmax: maximum x-axis energy value to display
        :param lambda_xmax: maximum x-axis lambda value to display
        :param transmission: boolean. False -> show resonance peaks
                                      True -> show resonance dips
        :param baseline: boolean. True -> remove baseline by detrend
        :param x_axis: string. x-axis type, must be either 'energy' or 'lambda' or 'time' or 'number'
        :param time_unit: string. Must be either 's' or 'us' or 'ns'
        :return: display raw data signals
        """
        if x_axis not in ['energy', 'lambda', 'time', 'number']:
            raise ValueError("Please specify the x-axis type using one from '['energy', 'lambda', 'time', 'number']'.")
        if time_unit not in ['s', 'us', 'ns']:
            raise ValueError("Please specify the time unit using one from '['s', 'us', 'ns']'.")
        x_axis_label = None
        x_exp_raw = None
        if 'offset_us' in kwargs.keys():
            self.offset_us = kwargs['offset_us']
        if 'source_to_detector_m' in kwargs.keys():
            self.source_to_detector_m = kwargs['source_to_detector_m']
        if baseline is True:
            if self.baseline is True:
                baseline = False

        """X-axis"""
        # determine values and labels for x-axis with options from
        # 'energy(eV)' & 'lambda(A)' & 'time(us)' & 'image number(#)'
        if x_axis in ['energy', 'lambda']:
            if x_axis == 'energy':
                x_axis_label = 'Energy (eV)'
                angstrom = False
                plt.xlim(xmin=0, xmax=energy_xmax)
            else:
                x_axis_label = u"Wavelength (\u212B)"
                angstrom = True
                if lambda_xmax is not None:
                    plt.xlim(xmin=0, xmax=lambda_xmax)
            x_exp_raw = self.x_raw(angstrom=angstrom,
                                   offset_us=self.offset_us,
                                   source_to_detector_m=self.source_to_detector_m)

        if x_axis in ['time', 'number']:
            x_exp_raw = self.x_raw(angstrom=False,
                                   offset_us=self.offset_us,
                                   source_to_detector_m=self.source_to_detector_m)
            if x_axis == 'time':
                if time_unit == 's':
                    x_axis_label = 'Time (s)'
                    x_exp_raw = self.spectra[0]
                if time_unit == 'us':
                    x_axis_label = 'Time (us)'
                    x_exp_raw = 1e6 * self.spectra[0]
                if time_unit == 'ns':
                    x_axis_label = 'Time (ns)'
                    x_exp_raw = 1e9 * self.spectra[0]

            if x_axis == 'number':
                x_axis_label = 'Image number (#)'
                x_exp_raw = self.data.index.values
        if x_axis_label is None:
            raise ValueError("x_axis_label does NOT exist, please check.")

        """Y-axis"""
        # Determine to plot transmission or attenuation
        # Determine to put transmission or attenuation words for y-axis
        if transmission:
            y_axis_label = 'Neutron Transmission'
        else:
            y_axis_label = 'Neutron Attenuation'
        y_exp_raw = self.y_raw(transmission=transmission, baseline=baseline)

        # Plot
        plt.plot(x_exp_raw, y_exp_raw, 'o', label=self.data_file, markersize=2)
        plt.ylim(ymax=1.01)
        plt.xlabel(x_axis_label)
        plt.ylabel(y_axis_label)
        plt.legend(loc='best')

    def export_raw(self, filename=None,
                   transmission=False, baseline=False,
                   x_axis='energy', time_unit='us', **kwargs):
        """
        Export the calculated signal from data and spectra files.
        :param filename: filename (with .csv suffix) you would like to save as
                                None -> export to clipboard
        :type filename: string.
        :param transmission: boolean. False -> show resonance peaks
                                      True -> show resonance dips
        :param baseline: boolean. True -> remove baseline by detrend
        :param x_axis: string. x-axis type, must be either 'energy' or 'lambda' or 'time' or 'number'
        :param time_unit: string. Must be either 's' or 'us' or 'ns'
        :return: display raw data signals
        """
        if x_axis not in ['energy', 'lambda', 'time', 'number']:
            raise ValueError("Please specify the x-axis type using one from '['energy', 'lambda', 'time', 'number']'.")
        if time_unit not in ['s', 'us', 'ns']:
            raise ValueError("Please specify the time unit using one from '['s', 'us', 'ns']'.")
        x_axis_label = None
        x_exp_raw = None
        df = pd.DataFrame()
        if 'offset_us' in kwargs.keys():
            self.offset_us = kwargs['offset_us']
        if 'source_to_detector_m' in kwargs.keys():
            self.source_to_detector_m = kwargs['source_to_detector_m']
        if baseline is True:
            if self.baseline is True:
                baseline = False

        """X-axis"""
        # determine values and labels for x-axis with options from
        # 'energy(eV)' & 'lambda(A)' & 'time(us)' & 'image number(#)'
        if x_axis in ['energy', 'lambda']:
            if x_axis == 'energy':
                x_axis_label = 'Energy (eV)'
                angstrom = False
            else:
                x_axis_label = u"Wavelength (\u212B)"
                angstrom = True
            x_exp_raw = self.x_raw(angstrom=angstrom, offset_us=self.offset_us,
                                   source_to_detector_m=self.source_to_detector_m)

        if x_axis in ['time', 'number']:
            x_exp_raw = self.x_raw(angstrom=False, offset_us=self.offset_us,
                                   source_to_detector_m=self.source_to_detector_m)
            if x_axis == 'time':
                if time_unit == 's':
                    x_axis_label = 'Time (s)'
                    x_exp_raw = self.spectra[0]
                if time_unit == 'us':
                    x_axis_label = 'Time (us)'
                    x_exp_raw = 1e6 * self.spectra[0]
                if time_unit == 'ns':
                    x_axis_label = 'Time (ns)'
                    x_exp_raw = 1e9 * self.spectra[0]

            if x_axis == 'number':
                x_axis_label = 'Image number (#)'
                x_exp_raw = np.array(range(1, len(self.data[0]) + 1))
        if x_axis_label is None:
            raise ValueError("x_axis_label does NOT exist, please check.")

        df[x_axis_label] = x_exp_raw

        """Y-axis"""
        # Determine to plot transmission or attenuation
        # Determine to put transmission or attenuation words for y-axis
        if transmission:
            y_axis_label = 'Neutron Transmission'
        else:
            y_axis_label = 'Neutron Attenuation'
        y_exp_raw = self.y_raw(transmission=transmission, baseline=baseline)
        df[y_axis_label] = y_exp_raw

        # Export
        if filename is None:
            df.to_clipboard(excel=True)
        else:
            df.to_csv(filename)
