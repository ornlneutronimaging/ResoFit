import numpy as np
import pandas as pd
import ImagingReso._utilities as reso_utils
import os
from scipy.interpolate import interp1d
import peakutils as pku
from ResoFit._utilities import load_txt_csv


class Experiment(object):
    source_to_detector_m = np.NaN
    offset_us = np.NaN
    spectra_path = ''
    data_path = ''
    spectra = pd.DataFrame
    data = pd.DataFrame
    slice_start = None
    slice_end = None

    def __init__(self, spectra_file, data_file, repeat=1, folder='data'):
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

        # Error for 'repeat' int & >=1
        if isinstance(repeat, int) is False:
            raise ValueError("Repeat value must be an integer!")
        if repeat < 1:
            raise ValueError("Repeat value must be an integer >= 1 !")

        self.spectra = load_txt_csv(self.spectra_path)
        self.data = load_txt_csv(self.data_path)
        self.repeat = repeat

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

    def x_raw(self, angstrom=False, offset_us=0., source_to_detector_m=15):
        """
        Get the 'x' in eV or angstrom with experimental parameters
        :param angstrom: bool to switch between eV and angstrom
        :param offset_us: offset_us for the actual measurement
        :param source_to_detector_m: detector position for the actual measurement
        :return: array
        """
        self.offset_us = offset_us
        self.source_to_detector_m = source_to_detector_m
        x_exp_raw = np.array(reso_utils.s_to_ev(self.spectra[0],  # x in seconds
                                                offset_us=offset_us,
                                                source_to_detector_m=source_to_detector_m))
        if angstrom is True:
            x_exp_raw = np.array(reso_utils.ev_to_angstroms(x_exp_raw))
        return x_exp_raw

    def y_raw(self, transmission=False):
        """
        Get the 'y' in eV or angstrom with experimental parameters
        :param transmission: bool to switch between transmission and attenuation
        :return: array
        """
        y_exp_raw = np.array(self.data[0]) / self.repeat
        if transmission is False:
            y_exp_raw = 1 - y_exp_raw
        # baseline = pku.baseline(y_exp_raw)
        # y_exp_raw = y_exp_raw - baseline
        return y_exp_raw

    def xy_scaled(self, energy_min, energy_max, energy_step, angstrom=False, transmission=False,
                  offset_us=0, source_to_detector_m=15):
        """
        Get interpolated x & y within the scaled range same as simulation
        :param energy_min:
        :param energy_max:
        :param energy_step:
        :param angstrom:
        :param transmission:
        :param offset_us:
        :param source_to_detector_m:
        :return:
        """
        self.offset_us = offset_us
        self.source_to_detector_m = source_to_detector_m
        x_exp_raw = reso_utils.s_to_ev(self.spectra[0],  # x in seconds
                                       offset_us=offset_us,
                                       source_to_detector_m=source_to_detector_m)
        _list = list(x_exp_raw)
        _x_max = _list[0]
        _x_min = _list[-1]
        if energy_min < _x_min:
            raise ValueError("'Energy min' ({} eV) used for interpolation is beyond 'data min' ({} eV) ".format(energy_min, _x_min))
        if energy_max > _x_max:
            raise ValueError("'Energy max' ({} eV) used for interpolation is beyond 'data max' ({} eV) ".format(energy_max, _x_max))

        y_exp_raw = np.array(self.data[0]) / self.repeat
        if transmission is False:
            y_exp_raw = 1 - y_exp_raw

        nbr_point = int((energy_max - energy_min) / energy_step + 1)
        x_interp = np.linspace(energy_min, energy_max, nbr_point)
        y_interp_function = interp1d(x=x_exp_raw, y=y_exp_raw, kind='cubic')
        y_interp = y_interp_function(x_interp)
        # baseline = pku.baseline(y_interp)
        # y_interp = y_interp - baseline

        if angstrom is True:
            x_interp = reso_utils.ev_to_angstroms(x_interp)
        return x_interp, y_interp

    def slice(self, slice_start=None, slice_end=None, reset_index=False):
        if slice_end is not None:
            if slice_end == slice_start:
                raise ValueError("The image number of 'start' ({}) and 'end' ({}) can not be the same.".format(slice_start, slice_end))
            self.data.drop(self.data.index[slice_end:], inplace=True)
            self.spectra.drop(self.spectra.index[slice_end:], inplace=True)
            # No 'index reset needed' after drop
            self.slice_end = slice_end
        if slice_start is not None:
            if slice_start == slice_end:
                raise ValueError("The image number of 'start' ({}) and 'end' ({}) can not be the same.".format(slice_start, slice_end))
            self.data.drop(self.data.index[:slice_start], inplace=True)
            self.spectra.drop(self.spectra.index[:slice_start], inplace=True)
            self.slice_start = slice_start
            if reset_index is True:
                self.spectra.reset_index(drop=True, inplace=True)
                self.data.reset_index(drop=True, inplace=True)
        return self.spectra[0], self.data[0]

    def norm_to(self, file, reset_index=False):
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
        self.data[0] = self.data[0] / df[0]

