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

    def __init__(self, spectra_file, data_file, repeat=1, folder='data'):
        """
        Load experiment data from 'YOUR_FILE_NAME.csv' or 'YOUR_FILE_NAME.txt' files
        :param spectra_file: data file stores the time-of-flight
        :param data_file: data file of neutron transmission
        :param repeat: input is needed only if the exp data is a summed result of multiple runs
        :param folder: folder name in str under /ResoFit directory

        """
        _file_path = os.path.abspath(os.path.dirname(__file__))
        _folder_path = os.path.join(_file_path, folder)

        # Error for 'folder' existence
        if os.path.isdir(_folder_path) is False:
            raise ValueError("Folder '{}' specified does not exist".format(folder))
        # Error for 'spectra_file' format and existence
        spectra_format = spectra_file[-4:]
        if spectra_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        self.spectra_path = os.path.join(_folder_path, spectra_file)
        if os.path.exists(self.spectra_path) is False:
            raise ValueError("Can not find spectra file '{}' in '{}' folder".format(spectra_file, _folder_path))
        # Error for 'data_file' format and existence
        date_format = data_file[-4:]
        if date_format not in ['.txt', '.csv']:
            raise ValueError("Spectra file must be in the format of '.txt' or '.csv'")
        self.data_path = os.path.join(_folder_path, data_file)
        if os.path.exists(self.data_path) is False:
            raise ValueError("Can not find data file '{}' in '{}' folder".format(data_file, _folder_path))
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
        if list(self.data[0][:4]) == [1, 2, 3, 4]:
            y_exp_raw = np.array(self.data[1]) / self.repeat
        else:
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
        if list(self.data[0][:4]) == [1, 2, 3, 4]:
            y_exp_raw = np.array(self.data[1]) / self.repeat
        else:
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

    def xy_sliced(self, slice_start, slice_end, baseline=True):
        x_exp_raw_sliced = self.spectra[0][slice_start:slice_end]
        if list(self.data[0][:4]) == [1, 2, 3, 4]:
            y_exp_raw_sliced = np.array(self.data[1]) / self.repeat
        else:
            y_exp_raw_sliced = np.array(self.data[0]) / self.repeat

        y_exp_raw_sliced = 1 - y_exp_raw_sliced[slice_start:slice_end]#/2.574063196
        if baseline is True:
            baseline = pku.baseline(y_exp_raw_sliced)
            y_exp_raw_sliced = y_exp_raw_sliced - baseline

        return x_exp_raw_sliced, y_exp_raw_sliced
